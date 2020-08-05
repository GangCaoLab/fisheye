import os
from os.path import join, exists
import random
import fire
import numpy as np
import math
from numpy import array
import pandas as pd
from pyfaidx import Fasta
import RNA
import primer3
import glob
import typing as t
import subprocess as subp
import pysam
import editdistance
from pathos.pools import ProcessPool
from fisheye.utils import get_logger
from fisheye.utils import reverse_complement
from fisheye.primer_design.coding import coding_llhc, coding_random


log = get_logger(__name__)


def get_tmp_dir(basename):
    i = 0
    dirname = lambda: f"{basename}.{i}"
    while exists(dirname()):
        i += 1
    os.mkdir(dirname())
    return dirname()


TMP_DIR = None


def read_gene(genelist):
    genelist = pd.read_csv(genelist, header=None, comment='#')
    if genelist.shape[1] == 1:
        genelist.columns = ['geneID']
    else:
        genelist = genelist.iloc[:, :2]
        genelist.columns = ['geneID', 'score']
    return genelist


def read_gtf(gtf):
    df = pd.read_csv(gtf, sep='\t', header=None, comment='#')
    df.columns = ['chr', 'source', 'type', 'start', 'end', 'X1', 'strand', 'X2', 'fields']
    df['gene_name'] = df.fields.str.extract("gene_name \"(.*?)\"")
    df['length'] = df['end'] - df['start']
    return df


Exon = t.Tuple[str, str, int]  # (name, seq, n_trans)


def extract_exons(df_gtf: pd.DataFrame, fa: Fasta,
                  genelist: pd.DataFrame, min_length: int=40
                  ) -> t.Mapping[str, t.List[Exon]]:
    """Extract all exons of each gene"""
    gene2exons = {}
    for item in genelist.iterrows():
        gene = item[1]['geneID']
        df_gene = df_gtf[df_gtf.gene_name == gene]
        if df_gene.shape[0] <= 0:
            raise ValueError(f"Gene {gene} not exists in GTF file.")

        df_exons = df_gene[df_gene.type == 'CDS'].copy()
        if df_exons.shape[0] == 0:
            df_exons = df_gene[df_gene.type == 'exon'].copy()
        df_exons = df_exons[df_exons.length > min_length]
        if df_exons.shape[0] == 0:
            raise ValueError(f"Gene {gene} can't found any exon records.")
        exon_cnts = df_exons.groupby(by=['chr', 'start', 'end', "length", "strand"], as_index=False).count()
        gene2exons[gene] = []
        for idx, row in exon_cnts.iterrows():
            chr_, start, end, strand = str(row['chr']), row['start'], row['end'], row['strand']
            name = '_'.join([chr_, str(start), str(end), strand])
            n_trans = row['type']
            seq = fa[chr_][start:end].seq.upper()
            if strand == '-':
                seq = reverse_complement(seq)
            exon = (name, seq, n_trans)
            gene2exons[gene].append(exon)
    return gene2exons


def self_match(probe, min_match = 4):
    length = len(probe)
    probe_re = reverse_complement(probe)
    match_pairs = 0
    for i in range(0,length-min_match+1):
        tem = probe[i:min_match+i]
        for j in range(0,length-min_match+1):
            tem_re = probe_re[j:min_match+j]
            if tem == tem_re and i + j + min_match - 2 != length:
                match_pairs = match_pairs + 1
    return match_pairs


def write_fastq(gene, seqs):
    fq = f'{TMP_DIR}/{gene}.fq'
    with open(fq, 'w') as f:
        for i, seq in enumerate(seqs):
            f.write(f"@{gene}_{i}\n")
            f.write(seq+"\n")
            f.write("+\n")
            f.write("~"*len(seq)+"\n")
    return fq


def align_se_sen(fq_path: str,
                 index: str,
                 sam_path: str,
                 threads: int = 10,
                 log_file: t.Optional[str] = 'bowtie2.log',
                 header: bool = True,
                 ) -> str:
    cmd = ["bowtie2", "-x", index, "-U", fq_path]
    if not header:
        cmd.append("--no-hd")
    cmd += ["-t", "-k", "100", "--very-sensitive-local"]
    cmd += ["-p", str(threads)]
    cmd += ["-S", sam_path]
    cmd = " ".join(cmd)
    if log_file:
        cmd += f" > {log_file} 2>&1"
    log.info(f"Call cmd: {cmd}")
    subp.check_call(cmd, shell=True)
    return sam_path


Aln = t.Tuple[str, int, int]
Block = t.Tuple[str, str, t.List[Aln]]


def read_align_blocks(
        sam_path: str
        ) -> t.Iterable[Block]:
    def yield_cond(old, rec, block, end=False):
        res = (old is not None)
        if res and not end:
            res &= rec.query_name != old.query_name
        return res
    with pysam.AlignmentFile(sam_path, mode='r') as sam:
        alns = []
        old = None
        rec = None
        for rec in sam.fetch():
            aln = rec.reference_name, rec.reference_start, rec.reference_end
            if yield_cond(old, rec, alns):
                yield old.query_name, old.query_sequence, alns
                alns = []
            if aln[0] is not None:
                alns.append(aln)
            old = rec
        if (rec is not None) and yield_cond(old, rec, alns, end=True):
            yield old.query_name, old.query_sequence, alns


def count_n_aligned_genes(sub_seqs, name, index_prefix, threads):
    fq_path = write_fastq(name, sub_seqs)
    sam_path = f"{TMP_DIR}/{name}.sam"
    if not os.path.exists(sam_path):
        align_se_sen(fq_path, index_prefix,
                     sam_path, threads=threads,
                     log_file=f"{TMP_DIR}/{name}.bowtie2.log")
    else:
        log.info("{} exists".format(sam_path))
    n_mapped_genes = []
    for _, _, alns in read_align_blocks(sam_path):
        n_genes = len(set([chr_.split("_")[0] for chr_, s, e in alns]))
        n_mapped_genes.append(n_genes)
    return n_mapped_genes


def get_sub_seq_params(tem, offset, whole_fold, barcode):
    # small is better, means RNA more unlikely to fold
    target_fold_score = -RNA.fold_compound(tem).mfe()[1]
    tem1 = tem[0:13]
    tem2 = tem[13:26]
    tem3 = tem[27:40]
    tem4 = tem[26]
    tm1 = primer3.calcTm(tem1)
    tm2 = primer3.calcTm(tem2)
    tm3 = primer3.calcTm(tem3)
    region = max(tm1,tm2,tm3) - min(tm1,tm2,tm3)
    tem1_re = reverse_complement(tem1)
    tem2_re = reverse_complement(tem2)
    tem3_re = reverse_complement(tem3)
    code1, code2 = split_barcode(barcode).split(" ")
    pad_probe = tem1_re+"CC"+code1+"TGCGTCTATTTGT"+code2+"TAGTGGAGCCT"+tem2_re
    amp_probe = tem3_re+tem4+"AGGCTCCACTA"
    match_pairs_pad = self_match(pad_probe)
    match_pairs_amp = self_match(amp_probe)
    pad_fold_score = -RNA.fold_compound(pad_probe).mfe()[1]
    amp_fold_score = -RNA.fold_compound(amp_probe).mfe()[1]
    target_fold = whole_fold[0][offset:offset+len(tem)]
    target_blocks = len(target_fold) - target_fold.count('.')  # smaller is better
    return [
        offset, pad_fold_score, amp_fold_score,
        match_pairs_pad, match_pairs_amp,
        region, tm1, tm2, tm3,
        target_fold_score, target_blocks,
        pad_probe, amp_probe
    ]


def primer_design(name, exons, index_prefix, barcode, ori_len, threads=10, min_length=40):
    """Design primers for one gene"""
    df_rows = []
    sub_seqs = []
    pool = ProcessPool(ncpus=threads)
    map_ = map if threads <= 1 else pool.map

    def process_(exon):
        exon_name, seq, n_trans = exon
        seq_len = len(seq)
        whole_fold = RNA.fold_compound(seq).mfe()
        seqs = []; rows = []
        for i in range(0,len(seq)-min_length+1):
            tem = seq[i:min_length+i]
            params = get_sub_seq_params(tem, i, whole_fold, barcode)
            row = [exon_name, n_trans] + params + [tem, split_barcode(barcode), ori_len]
            seqs.append(tem)
            rows.append(row)
        return seqs, rows

    for seqs, rows in map_(process_, exons):
        sub_seqs.extend(seqs)
        df_rows.extend(rows)

    df = pd.DataFrame(df_rows)
    df.columns = ['exon_name', 'n_trans', 'offset',
                  'pad_fold_score', 'amp_fold_score',
                  'self_match_pad', 'self_match_amp', 
                  'tm_region', 'tm1', 'tm2', 'tm3',
                  'target_fold_score', 'target_blocks',
                  'primer_pad', 'primer_amp',
                  'target_seq',
                  'barcode', 'ori_len']
    n_mapped_genes = count_n_aligned_genes(sub_seqs, name, index_prefix, threads)
    df['n_mapped_genes'] = n_mapped_genes

    return df


def build_bowtie2_index(fasta_path, index_prefix, threads=10, log_file=f"{TMP_DIR}/build.bowtie2.log"):
    cmd = ["bowtie2-build", "--threads", str(threads), fasta_path, index_prefix]
    if log_file:
        cmd += [f' > {log_file} 2>&1']
    cmd = " ".join(cmd)
    log.info(f"Call cmd: {cmd}")
    subp.check_call(cmd, shell=True)


def split_barcode(barcode):
    """
    Split ordered barcode to insert shape.
    ATGTCG -> CGA GGT
    """
    code1 = []; code2 = []
    for i in range(0, len(barcode), 2):
        c1, c2 = sorted([barcode[i], barcode[i+1]])
        code1.append(c1)
        code2.append(c2)
    code1.reverse(); code2.reverse()
    barcode_ins = " ".join([''.join(code1), ''.join(code2)])
    return barcode_ins


def parse_barcode(barcode_ins):
    """
    Parse insert shape barcode to ordered.
    CGA GGT -> ATGTCG
    """
    code1, code2 = barcode_ins.split(" ")
    code = []
    assert len(code1) == len(code2)
    for i in range(len(code1)-1, -1, -1):
        code.append("".join(sorted([code1[i],code2[i]])))
    barcode = "".join(code)
    return barcode


def read_existing_codes(path):
    codes = set()
    with open(path) as f:
        for line in f:
            code = line.strip()
            if " " in code:
                code = parse_barcode(code)
            codes.add(code)
    return codes


def filter_res(res_df, tm_range, target_fold_thresh, n_mapped_genes_thresh=4):
    for i in range(1, 4):  # tm1 tm2 tm3
        res_df = res_df[(tm_range[0] <= res_df[f"tm{i}"]) & (res_df[f"tm{i}"] <= tm_range[1])]
    res_df = res_df[res_df['n_mapped_genes'] <= n_mapped_genes_thresh]
    res_df = res_df[res_df['target_fold_score'] <= target_fold_thresh]
    return res_df

def pick_bests(res_df, best_num, gene, overlap_thresh=40, edit_dist_thresh=10):
    rows = []
    i = 0
    if res_df.shape[0] < best_num:
        log.warning( "{gene} only has {} rows, small than {}".format(
            gene, res_df.shape[0], best_num))
        return [ [gene] + list(row) for _, row in res_df.iterrows()]
    while len(rows) < best_num:
        try:
            row = res_df.iloc[i, :]
        except IndexError:
            log.warning("Failed to select best for {}, select top probes".format(gene))
            rows = [ row for _, row in res_df.iloc[:best_num, :].iterrows() ]
            break
        for o in rows:
            is_overlap = (row.exon_name == o.exon_name) and (abs(row.offset - o.offset) <= overlap_thresh)
            is_similar = editdistance.distance(row.target_seq, o.target_seq) < edit_dist_thresh
            if is_overlap or is_similar:
                break
        else:
            # not overlap with any row in rows
            rows.append(row)
        i += 1
    bests = [[gene] + list(row) for row in rows]
    return bests

def sort_res(res_df):
    df = res_df
    df['n_trans'] = - df['n_trans']
    df.sort_values(['n_mapped_genes',
                    'pad_fold_score', 'self_match_pad',
                    'amp_fold_score', 'self_match_amp',
                    'target_blocks', 'target_fold_score',
                    'tm_region', 'n_trans'],
                   inplace=True)
    df['n_trans'] = - df['n_trans']
    df = df.reset_index(drop=True)
    return df

def main(genelist, gtf, fasta,
         barcode_length=3, threads=10,
         index_prefix=None,
         existing_codes=None,
         tm_range=(36, 44),
         target_fold_thresh=10,
         output_dir="primers",
         best_num=2,
         output_raw=False,
         tmp_dir=None,
         ):
    """
    input: genelist gtf fasta
    parameters:
        barcode_length: Length of barcode.
        threads: Number of threads for bowtie2 align.
        index_prefix: Index prefix for bowtie2 align.
        existing_codes: Path to list of existing barcodes.
    output: output_dir
    """
    global TMP_DIR
    if (tmp_dir is not None) and (os.path.exists(tmp_dir)):
        TMP_DIR = tmp_dir
    else:
        TMP_DIR = get_tmp_dir("./.primer_tmp")

    if index_prefix is None:
        log.info("No bowtie2 index input, will build it.")
        from fisheye.primer_design.extract_tran_seq import extract_trans_seqs
        trans_fasta_path = f"{TMP_DIR}/transcript.fa"
        extract_trans_seqs(gtf, fasta, trans_fasta_path)
        index_prefix = f"{TMP_DIR}/transcript"
        build_bowtie2_index(trans_fasta_path, index_prefix, threads)

    fa = Fasta(fasta)
    genelist = read_gene(genelist)

    log.info("Reading gtf: " + gtf)
    df_gtf = read_gtf(gtf)
    log.info("pickup seqences..")
    log.info(f"Create tmp dir: {TMP_DIR}, fastq, sam... files will save to it")
    gene2exons = extract_exons(df_gtf, fa, genelist, min_length=40)
    if not exists(output_dir):
        os.mkdir(output_dir)

    # generate barcodes
    if existing_codes is None:
        coding_func = coding_llhc if 'score' in genelist.columns else coding_random
        barcodes, ori_lens = coding_func(genelist, barcode_length)
    else:
        log.info(f"Consider exisiting barcodes in path: {existing_codes}")
        if 'score' in genelist.columns:
            log.warning("Ignore score of genelist, due to existing barcodes.")
        existing_codes = read_existing_codes(existing_codes)
        barcodes, ori_lens = coding_random(genelist, barcode_length, existing_codes=existing_codes)

    # design probe
    best_rows = []
    idx = 0
    n_genes = len(gene2exons)
    for name, exons in gene2exons.items():
        idx += 1
        log.info(f"Designing primer for gene {name}({idx}/{n_genes}):")
        res_df = primer_design(name, exons, index_prefix, barcodes[name], ori_lens[name], threads)
        if (res_df.shape[0] > 0) and output_raw:
            out_path = join(output_dir, f"{name}.raw.csv")
            log.info("Save raw results(unfiltered) to: " + out_path)
            res_df.to_csv(out_path, index=False)
        res_df = filter_res(res_df, tm_range, target_fold_thresh)
        res_df = sort_res(res_df)
        if res_df.shape[0] > 0:
            try:
                best_rows.extend(pick_bests(res_df, best_num, name))
            except Exception:
                import ipdb; ipdb.set_trace()
            out_path = join(output_dir, f"{name}.csv")
            log.info("Save results to: " + out_path)
            res_df.to_csv(out_path, index=False)
        else:
            log.warning(f"{name} no rows pass selection.")
    best = pd.DataFrame(best_rows)
    best.columns = ['gene'] + list(res_df.columns)
    out_path = join(output_dir, "best_primer.csv")
    log.info(f"Store best primers to: {out_path}")
    best.to_csv(out_path, index=False)


fire.Fire(main)

