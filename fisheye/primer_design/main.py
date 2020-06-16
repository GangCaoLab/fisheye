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


TMP_DIR = get_tmp_dir("./primer_tmp")


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


def sequence_pickup(df_gtf, fa, genelist, min_length=40 ):
    """Pickup most represented cds"""
    seq_lst = {}
    for item in genelist.iterrows():
        name = item[1]['geneID']
        df_gene = df_gtf[df_gtf.gene_name == name]
        if df_gene.shape[0] <= 0:
            raise ValueError(f"Gene {name} not exists in GTF file.")
        df_cds = df_gene[df_gene.type == 'CDS'].copy()
        if df_cds.shape[0] == 0:
            df_cds = df_gene[df_gene.type == 'exon'].copy()
            if df_cds.shape[0] == 0:
                raise ValueError(f"Gene {name} can't found any CDS or exon records.")
        df_cds = df_cds[df_cds.length > min_length]
        cds_cnts = df_cds.groupby(by=['chr', 'start', 'end', "length", "strand"], as_index=False).count()
        cnts_max = cds_cnts[cds_cnts['type'] == cds_cnts['type'].max()]
        cds_select = cnts_max[cnts_max['length'] == cnts_max['length'].max()]
        row_ = cds_select.iloc[0]
        chr_, start, end, strand = row_.chr, row_.start, row_.end, row_.strand
        seq = fa[str(chr_)][start:end].seq
        if strand == '-':
            seq = reverse_complement(seq)
        seq_lst[name] = seq
    return seq_lst


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
        for rec in sam.fetch():
            aln = rec.reference_name, rec.reference_start, rec.reference_end
            if yield_cond(old, rec, alns):
                yield old.query_name, old.query_sequence, alns
                alns = []
            if aln[0] is not None:
                alns.append(aln)
            old = rec
        if yield_cond(old, rec, alns, end=True):
            yield old.query_name, old.query_sequence, alns


def count_n_aligned_genes(sub_seqs, name, index_prefix, threads):
    fq_path = write_fastq(name, sub_seqs)
    sam_path = align_se_sen(fq_path, index_prefix,
                            f"{TMP_DIR}/{name}.sam", threads=threads,
                            log_file=f"{TMP_DIR}/{name}.bowtie2.log")
    n_mapped_genes = []
    for _, _, alns in read_align_blocks(sam_path):
        n_genes = len(set([chr_.split("_")[0] for chr_, s, e in alns]))
        n_mapped_genes.append(n_genes)
    return n_mapped_genes


def get_sub_seq_params(tem, idx, whole_fold, barcode):
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
    target_fold = whole_fold[0][idx:idx+len(tem)]
    target_blocks = len(target_fold) - target_fold.count('.')  # smaller is better
    return [
        idx, pad_fold_score, amp_fold_score,
        match_pairs_pad, match_pairs_amp,
        region, tm1, tm2, tm3,
        target_fold_score, target_blocks,
        pad_probe, amp_probe
    ]


def primer_design(name, seq, index_prefix, barcode, ori_len, threads=10, min_length=40):
    """Design primers for one gene"""
    df_rows = []
    seq_len = len(seq)
    sub_seqs = []
    whole_fold = RNA.fold_compound(seq).mfe()
    for i in range(0,len(seq)-min_length+1):
        tem = seq[i:min_length+i]
        sub_seqs.append(tem)
        params = get_sub_seq_params(tem, i, whole_fold, barcode)
        row = params + [split_barcode(barcode), ori_len]
        df_rows.append(row)

    df = pd.DataFrame(df_rows)
    df.columns = ['offset', 'pad_fold_score', 'amp_fold_score',
                  'self_match_pad', 'self_match_amp', 
                  'tm_region', 'tm1', 'tm2', 'tm3',
                  'target_fold_score', 'target_blocks',
                  'primer_pad', 'primer_amp',
                  'barcode', 'ori_len']
    n_mapped_genes = count_n_aligned_genes(sub_seqs, name, index_prefix, threads)
    df['n_mapped_genes'] = n_mapped_genes

    df.sort_values(['n_mapped_genes',
                    'pad_fold_score', 'self_match_pad',
                    'amp_fold_score', 'self_match_amp',
                    'target_blocks', 'target_fold_score',
                    'tm_region'],
                   inplace=True)
    df = df.reset_index(drop=True)
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


def main(genelist, gtf, fasta,
         barcode_length=3, threads=10,
         index_prefix=None,
         existing_codes=None,
         output_dir="primers"):
    """
    input: genelist gtf fasta
    parameters:
        barcode_length: Length of barcode.
        threads: Number of threads for bowtie2 align.
        index_prefix: Index prefix for bowtie2 align.
        existing_codes: Path to list of existing barcodes.
    output: output_dir
    """
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
    seq_lst = sequence_pickup(df_gtf, fa, genelist, min_length=40)
    if not exists(output_dir):
        os.mkdir(output_dir)

    if existing_codes is None:
        coding_func = coding_llhc if 'score' in genelist.columns else coding_random
        barcodes, ori_lens = coding_func(genelist, barcode_length)
    else:
        log.info(f"Consider exisiting barcodes in path: {existing_codes}")
        if 'score' in genelist.columns:
            log.warning("Ignore score of genelist, due to existing barcodes.")
        existing_codes = read_existing_codes(existing_codes)
        barcodes, ori_lens = coding_random(genelist, barcode_length, existing_codes=existing_codes)

    best_rows = []
    for name, seq in seq_lst.items():
        log.info("Designing primer for gene " + name + ":")
        res_df = primer_design(name, seq, index_prefix, barcodes[name], ori_lens[name], threads)
        best_rows.append(res_df.iloc[0, :])
        out_path = join(output_dir, f"{name}.csv")
        log.info("Save results to: " + out_path)
        res_df.to_csv(out_path)
    best = pd.DataFrame(best_rows)
    out_path = join(output_dir, "best_primer.csv")
    log.info(f"Store best primers to: {out_path}")
    best.to_csv(out_path, index=False)


fire.Fire(main)

