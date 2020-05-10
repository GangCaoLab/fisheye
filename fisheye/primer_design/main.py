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
    genelist = pd.read_csv(genelist,header = None, comment='#')
    genelist.columns = ['geneID']
    return genelist


def read_gtf(gtf):
    df = pd.read_csv(gtf, sep='\t', header=None, comment='#')
    df.columns = ['chr', 'source', 'type', 'start', 'end', 'X1', 'strand', 'X2', 'fields']
    df['gene_name'] = df.fields.str.extract("gene_name \"(.*?)\"")
    df['length'] = df['end'] - df['start']
    return df


def sequence_pickup(df_gtf, fa, genelist, min_length=40 ):
    seq_lst = {}
    for item in genelist.iterrows():
        name = item[1]['geneID']
        df_gene = df_gtf[df_gtf.gene_name == name]
        df_cds = df_gene[df_gene.type == 'CDS'].copy()
        df_cds = df_cds[df_cds.length > min_length]
        cds_cnts = df_cds.groupby(by=['chr', 'start', 'end', "length", "strand"], as_index=False).count()
        cnts_max = cds_cnts[cds_cnts['type'] == cds_cnts['type'].max()]
        cds_select = cnts_max[cnts_max['length'] == cnts_max['length'].max()]
        chr_, start, end, strand = cds_select.iloc[0].chr, cds_select.iloc[0].start, cds_select.iloc[0].end, \
                                   cds_select.iloc[0].strand
        seq = fa[chr_][start:end].seq
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


def cal_weight(df):
    x = df.apply(lambda x: ((x - np.min(x)) / (np.max(x) - np.min(x))))
    rows = x.index.size
    cols = x.columns.size
    k = 1.0 / math.log(rows)
    lnf = [[None] * cols for i in range(rows)]
    x = array(x)
    lnf = [[None] * cols for i in range(rows)]
    lnf = array(lnf)
    for i in range(0, rows):
        for j in range(0, cols):
            if x[i][j] == 0:
                lnfij = 0.0
            else:
                p = x[i][j] / x.sum(axis=0)[j]
                lnfij = math.log(p) * p * (-k)
            lnf[i][j] = lnfij
    lnf = pd.DataFrame(lnf)
    E = lnf
    d = 1 - E.sum(axis=0)
    w = [[None] * 1 for i in range(cols)]
    for j in range(0, cols):
        wj = d[j] / sum(d)
        w[j] = wj
    w = pd.DataFrame(w).T
    w.columns = list(df.columns)
    return w


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
                 log: t.Optional[str] = 'bowtie2.log',
                 header: bool = True,
                 ) -> str:
    cmd = ["bowtie2", "-x", index, "-U", fq_path]
    if not header:
        cmd.append("--no-hd")
    cmd += ["-t", "-k", "100", "--very-sensitive-local"]
    cmd += ["-p", str(threads)]
    cmd += ["-S", sam_path]
    cmd = " ".join(cmd)
    if log:
        cmd += f" 2> {log}"
    subp.check_call(cmd, shell=True)
    return sam_path


Aln = t.Tuple[str, int, int]
Block = t.Tuple[str, str, t.List[Aln]]


def read_align_blocks(
        sam_path: str
        ) -> t.Iterable[Block]:
    def yield_cond(old, rec, block, end=False):
        res = (old is not None) and (len(block) > 0)
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


def primer_design(name, seq, index_prefix, barcode, threads=10, min_length=40):
    df_lst = []
    seq_len = len(seq)
    sub_seqs = []
    for i in range(0,len(seq)-min_length+1):
        tem = seq[i:min_length+i]
        sub_seqs.append(tem)
        fold_score = round(RNA.fold_compound(tem).mfe()[1],2)
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
        pad_probe = tem1_re+barcode[0:3]+"CCAGTGCGTCTATTTAGTGGAGCCTGCAGT"+barcode[3:6]+tem2_re
        amp_probe = tem3_re+tem4+"ACTGCAGGCTCCA"
        match_pairs_pad = self_match(pad_probe)
        match_pairs_amp = self_match(amp_probe)
        df_lst.append([match_pairs_pad, match_pairs_amp, region, tm1, tm2, tm3, fold_score, pad_probe, amp_probe])


    fq_path = write_fastq(name, sub_seqs)
    sam_path = align_se_sen(fq_path, index_prefix,
                            f"{TMP_DIR}/{name}.sam", threads=threads,
                            log=f"{TMP_DIR}/{name}.bowtie2.log")
    n_mapped_genes = []
    for _, _, alns in read_align_blocks(sam_path):
        n_genes = len(set([chr_.split("_")[0] for chr_, s, e in alns]))
        n_mapped_genes.append(n_genes)

    df = pd.DataFrame(df_lst)
    df.columns = ['point1', 'point2', 'tm_region', 'tm1', 'tm2', 'tm3', 'RNAfold_score', 'primer1', 'primer2']
    df['n_mapped_genes'] = n_mapped_genes

    weight_df = cal_weight(df[['point1','point2','tm_region','RNAfold_score','n_mapped_genes']])
    scores = df['point1'] * weight_df['point1'].iloc[0] +\
             df['point2'] * weight_df['point2'].iloc[0] +\
             df['tm_region'] * weight_df['tm_region'].iloc[0] +\
             df['RNAfold_score'] * weight_df['RNAfold_score'].iloc[0] + \
             df['n_mapped_genes'] * weight_df['n_mapped_genes'].iloc[0]
    df['score'] = scores
    df.sort_values('score', inplace=True)
    df = df.reset_index(drop=True)
    return df


def build_bowtie2_index(fasta_path, index_prefix, threads=10):
    subp.check_call(["bowtie2-build", "--threads", str(threads),
                     fasta_path, index_prefix])

def main(genelist, gtf, fasta, barcode_length=3, threads=10, index_prefix=None, output_dir="primers"):
    """
    input: genelist gtf fasta
    output: results/{gene}.csv
    """
    if index_prefix is None:
        log.info("No bowtie2 index input, will build it.")
        from fisheye.primer_design.extract_tran_seq import extract_trans_seqs
        trans_fasta_path = f"{TMP_DIR}/index.fa"
        extract_trans_seqs(gtf, fasta, trans_fasta_path)
        index_prefix = f"{TMP_DIR}/index"
        build_bowtie2_index(trans_fasta_path, index_prefix, threads)

    fa = Fasta(fasta)
    log.info("Reading gtf: " + gtf)
    genelist = read_gene(genelist)
    df_gtf = read_gtf(gtf)
    log.info("pickup seqences..")
    log.info(f"Create tmp dir: {TMP_DIR}, fastq, sam... files will save to it")
    seq_lst = sequence_pickup(df_gtf, fa, genelist, min_length=40)
    if not exists(output_dir):
        os.mkdir(output_dir)

    import fisheye.primer_design.coding
    barcodes = coding(genelist)

    for name1, seq in seq_lst.items():
        for name2, barcode in barcodes.items():
            log.info("Designing primer for gene " + name + ":")
            res_df = primer_design(name, seq, index_prefix, barcode, threads)
            best.append(res_df.iloc[0, :])
            out_path = join(output_dir, f"{name}.csv")
            log.info("Save results to: " + out_path)
            res_df.to_csv(out_path)
    best = pd.DataFrame(best)
    out_path = join(output_dir, "best_primer.csv")
    log.info(f"Store best primers to: {out_path}")
    best.to_csv(out_path, index=False)


fire.Fire(main)

