import os
from os.path import join, exists
import fire
import numpy as np
import math
from numpy import array
import pandas as pd
from pyfaidx import Fasta
import RNA
import primer3
import glob
from fisheye.utils import get_logger


log = get_logger(__name__)


def read_gene(genelist):
    genelist = pd.read_csv(genelist,header = None)
    genelist.columns = ['geneID']
    return genelist

def read_gtf(gtf):
    df = pd.read_csv(gtf, sep='\t', header=None)
    df.columns = ['chr', 'source', 'type', 'start', 'end', 'X1', 'strand', 'X2', 'fields']
    df['gene_name'] = df.fields.str.extract("gene_name \"(.*?)\"")
    df['length'] = df['end'] - df['start']
    return df

###################################################
def get_base_map():
    base_map = [b'\0' for i in range(256)]
    base_map[ ord('A') ] = b'T'
    base_map[ ord('T') ] = b'A'
    base_map[ ord('C') ] = b'G'
    base_map[ ord('G') ] = b'C'
    base_map[ ord('a') ] = b't'
    base_map[ ord('t') ] = b'a'
    base_map[ ord('c') ] = b'g'
    base_map[ ord('g') ] = b'c'
    base_map[ ord('N') ] = b'N'
    base_map[ ord('n') ] = b'n'
    base_map = bytes(b''.join(base_map))
    return base_map

BASEMAP = get_base_map()
def reverse_complement(seq):
    res = seq[::-1]
    res = res.translate(BASEMAP)
    return res

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

#####################################################
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
    '''熵值法计算变量的权重'''
    # 标准化
    x = df.apply(lambda x: ((x - np.min(x)) / (np.max(x) - np.min(x))))
    # 求k
    rows = x.index.size  # 行
    cols = x.columns.size  # 列
    k = 1.0 / math.log(rows)

    lnf = [[None] * cols for i in range(rows)]

    # 信息熵
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

    # 计算冗余度
    d = 1 - E.sum(axis=0)
    # 计算各指标的权重
    w = [[None] * 1 for i in range(cols)]
    for j in range(0, cols):
        wj = d[j] / sum(d)
        w[j] = wj
    w = pd.DataFrame(w).T
    w.columns = ['point1', 'point2', 'tm_region', 'RNAfold_score']
    return w

def primer_design(name, seq, min_length=40):
    df_lst = []
    seq_len = len(seq)
    for i in range(0,len(seq)-min_length+1):
        tem = seq[i:min_length+i]
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
        pad_probe = tem1_re+"CCAGTGCGTCTATTTAGTGGAGCCTGCAGT"+tem2_re
        amp_probe = tem3_re+tem4+"ACTGCAGGCTCCA"
        match_pairs_pad = self_match(pad_probe)
        match_pairs_amp = self_match(amp_probe)
        with open(f'{name}.fa', 'a+') as f:
            f.write('>' + name + '\t' + str(i) + '\n' + tem + '\n')
        df_lst.append([match_pairs_pad, match_pairs_amp, region, tm1, tm2, tm3, fold_score, pad_probe, amp_probe])
    df = pd.DataFrame(df_lst)
    df.columns = ['point1', 'point2', 'tm_region', 'tm1', 'tm2', 'tm3', 'RNAfold_score', 'primer1', 'primer2']

    weight_df = cal_weight(df[['point1','point2','tm_region','RNAfold_score']])
    scores = df['point1'] * weight_df['point1'].iloc[0] +\
             df['point2'] * weight_df['point2'].iloc[0] +\
             df['tm_region'] * weight_df['tm_region'].iloc[0] +\
             df['RNAfold_score'] * weight_df['RNAfold_score'].iloc[0]
    df['score'] = scores
    df.sort_values('score', inplace=True)
    df = df.reset_index(drop=True)
    return df

def write_fq():
    for i in glob.glob('*.fa'):
        with open(i, 'r') as f1:
            N = 0
            with open("%s.fq" % i[:-3], 'w') as f2:
                for line in f1:
                    N += 1
                    if N % 2 == 1:
                        f2.write('@' + line[1:])
                    else:
                        length = len(line[5:-5])
                        f2.write(line + '+' + '\n' + '>' * 5 + 'G' * length + 'B' * 4 + '\n')

def main(genelist, gtf, fasta, output_dir="primers"):
    """
    input: genelist gtf fasta
    output: results/{gene}.csv
    """
    fa = Fasta(fasta)
    log.info("Reading gtf: " + gtf)
    genelist = read_gene(genelist)
    df_gtf = read_gtf(gtf)
    log.info("pickup seqences..")
    seq_lst = sequence_pickup(df_gtf, fa, genelist, min_length=40)
    # TODO: calculate specificity:
    # 1. save all gene's all transcripts to fasta
    #    make aligner index
    # 2. save all sub seqs(40bp) of seq in seq_lst to fastq
    # 3. align fastq to aligner index get sam file
    # 4. iterate sam file, count out-of region aligns.
    if not exists(output_dir):
        os.mkdir(output_dir)
    best = []
    for name, seq in seq_lst.items():
        log.info("Designing primer for gene " + name + ":")
        res_df = primer_design(name,seq)
        write_fq()
        best.append(res_df.iloc[0, :])
        out_path = join(output_dir, f"{name}.csv")
        log.info("Save results to: " + out_path)
        res_df.to_csv(out_path)
    best = pd.DataFrame(best)
    out_path = join(output_dir, "best_primer.csv")
    log.info(f"Store best primers to: {out_path}")
    best.to_csv(out_path, index=False)


fire.Fire(main)

