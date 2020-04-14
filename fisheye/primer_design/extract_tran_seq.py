import pandas as pd
from pyfaidx import Fasta

gtf_file = 'gencode.vM24.annotation.gtf'
fa = Fasta('GRCm38.primary_assembly.genome.fa')

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

def extract_exons(gtf_file,fa):
    gtf_df = pd.read_csv(gtf_file,sep='\t',header=None)
    gtf_df.columns = ['chrom', 'source', 'feature', 'left', 'right', 'score', \
                    'strand', 'frame', 'info']
    exons_df = gtf_df[gtf_df.feature == 'exon']
    exons_df = exons_df[exons_df.left < exons_df.right]
    exons_df['gene_id'] = exons_df['info'].str.split(';',expand = True)[0].str.split(' ',expand = True)[1].str.split('"',expand = True)[1]
    exons_df['transcript_id'] = exons_df['info'].str.split(';',expand = True)[1].str.split(' ',expand = True)[2].str.split('"',expand = True)[1]
    exons_df = exons_df[['chrom','left','right','strand','gene_id','transcript_id']].dropna(axis=0, how="any", subset=['transcript_id'])
    trans = {}
    for item in exons_df.iterrows():
        transcript_id = item[1]['transcript_id']
        chrom, strand, left, right = item[1]['chrom'],item[1]['strand'],item[1]['left'],item[1]['right']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
        else:
            trans[transcript_id][2].append([left, right])
    for tran, [chrom, strand, exons] in list(trans.items()):
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            trans[tran] = [chrom, strand, tmp_exons]
    seq_dict = {}
    for tran, [chrom, strand, exons] in list(trans.items()):
        if len(exons) == 1:
            seq = fa[chrom][exons[0][0]:exons[0][1]].seq
            if strand == '-':
                seq = reverse_complement(seq)
                seq_dict[tran] = seq
            else:
                seq_dict[tran] = seq
        else:
            seq_lst = []
            for i in range(0, len(exons)):
                seq = fa[chrom][exons[i][0]:exons[i][1]].seq
                if strand == '-':
                    seq = reverse_complement(seq)
                    seq_lst.append(seq)
                else:
                    seq_lst.append(seq)
            seq = "".join(seq_lst)
            seq_dict[tran] = seq
    with open('GRCm38.primary_assembly.transcript.fa', 'w') as f:
        for tran, seq in seq_dict.items():
            f.write('>' + tran + '\n' + seq + '\n')
    pass