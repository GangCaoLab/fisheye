import pandas as pd
from pyfaidx import Fasta
from fire import Fire

from fisheye.utils import reverse_complement
from fisheye.utils import get_logger


log = get_logger(__name__)


def extract_trans_seqs(gtf_path, fa_path, output_fa_path):
    """Extract transcript's sequences.
    """
    log.info(f"Extract transcript sequences from: {gtf_path}, {fa_path}")
    fa = Fasta(fa_path)
    gtf_df = pd.read_csv(gtf_path,sep='\t',header=None,comment='#')
    gtf_df.columns = ['chrom', 'source', 'feature', 'left', 'right', 'score', \
                    'strand', 'frame', 'info']
    exons_df = gtf_df[gtf_df.feature == 'exon']
    exons_df = exons_df[exons_df.left < exons_df.right]
    exons_df['gene_id'] = exons_df['info'].str.split(';',expand = True)[0].str.split(' ',expand = True)[1].str.split('"',expand = True)[1]
    exons_df['transcript_id'] = exons_df['info'].str.split(';',expand = True)[1].str.split(' ',expand = True)[2].str.split('"',expand = True)[1]
    exons_df = exons_df[['chrom','left','right','strand','gene_id','transcript_id']].dropna(axis=0, how="any", subset=['transcript_id'])
    transid_to_geneid = {}
    trans = {}  # trans_id -> [chrom, strand, exons],  exons: (left, right)
    for (_, row) in exons_df.iterrows():
        transcript_id = row['transcript_id']
        chrom, strand, left, right = row['chrom'],row['strand'],row['left'],row['right']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
        else:
            trans[transcript_id][2].append([left, right])
        transid_to_geneid[row['transcript_id']] = row['gene_id']
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

    log.info(f"Save results to {output_fa_path}")
    with open(output_fa_path, 'w') as f:
        for tran, seq in seq_dict.items():
            gene_id = transid_to_geneid[tran]
            f.write(f">{gene_id}_{tran}\n")
            f.write(f"{seq}\n")


if __name__ == "__main__":
    fire.Fire(extract_trans_seqs)