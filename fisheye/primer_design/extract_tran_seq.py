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
    exons_df['gene_id'] = exons_df['info'].str.extract("gene_id \"(.*?)\"")
    exons_df['transcript_id'] = exons_df['info'].str.extract("transcript_id \"(.*?)\"")
    exons_df = exons_df[['chrom','left','right','strand','gene_id','transcript_id']].dropna(axis=0, how="any", subset=['transcript_id'])
    trans = {}  # (gene_id, trans_id) -> [chrom, strand, exons],  exons: (left, right)
    for (_, row) in exons_df.iterrows():
        key_ = (row['gene_id'], row['transcript_id'])
        chrom, strand, left, right = str(row['chrom']), row['strand'], row['left'], row['right']
        if key_ not in trans:
            trans[key_] = [chrom, strand, [[left, right]]]
        else:
            trans[key_][2].append([left, right])
    for key_, [chrom, strand, exons] in list(trans.items()):
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            trans[key_] = [chrom, strand, tmp_exons]
    seq_dict = {}
    for key_, [chrom, strand, exons] in list(trans.items()):
        if len(exons) == 1:
            seq = fa[chrom][exons[0][0]:exons[0][1]].seq
            if strand == '-':
                seq = reverse_complement(seq)
                seq_dict[key_] = seq
            else:
                seq_dict[key_] = seq
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
            seq_dict[key_] = seq

    log.info(f"Save results to {output_fa_path}")
    with open(output_fa_path, 'w') as f:
        for (gene_id, tran_id), seq in seq_dict.items():
            f.write(f">{gene_id}_{tran_id}\n")
            f.write(f"{seq}\n")


if __name__ == "__main__":
    fire.Fire(extract_trans_seqs)