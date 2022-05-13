import pandas as pd
import os

hlas = snakemake.params.HLAS

print(os.getcwd())
report_h = snakemake.input[0]
print(report_h)
report = pd.read_table(report_h,index_col=0)
output = snakemake.output
out_dir = snakemake.params.out_dir

os.makedirs(out_dir, exist_ok=True)

def getKmers(sequence, size):
    return [sequence[x:x+size].upper() for x in range(len(sequence) - size + 1)]


def generate_hla_query(hlas, report):
    short_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() < 9)]
    valid_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() == 9) | (report['Stripped.Sequence'].str.len() == 10)]
    long_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() > 10)]
    seq_set = set()
    for seq in long_seqs:
        seq_set.update(getKmers(seq,9))
    for seq in valid_seqs:
        seq_set.update(getKmers(seq,9))
    tiled_seq = pd.Series(list(seq_set)).to_frame()
    tiled_seq['HLA'] = None
    tiled_seq.loc[:,'HLA'] = [hlas] * tiled_seq.shape[0]
    tiled_seq = tiled_seq.explode('HLA').drop_duplicates()
    groups = tiled_seq.groupby('HLA')
    for g,frame in groups:
        out = os.path.join(out_dir,'{}.tsv'.format(g))
        frame.to_csv(out, header=None, index=False)

generate_hla_query(hlas, report)
