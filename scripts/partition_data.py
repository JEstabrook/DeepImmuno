import pandas as pd
import os

hlas = snakemake.params.HLAS

report_h = snakemake.input[0]
print(report_h)
report = pd.read_table(report_h,index_col=0)
output = snakemake.output
out_dir = snakemake.params.out_dir

os.makedirs(out_dir, exist_ok=True)

def generate_hla_query(hlas, report):
    short_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() < 9)]
    valid_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() == 9) | (report['Stripped.Sequence'].str.len() == 10)].to_frame()
    long_seqs = report['Stripped.Sequence'][(report['Stripped.Sequence'].str.len() > 10)]
    valid_seqs['HLA'] = None
    valid_seqs.loc[:,'HLA'] = [hlas] * valid_seqs.shape[0]
    valid_seqs = valid_seqs.explode('HLA').drop_duplicates()
    groups = valid_seqs.groupby('HLA')
    for g,frame in groups:
        frame.to_csv(os.path.join(out_dir,'{}.tsv'.format(g)),header = None,index=False)

generate_hla_query(hlas, report)
