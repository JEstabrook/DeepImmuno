import pandas as pd

#HLAS, = glob_wildcards("samples/input_data/{hla}.tsv")
hla_types = pd.read_table('data/hla2paratopeTable_aligned.txt')
HLAS = hla_types['HLA'].unique().tolist()
HLAS.remove('HLA-B*0602')
HLAS.remove('HLA-B*9234')

rule all:
    input:
        expand("samples/results/{hla}/deepimmuno-cnn-result.csv", hla=HLAS),
        expand("samples/results/{hla}/mhcflurry_result.tsv", hla=HLAS),
        expand("samples/tiled_results/{hla}/mhcflurry_result.tsv", hla=HLAS)

rule partition_data:
    input:
        "input/Final_comparisons.tsv"
    output:
        expand("samples/input_data/{hla}.tsv",hla=HLAS)
    params:
        HLAS = HLAS,
        out_dir = "samples/input_data"
    script:
         "./scripts/partition_data.py" 

rule run_deepimmuno:
    input:
        "samples/input_data/{hla}.tsv",
    output:
        "samples/results/{hla}/deepimmuno-cnn-result.txt"
    params:
        out_dir = lambda wildcards: "samples/results/{hla}/".format(hla=wildcards.hla),
    shell:
        """
        set +u
        source deactivate
        conda activate deepimmuno_env
        mkdir -p {params.out_dir}
        python deepimmuno-cnn.py --mode "multiple" --intdir "{input}" --outdir "{params.out_dir}"
        """

rule convert_to_csv:
    input:
        "samples/results/{hla}/deepimmuno-cnn-result.txt"
    output:
        "samples/results/{hla}/deepimmuno-cnn-result.csv"
    shell:
        """
        cat {input} | tr "\\t" "," > {output}
        """


rule run_mhcflurry:
    input:
        "samples/results/{hla}/deepimmuno-cnn-result.csv"
    output:
        "samples/results/{hla}/mhcflurry_result.tsv"
    shell:
        """
        set +u
        source deactivate
        conda activate mhcflurry_test_env 
        mhcflurry-predict {input} --allele-column HLA --peptide-column peptide --out {output} --output-delimiter '\t'
        """

#### Adding tiling

rule partition_tiled_data:
    input:
        "input/Final_comparisons.tsv"
    output:
        expand("samples/input_tiled_data/{hla}.tsv",hla=HLAS)
    params:
        HLAS = HLAS,
        out_dir = "samples/input_tiled_data"
    script:
         "./scripts/partition_tiled_data.py"

rule run_tiled_deepimmuno:
    input:
        "samples/input_tiled_data/{hla}.tsv",
    output:
        "samples/tiled_results/{hla}/deepimmuno-cnn-result.txt"
    params:
        out_dir = lambda wildcards: "samples/tiled_results/{hla}/".format(hla=wildcards.hla),
    shell:
        """
        set +u
        source deactivate
        conda activate deepimmuno_env
        mkdir -p {params.out_dir}
        python deepimmuno-cnn.py --mode "multiple" --intdir "{input}" --outdir "{params.out_dir}"
        """

rule convert_tiled_to_csv:
    input:
        "samples/tiled_results/{hla}/deepimmuno-cnn-result.txt"
    output:
        "samples/tiled_results/{hla}/deepimmuno-cnn-result.csv"
    shell:
        """
        cat {input} | tr "\\t" "," > {output}
        """


rule run_tiled_mhcflurry:
    input:
        "samples/tiled_results/{hla}/deepimmuno-cnn-result.csv"
    output:
        "samples/tiled_results/{hla}/mhcflurry_result.tsv"
    shell:
        """
        set +u
        source deactivate
        conda activate mhcflurry_test_env
        mhcflurry-predict {input} --allele-column HLA --peptide-column peptide --out {output} --output-delimiter '\t'
        """


