"""
Workflow to map sequencing data from 4 runs to the genome and mark duplicates
inputs: single end fastq files in folder named fastq/
outputs:
    bam files: sorted, duplicate marked bam file for each sample
    plot: a bar plot of percentage duplication for all samples
"""
configfile: "config.yaml"

from snakemake.utils import min_version
min_version("5.26")

SAMPLES = config["SAMPLES"]

PROJECT = config["PROJECT"]

print(SAMPLES)
onsuccess: print("finished successfully") # insert more useful code here, e.g. send email to yourself
onerror: print("finished with errors") # insert more useful code here, e.g. send email to yourself

rule all:
    input:
        # merged_bam = "data/merged_bam/merged_bam.bam",
        # dupmarked_bam = expand("data/dupmarked_bam/{sample}_dupmarked.bam", sample = SAMPLES),
        # metrics = expand("data/dupmarked_bam/{sample}_dupmetrics.txt", sample = SAMPLES),
        plot = "data/plots/dups.pdf",
        python_plot = "data/plots/dups_python.pdf"

#-----------------------------------------------------------------------------------------------------------------------
def get_ref(wildcards):
    """
    Return the correct genome file based on REF_VERSION value set in config.yaml
    :param wildcards:
    :return:
    """
    if config["REF_VERSION"] == 37:
        return [config["REF37"]]
    elif config["REF_VERSION"] == 38:
        return [config["REF38"]]
    else:
        print("incorrect value for reference!")
#-----------------------------------------------------------------------------------------------------------------------

rule align:
    input:
        fastq = "data/fastq/{sample}.fastq.gz",
        ref = get_ref
    output: temp("data/sam/{sample}.sam") # sam files will be deleted after workflow finishes
    params:
        RG = "'@RG\\tID:{sample}\\tSM:{sample}_subset\\tLB:{sample}'"
    conda: "envs/bwa.yaml"
    envmodules: "BWA/0.7.17-foss-2018b"
    threads: 1
    shell: "which bwa && bwa mem -t {threads} -R {params.RG} {input.ref} {input.fastq} -o {output}"

rule sort_bams:
    input: rules.align.output
    output: "data/bam/{sample}.bam"
    conda: "envs/samtools.yaml"
    envmodules: "samtools/1.8-gcc5.4.0"
    threads: 1
    shell: "which samtools && samtools view -b {input} | samtools sort -o {output} && samtools index {output}"

rule mark_dups:
    input: rules.sort_bams.output
    output:
        bam = "data/dupmarked_bam/{sample}_dupmarked.bam",
        metrics = "data/dupmarked_bam/{sample}_dupmetrics.txt"
    conda: "envs/gatk.yaml"
    envmodules: "GATK/4.1.7.0-GCCcore-8.3.0-Java-11"
    threads: 1
    shell: "which gatk && gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

rule plot_dupmetrics:
    input: expand("data/dupmarked_bam/{sample}_dupmetrics.txt", sample = SAMPLES)
    output: "data/plots/dups.pdf"
    conda: "envs/r.yaml"
    envmodules: "R/default"
    script: "scripts/plot.R"

rule plot_dupmetrics_python:
    input: expand("data/dupmarked_bam/{sample}_dupmetrics.txt", sample = SAMPLES)
    output: "data/plots/dups_python.pdf"
    conda: "envs/python.yaml"
    script: "scripts/plot.py"
