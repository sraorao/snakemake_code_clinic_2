"""
Workflow to merge bam files
inputs: 4 bam files from folder /data/bam/
outputs: 1 merged bam file in folder /data/merged_bam/
"""
configfile: "config.yaml"

SAMPLES = config["SAMPLES"]

print(SAMPLES)
onsuccess: print("finished successfully") # insert more useful code here, e.g. send email to yourself
onerror: print("finished with errors") # insert more useful code here, e.g. send email to yourself


rule bedtools_version:
    input: expand("data/bam/{sample}.bam", sample = SAMPLES)
    output: "data/version.txt"
    container:
        "docker://biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    threads: 1
    shell: "bedtools --version > {output}"