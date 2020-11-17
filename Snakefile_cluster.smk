"""
A simple workflow to submit jobs to the cluster
inputs: Snakefile_cluster.smk (dummy input)
output: hostname.txt
"""

from snakemake.utils import min_version
min_version("5.26")

onsuccess: print("finished successfully") # insert more useful code here, e.g. send email to yourself
onerror: print("finished with errors") # insert more useful code here, e.g. send email to yourself

rule get_hostname:
    input: "Snakefile_cluster.smk"
    output: "hostname.txt"
    shell: "hostname > {output}"

