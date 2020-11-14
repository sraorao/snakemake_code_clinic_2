# Title     : Plot Dupmetrics with python
# Objective : Script to collate duplication data from GATK MarkDuplicates and plot
# Created by: srao
# Created on: 07/11/2020

from pprint import pprint
import pandas
import matplotlib.pyplot as plt

# 'snakemake' object contains 'input', 'output', 'params', 'config'
print("----------------------------Variables in the snakemake object----------------------------")
pprint(vars(snakemake))
print("----------------------------------------------------------------------------------------")


metrics_files = snakemake.input  # a list of 4 filenames
plot_pdf = snakemake.output[0]  # extract output filename from a list of 1

# read all the dupmetrics files and concatenate them
dups = [pandas.read_csv(x, sep="\t", header=0, comment="#") for x in metrics_files]
dups = pandas.concat(dups)

# plot the percentage duplication for each sample and save to pdf
# note that MarkDuplicates reports a fraction for percentage, hence multiply with 100
dups.iloc[:, 8] = dups.iloc[:, 8] * 100
print(dups)

plt.figure()
dups.iloc[:, 8].plot.bar()
axes = plt.gca()
axes.set_ylim([0, 100])
plt.savefig(plot_pdf)

