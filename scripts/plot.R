# Title     : Plot Dupmetrics
# Objective : Script to collate duplication data from GATK MarkDuplicates and plot
# Created by: srao
# Created on: 07/11/2020

# 'snakemake' object contains 'input', 'output', 'params', 'config'
print("----------------------------Slots in the snakemake S4 object----------------------------")
slotNames(snakemake)
print("--------------------------Structure of the snakemake S4 object--------------------------")
str(snakemake)
print("----------------------------------------------------------------------------------------")

metrics_files = snakemake@input # a list of 4 filenames
plot_pdf = snakemake@output[[1]] # extract output filename from a list of 1

# print(metrics_files)

# read all the dupmetrics files and concatenate them
dups = lapply(metrics_files, read.table, header = TRUE, fill = TRUE)
dups = do.call(rbind, dups)

# plot the percentage duplication for each sample and save to pdf
# note that MarkDuplicates reports a fraction for percentage, hence multiply with 100
pdf(plot_pdf)
barplot(dups$PERCENT_DUPLICATION * 100, names.arg = dups$LIBRARY, ylim = c(0, 100),
        main = "Percent duplication in reads", sub = "Samples",
        col = c("lightgreen", "lightpink", "lightblue", "lightyellow"))
dev.off()
