#!/usr/bin/env Rscript 
# Runs facets in patients dir and never is launched by user, so for simplicity does not contain any arg checks. Takes patient dir as first and only argument

args = commandArgs(trailingOnly=TRUE) 
patient_folder = args[1]

library('facets')

# Run facets 

set.seed(42)
rcmat = readSnpMatrix(file.path(patient_folder, 'merged-pileup.gz'))
xx = preProcSample(rcmat)

oo = procSample(xx, cval=315)
fit = emcncf(oo)

a = c(fit$purity, fit$ploidy)
names(a) = c('purity', 'ploidy')

# write results

write.table(a, file.path(patient_folder, 'purity_ploidy.txt'), sep='\t')
write.table(fit$cncf, file.path(patient_folder, 'segments.txt'), sep='\t')

# plot and save picture 

pdf(file.path(patient_folder, 'facets_plot.pdf'))
plotSample(x=oo, emfit=fit)
dev.off()