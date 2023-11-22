#!/bin/sh

######
# figure 1: reanalyze wildtype reporters from Tunney et al
######

R CMD BATCH analysis_code/WT/WT_mRNA_analysis.R
R CMD BATCH analysis_code/WT/WT_protein_analysis.R

R CMD BATCH figure_code/citrine_speeds_figure.R
R CMD BATCH figure_code/WT/WT_TE_plot.R
R CMD BATCH figure_code/WT/WT_mRNA_plot.R
R CMD BATCH figure_code/WT/WT_protein_plot.R

######
# figure 2: TASEP model and polysome profiles
######
# (todo: add polysome RNA-seq analysis code)
# (todo: add TASEP code and plot)

R CMD BATCH figure_code/polysome_plots.R

######
# figure 3: RQC knockouts and E2A reporter
######

# the knockouts and chimeras are in the same flow run - process the data for both
R CMD BATCH analysis_code/knockouts_chimeras_analysis.R
R CMD BATCH analysis_code/E2A_normgate.R

R CMD BATCH figure_code/knockouts_plot.R
R CMD BATCH figure_code/E2A_plot.R

######
# figure 4: start codon occlusion
######
# (todo: add TASEP code and plot)

R CMD BATCH analysis_code/5pstandardized_flow_processing_normgate.R
R CMD BATCH figure_code/5prime_standardized_plot.R

R CMD BATCH figure_code/chimera_plot.R

######
# figure 5: CRISPRi
######
R CMD BATCH analysis_code/individual_CRISPRi/crispri_growthdefectscontrols_protein_normgate.R
R CMD BATCH analysis_code/individual_CRISPRi/crispri_protein_normgate.R
R CMD BATCH analysis_code/individual_CRISPRi/crispri_te_with_errorprop.R

R CMD BATCH figure_code/individual_CRISPRi/crispri_growthdefectcontrols_plot.R
R CMD BATCH figure_code/individual_CRISPRi/crispri_protein_plot.R
R CMD BATCH figure_code/individual_CRISPRi/crispri_te_plot.R

######
# figure 6: stem loops
######


