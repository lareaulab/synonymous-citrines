#!/bin/sh

######
# figure 1: reanalyze wildtype reporters from Tunney et al
######

# codon speeds of each construct
R CMD BATCH figure_code/citrine_speeds_figure.R

# flow of citrine constructs (Tunney reanalysis)
R CMD BATCH analysis_code/WT/WT_protein_analysis.R
R CMD BATCH figure_code/WT/WT_protein_plot.R

# qPCR of citrine constructs (Tunney reanalysis)
R CMD BATCH analysis_code/WT/WT_mRNA_analysis.R
R CMD BATCH figure_code/WT/WT_mRNA_plot.R

# TE of citrine constructs (Tunney reanalysis)
R CMD BATCH figure_code/WT/WT_TE_plot.R

# (todo: supp figure with re-made strains)

######
# figure 2: TASEP model and polysome profiles
######

# (todo: add TASEP code and plot)

# polysome gradients
# (todo: add polysome RNA-seq analysis code)
R CMD BATCH figure_code/polysome_plots.R

######
# figure 3: RQC knockouts and E2A reporter
######

# RQC knockouts
# (the knockouts and chimeras are in the same flow run - process the data for both)
R CMD BATCH analysis_code/knockouts_chimeras_analysis.R
R CMD BATCH figure_code/knockouts_plot.R

# 2A reporter construct
R CMD BATCH analysis_code/E2A_normgate.R
R CMD BATCH figure_code/E2A_plot.R

######
# figure 4: start codon occlusion
######

# (todo: add TASEP occlusion plot)

# 5' standardized constructs
R CMD BATCH analysis_code/5pstandardized_flow_processing_normgate.R
R CMD BATCH figure_code/5prime_standardized_plot.R

# chimera constructs
# (chimera data was analyzed in the same flow run as RQC knockouts, above)
R CMD BATCH figure_code/chimera_plot.R

######
# figure 5: CRISPRi
######

# main ciber-seq results
# (todo: add RNA-seq analysis code)
R CMD BATCH figure_code/ciber-seq_plot.R

# individual crispri confirmations (protein)
R CMD BATCH analysis_code/individual_CRISPRi/crispri_protein_normgate.R
R CMD BATCH figure_code/individual_CRISPRi/crispri_protein_plot.R

# individual crispri confirmations (mRNA + protein)
R CMD BATCH analysis_code/individual_CRISPRi/crispri_te_with_errorprop.R
R CMD BATCH figure_code/individual_CRISPRi/crispri_te_plot.R

# supp figure - ZEM citrine fusion confirming slow makes less protein

# supp figure - growth defect controls
R CMD BATCH analysis_code/individual_CRISPRi/crispri_growthdefectscontrols_protein_normgate.R
R CMD BATCH figure_code/individual_CRISPRi/crispri_growthdefectcontrols_plot.R

######
# figure 6: stem loops
######

# protein
R CMD BATCH analysis_code/stemloop/stemloop_protein_normgate.R
R CMD BATCH figure_code/stemloop/SL_protein_plot.R

# mRNA
R CMD BATCH analysis_code/stemloop/SL_mRNA_deltaCq.R
R CMD BATCH figure_code/stemloop/SL_mRNA_plots.R # makes both the mrna plot and the mrna ratio plot

# TE
R CMD BATCH figure_code/stemloop/SL_translationefficiency_plots.R # makes both the TE plot and the TE ratio plot
# this one looks a little different after my changes?


