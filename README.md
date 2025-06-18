# ClipMixcure

R language statistical function package for calculating profile likelihood based inference or Wald statistic based inference (Rubin's Rule) for multiply imputed data (obtained via mice) under mixture cure model with Firth penalized likelihood or maximum likelihood for point estimation. The package consists of four functions for i) generating multiple parameter estimates (maximum likelihood or penalized likelihood) by fitting mixture cure model to individual multiply imputed datasets ('mixcure.penal.mi.R'); ii) generating Rubin's Rule based inference for pooling multiple estimates ('pool.mixcure.R'); and iii) generating combined profile likelihood inference for pooling multiple estimates ('clip.mixcure.R' and 'mixcure.clip.pdf.r'). 

Created by Changchang Xu

Contact:changchang.xu@alumni.utoronto.ca

This package can be installed via the following R code:

devtools::install_github("ChangchangXu-LTRI/ClipMixcure", build = TRUE, build_opts = c())

library(ClipMixcure)
