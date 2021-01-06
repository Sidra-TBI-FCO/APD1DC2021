#!/bin/bash

Rscript --quiet --vanilla /scripts/normalization.R
#Rscript --quiet --vanilla /scripts/szabo_inflammation_signature.R
#Rscript --quiet --vanilla /scripts/icr_signature.R
Rscript --quiet --vanilla /scripts/pathways_signature.R
#Rscript --quiet --vanilla /scripts/MiracleScore.R
Rscript --quiet --vanilla /scripts/Prediction_dummy.R
