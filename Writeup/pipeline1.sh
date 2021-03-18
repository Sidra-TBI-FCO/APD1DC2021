#!/bin/bash

Rscript --quiet --vanilla /scripts/revised_normalization.R
Rscript --quiet --vanilla /scripts/tmb_pdl1_submission.R
