#!/bin/bash

Rscript --quiet --vanilla /scripts/normalization.R
Rscript --quiet --vanilla /scripts/dockerizable_ml_model.R
