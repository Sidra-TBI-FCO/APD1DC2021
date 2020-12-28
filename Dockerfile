FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.7 python3-pip python3-setuptools python3-dev python3-pip libcurl4-gnutls-dev libxml2-dev

RUN R -e "install.packages('data.table', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('NOISeq')"
RUN R -e "BiocManager::install('GSVA')"
RUN R -e "BiocManager::install('gclus')"
RUN R -e "BiocManager::install('ConsensusClusterPlus')"
RUN R -e "BiocManager::install('clue')"

RUN pip3 install -U scikit-learn notebook pandas numpy jupytext seaborn xgboost shap scipy keras

RUN  mkdir -p /scripts/ && mkdir -p /data/ && mkdir -p /output/&& mkdir -p /sub_pipelines/

COPY Model_building/Scripts/szabo_inflammation_signature.R /scripts/
COPY Model_building/Scripts/icr_signature.R /scripts/
COPY Model_building/Scripts/pathways_signature.R /scripts/
COPY Model_building/Required_Files/ICR_genes.RData /data/
COPY Model_building/Required_Files/Selected.pathways.3.4.RData /data/
COPY Synthetic_Data/GRCh37ERCC_refseq105_genes_count.csv /data/

COPY Model_building/Docker_files/baseline_pipeline.sh /sub_pipelines/
COPY Model_building/Docker_files/ICRscore_pipeline.sh /sub_pipelines/
COPY Model_building/Docker_files/pathwayssignature_pipeline.sh /sub_pipelines/
COPY Model_building/Docker_files/pipeline.sh /

ENTRYPOINT ["/pipeline.sh"]