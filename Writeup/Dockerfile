FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.7 python3-pip python3-setuptools python3-dev python3-pip libcurl4-gnutls-dev libxml2-dev libpng-dev libjpeg-dev libblas-dev liblapack-dev r-bioc-genefilter liblzma-dev libssl-dev libbz2-dev gfortran r-cran-devtools

RUN R -e "install.packages('data.table', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('scales', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
#RUN R -e "BiocManager::install('NOISeq')"
#RUN R -e "BiocManager::install('gclus')"
RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('caret', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('AppliedPredictiveModeling', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('MLmetrics', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('gbm', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('randomForest', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('kernlab', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('bst', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('e1071', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('xgboost',repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_github('miccec/yaGST')"
RUN R -e "devtools::install_github('tolgaturan-github/Miracle')"

#RUN pip3 install -U scikit-learn notebook pandas numpy jupytext seaborn xgboost shap scipy keras

RUN  mkdir -p /scripts/ && mkdir -p /data1/ && mkdir -p /output/

#RUN R -e "install.packages('import', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('klaR',repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('matrixStats', repos = 'http://cran.us.r-project.org')"

COPY Input_files/clinical_data.csv /data1/
COPY Input_files/GRCh37ERCC_refseq105_genes_count.csv /data1/
COPY Input_files/Selelected_Path_VariousGeneIDs.RData /data1/
COPY Input_files/geneInfo.July2017.RData /data1/
COPY Input_files/ICR_genes.RData /data1/
COPY Input_files/revised_normalization.R /scripts/
COPY tmb_pdl1_submission.R /scripts/

COPY pipeline1.sh /

ENTRYPOINT ["bash","/pipeline1.sh"]
