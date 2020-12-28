# APD1DC2021 
Anti-PD1 Dream Challenge 2020-2021

Folder stucture Github :

-Model_building <br />
--Proccessed_data <br />
---Riaz <br />
---... <br />
--Required_files (gene list and gene info files, in the Rdata format if possible) <br />
--Scripts <br />
--Docker_files (the pipelines and the alterantive docker build files) <br />
-Synthetic_data (fake/synthetis data with the exact file names and structure of the input data) <br />
-Output (our prediction scores) <br />
Dockerfile (the file to build the docker)

Folder stucture Docker :

-scripts (the script) <br />
-data (input data & required files) <br />
-output (the output) <br />
-sub_pipelines (modules of the pipeline) <br />
pipeline.sh (the entry point, gets run when the docker is run) <br />
