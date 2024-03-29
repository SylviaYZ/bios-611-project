Exploratory analysis of multi-omics  data
================================

The datasets of interest for analysis are from National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) with accession GSE155673. The study is on "Systems biological assessment of immunity to severe and mild COVID-19 infections". The study uses a new sequencing technology (CITE-seq) that enables simultaneous sequencing of both RNA and surface protein abundance of single cells. With both measurements, it could strengthen our understanding of cells, such as differentially expressed genes across different cell types and states.

The data provided by Arunachalam et al. are sequenced on peripheral blood mononuclear cells. Exploratory analysis will be conducted on this new type of sequencing data. One of the main interests is to explore the additional information antibody data provides as comparing to single-cell RNA-seq. Detail analyses include cell classification and differentially expressed gene identification. Analyses are performed on each type of measurement, and results are compared. Last, most analyses will be using currently available methods and packages in R.

# Using this repository

To build the final report, user needs to follow the code provided before to build a docker image and create a container. 

Step 1: build docker image

The current version of docker image is built upon an R built for apple silicon processor.

For intel processor users, please change the first line in Dockerfile to "FROM rocker/verse" to build the docker image. 

```
docker build . -t 611proj
```

Step 2: create password for launching docker container.

Create a .password file containing the password user wishes to set.

Step 3: create docker container

```
docker run -v $(pwd):/home/rstudio/work\
           -p 8787:8787\
           -e PASSWORD="$(cat .password)"\
           -it 611proj
```

Step 4: open any browser and visit http://localhost:8787. Then log in to the RStudio with the password user set.

Step 5: Go to terminal inside RStudio, change directory to work.

Step 6: Build report using Makefile with code provided below.

```
make clean
make report.html
```

# Conclusion


After pre-processing, celltype specific markers in RNA and antibody data behave consistently. Antibody data can marginally improve cell type identification of single cell data.

Figures below show annotation of single-cell data using RNA and antibody information. Given the uncertainty in UMAP visualization, these two plots are not replicable using Makefile. 

![UMAP celltype RNA](./UMAP_celltypeRNA.png)

![UMAP celltype RNA](./UMAP_celltypeADT.png)

