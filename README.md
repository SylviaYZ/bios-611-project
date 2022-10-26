Analysis of Multi-omics genomics data
===============================

The datasets of interest for analysis are from National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) with accession GSE155673. The study is on "Systems biological assessment of immunity to severe and mild COVID-19 infections". The study uses a new sequencing technology (CITE-seq) that enables simultaneous sequencing of both RNA and surface protein abundance of single cells. With both measurements, it could strengthen our understanding of cells, such as differentially expressed genes across different cell types and states.

The data provided by Arunachalam et al. are sequenced on peripheral blood mononuclear cells. Exploratory analysis will be conducted on this new type of sequencing data. One of the main interests is analyzing the integrated data, that is, integrating the RNA expression with protein abundance data. Detail analyses include cell classification and differentially expressed gene identification. Further, analysis can also be performed on each measurement, and results can be compared. Last, most analyses will be using currently available methods and packages.

# Using this repository

To build the final report, one needs to follow the code provided before to build a docker image and create a container. The current version of docker image is built upon an R built for apple silicon chip. 

Step 1: build docker image

```
docker build . -t 611proj
```

Step 2: create docker container

```
docker run -v $(pwd):/home/rstudio/work\
           -p 8787:8787\
           -e PASSWORD="$(cat .password)"\
           -it 611proj
```

Step 3: open any browser and visit http://localhost:8787. Then log in to the RStudio with the password you set.

Step 4: Go to terminal inside the RStudio, change directory to work.

Step 5: Build report using Makefile. Easier just to use the code below in terminal.

```
make clean
make report.html
```