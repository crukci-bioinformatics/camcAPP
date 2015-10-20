
FROM bioconductor/release_base
MAINTAINER Mark Dunning<mark.dunning@cruk.cam.ac.uk>

RUN sudo apt-get update
RUN sudo apt-get install -y git libssl-dev
###Get repository of the course. Install data and R packages
RUN git clone https://github.com/bioinformatics-core-shared-training/camcAPP.git /home/rstudio/
WORKDIR /home/rstudio
RUN sudo R -f installBioCPkgs.R

WORKDIR /home/rstudio

