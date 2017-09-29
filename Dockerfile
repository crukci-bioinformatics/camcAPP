FROM bioconductor/release_base
MAINTAINER Mark Dunning<mark.dunning@cruk.cam.ac.uk>
RUN rm -rf /var/lib/apt/lists/*
RUN sudo apt-get update
RUN sudo apt-get install -y git libssl-dev
###Get repository of the app. Install data and R packages
WORKDIR /home/
RUN git clone https://github.com/bioinformatics-core-shared-training/camcAPP.git
WORKDIR camcAPP
RUN sudo R -f installBioCPkgs.R
RUN sudo R -f createSQL.R
