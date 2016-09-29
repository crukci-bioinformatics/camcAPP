
FROM bioconductor/release_base
MAINTAINER Mark Dunning<mark.dunning@cruk.cam.ac.uk>

RUN sudo apt-get update
RUN sudo apt-get install -y git libssl-dev
###Get repository of the course. Install data and R packages
RUN git clone https://github.com/bioinformatics-core-shared-training/camcAPP.git /tmp/camcAPP
WORKDIR /tmp/camcAPP
RUN sudo R -f installBioCPkgs.R
RUN sudo R -f createSQL.R
RUN cp server.R /home/rstudio
RUN cp UI.R /home/rstudio
RUN *.sqlite3 /home/rstudio
WORKDIR /home/rstudio

