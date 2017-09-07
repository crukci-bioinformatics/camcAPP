#!/bin/bash

git clone https://github.com/crukci-bioinformatics/camcAPP.git camcAPP

chmod ugo+r camcAPP/*.R
chmod ugo+r camcAPP/www/*
chmod ugo+r shiny-server.conf

# build docker image
docker build --tag="crukci-bioinformatics/camcapp" Dockerfile-server

# remove dangling/untagged images
#docker rmi $(docker images --filter "dangling=true" -q --no-trunc)


## To test on local machine, go to:-
## localhost/G4Hunter