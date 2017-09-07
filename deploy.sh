#!/bin/bash

mkdir -p /opt/shiny/camcapp/logs
chmod ugo+rwx /opt/shiny/camcapp/logs
docker run -u shiny --rm -p 80:3838 crukci-bioinformatics/camcapp:latest &
sleep 2
docker ps -a
