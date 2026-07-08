# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:devel
LABEL Name=bgnormr Version=0.99.3
RUN mkdir /root/rpkg
COPY . /root/rpkg
RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
