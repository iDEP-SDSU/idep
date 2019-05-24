FROM rocker/shiny:latest
#FROM debian # for testing

MAINTAINER Ge lab "xijin.ge@sdstate.edu"
RUN apt-get update -qq && apt-get install -y \
  git-core \
  libcurl4-openssl-dev \
  libxml2-dev \
  libxml2  \
  libssl-dev \
  libudunits2-dev \
  libmariadbclient-dev \
  libpng-dev \
  wget \
  unzip

# COPY ./RSet /usr/local/src/myscripts
COPY ./classes /usr/local/src/myscripts
COPY ./shinyapps /srv/shiny-server

RUN mkdir -p /srv/data/geneInfo
RUN mkdir -p /srv/data/gmt
RUN mkdir -p /srv/data/motif
RUN mkdir -p /srv/data/pathwayDB
RUN mkdir -p /srv/data/data_go

WORKDIR /usr/local/src/myscripts

EXPOSE 3838

# Install required R libraries
RUN Rscript librarySetup.R
