FROM rocker/shiny:latest

MAINTAINER Kevin Son "eunwoo.son@sdstate.edu"

RUN apt-get update -qq && apt-get install -y \
  git-core \
  libcurl4-openssl-dev \
  libxml2-dev \
  libxml2  \
  libssl-dev

# install additional packages
COPY ./RSet /usr/local/src/myscripts
COPY ./classes /usr/local/src/myscripts
COPY ./shinyapps /srv/shiny-server
WORKDIR /usr/local/src/myscripts

CMD ["Rscript","/usr/local/src/myscripts/librarySetup.R"]
#CMD ["/usr/bin/shiny-server.sh"]
