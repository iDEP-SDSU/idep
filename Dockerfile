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


# Download GeneInfo Data
RUN curl -sL --retry 3 \
  "https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb" \
  | unzip \
  | tar -x -C /srv/shiny-server/idep/data/geneInfo \

# gmt
RUN curl -sL --retry 3 \
  "https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fgmt.zip?alt=media&token=5ee100d1-e645-41ef-a591-7a9ba208ce3c" \
  | unzip \
  | tar -x -C /srv/shiny-server/idep/data/gmt \

# motif
RUN curl -sL --retry 3 \
  "https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fmotif.zip?alt=media&token=dc1e5972-ffd9-43a1-bbcc-49b78da6f047" \
  | unzip \
  | tar -x -C /srv/shiny-server/idep/data/motif \

# pathwayDB
RUN curl -sL --retry 3 \
  "https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FpathwayDB.zip?alt=media&token=e602f2f7-102a-4cc4-8412-be2b05997daa" \
  | unzip \
  | tar -x -C /srv/shiny-server/idep/data/pathwayDB \

# Download GeneInfo Data


WORKDIR /usr/local/src/myscripts

CMD ["Rscript","/usr/local/src/myscripts/librarySetup.R"]
#CMD ["/usr/bin/shiny-server.sh"]
