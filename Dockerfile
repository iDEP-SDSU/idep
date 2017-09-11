FROM rocker/shiny:latest
#FROM debian # for test cases

MAINTAINER Kevin Son "eunwoo.son@sdstate.edu"
RUN apt-get update -qq && apt-get install -y \
  git-core \
  libcurl4-openssl-dev \
  libxml2-dev \
  libxml2  \
  libssl-dev \
  wget \ 
  unzip

COPY ./RSet /usr/local/src/myscripts
COPY ./classes /usr/local/src/myscripts
COPY ./shinyapps /srv/shiny-server

RUN mkdir -p /srv/shiny-server/idep/data/geneInfo
RUN mkdir -p /srv/shiny-server/idep/data/gmt
RUN mkdir -p /srv/shiny-server/idep/data/motif
RUN mkdir -p /srv/shiny-server/idep/data/pathwayDB

# Download Required Data
RUN wget -qO- -O tmp.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb'\
  && unzip tmp.zip -d /srv/shiny-server/idep/data && rm tmp.zip
RUN wget -qO- -O tmp2.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fgmt.zip?alt=media&token=5ee100d1-e645-41ef-a591-7a9ba208ce3c' \
  && unzip tmp2.zip -d /srv/shiny-server/idep/data && rm tmp2.zip
RUN wget -qO- -O tmp3.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fmotif.zip?alt=media&token=dc1e5972-ffd9-43a1-bbcc-49b78da6f047' \
  && unzip tmp3.zip -d /srv/shiny-server/idep/data && rm tmp3.zip
RUN wget -qO- -O tmp4.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FpathwayDB.zip?alt=media&token=e602f2f7-102a-4cc4-8412-be2b05997daa' \
  && unzip tmp4.zip -d /srv/shiny-server/idep/data && rm tmp4.zip
RUN wget -qO- -O tmp5.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FconvertIDs.db.zip?alt=media&token=55c80e8c-d5c1-43a5-995f-5e0c56242013' \
  && unzip tmp5.zip -d /srv/shiny-server/idep/data && rm tmp5.zip

WORKDIR /usr/local/src/myscripts

# Install required R libraries
CMD ["rscript","/usr/local/src/myscripts/librarySetup.R"] 

#CMD ["/usr/bin/shiny-server.sh"] #If you don't use docker-compose need to comment out
