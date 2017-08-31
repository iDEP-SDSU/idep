FROM rocker/shiny:latest

MAINTAINER Kevin Son "eunwoo.son@sdstate.edu"

# install additional packages

COPY ./RSet /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

CMD ["Rscript", "packages.R"]
CMD ["Rscript", "bioLitePackages.R"]

CMD ["/usr/bin/shiny-server.sh"]
