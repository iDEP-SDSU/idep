# iDEP: Integrated Differential Expression and Pathway analysis



# What is iDep

Integrated Differential Expression and Pathway analysis (iDEP) of transcriptomic data. See documentation and manuscript. Based on annotation of 69 metazoa and 42 plant genomes in Ensembl BioMart as of 6/4/2017. Additional data from KEGG, Reactome, MSigDB (human), GSKB (mouse) and araPath (arabidopsis). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.

## iDEP is a web application hosted at http://bioinformatics.sdstate.edu/idep/ 
Local installation of this software is possible through steps below. But it is not supported or updated freqently. 

## Requirements
+ Storage should be more than 50GB
+ Memory should be more than 2GB



## Prerequisites

+ git [Installation Guide](https://gist.github.com/derhuerst/1b15ff4652a867391f03)

+ docker community edition (win10 above and MacOS)[Installation Guide]()
https://docs.docker.com/engine/installation/linux/docker-ce/debian/#install-using-the-repository

+ + docker toolbox (win7)


### Install `docker` and `docker-compose`

More detail information

+ [Installation](https://github.com/iDEP-SDSU/idep/wiki/Install-Docker-and-Docker-Compose)

### Download database files for running iDEP locally can be downloaded at https://doi.org/10.5281/zenodo.1451847.
## Quick Start

From your local machine Make clone

```
> git clone https://github.com/iDEP-SDSU/idep.git
```

The following command starts containers from docker-compose.yml. The R server listening for HTTP connections on port 3838

```
> docker-compose up -d
```

Check [localhost:3838] (http://localhost:3838/idep)


### Advanced Use (docker container base)

+ [Manage Docker Compose]()
+ [Manage Docker Image and Container]()

#### Build Docker image

You can build docker image of

```
docker build -t {name} .
```

#### Run Shiny Server Docker Image

```
docker run --rm -p 3838:3838 \
    -v $(pwd)/shinyapps/:/srv/shiny-server/ \
    -v $(pwd)/shinylog/:/var/log/ \
    idep/early
```

docker run --rm -p 3838:3838 idep/early


You can bring everything down, removing the containers entirely, with the down command. Pass `--voluems` to also remove the data volume.

```
> docker-compose down --volumes
```


## Documentation
https://idepsite.wordpress.com/
http://docs.rstudio.com/shiny-server/


## Resources
### Docker-Compoer documentation
https://docs.docker.com/compose/reference/overview/

### Shiny Server Log
https://support.rstudio.com/hc/en-us/articles/115003717168-Shiny-Server-Error-Logs

## Install required packages

open docker terminal
```
> Rscript librarySetup.R
```

## production server
docker-compose -f docker-compose.yml -f docker-compose.do.yml up -d


Proxy
https://gist.github.com/netdesk/c1db2985b542f9916995139318e5a7ce
