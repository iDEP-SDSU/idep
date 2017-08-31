# idep

[![](https://www.r-pkg.org/badges/version/plumber)](https://www.r-pkg.org/pkg/plumber)

Integrated Differential Expression and Pathway analysis (iDEP) of transcriptomic data. See documentation and manuscript. Based on annotation of 69 metazoa and 42 plant genomes in Ensembl BioMart as of 6/4/2017. Additional data from KEGG, Reactome, MSigDB (human), GSKB (mouse) and araPath (arabidopsis). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.


## Installation

> Docker Awesome List: https://github.com/veggiemonk/awesome-docker

You can install the idep by docker

### Docker-Compose

https://github.com/docker/compose

```
docker-compose up -d
```

### Build Shiny Server Docker Image

```
docker build -t ge-lab/shiny .
```


### Run Shiny Server Docker Image

```
docker run --rm -p 3838:3838 \
    -v $(pwd)/shinyapps/:/srv/shiny-server/ \
    -v $(pwd)/shinylog/:/var/log/ \
    ge-lab/shiny
```

You can bring everything down, removing the containers entirely, with the down command. Pass `--voluems` to also remove the data volume.

```
docker-compose down --volumes
```

## Documentation
https://idepsite.wordpress.com/
