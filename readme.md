# iDEP: Integrated Differential Expression and Pathway analysis


[iDEP](http://ge-lab.org/idep/) is a Shiny app for analyzing RNA-seq or other transcriptomic data. See [documentation](https://idepsite.wordpress.com/) and [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6). Based on annotation of 220 animal and plant genomes in Ensembl BioMart as of 6/4/2018. Additional data from KEGG, Reactome, MSigDB (human), GSKB (mouse) and araPath (arabidopsis). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.

## iDEP is a web application hosted at http://bioinformatics.sdstate.edu/idep/ 
Local installation of this software is possible through steps below. But it is not supported or updated freqently. 

## Requirements
+ Storage should be more than 200GB
+ Memory should be more than 2GB

## Setup

1. Download following script based on your system:
+ Ubuntu: [Ubuntu](https://idep-sdsu.github.io/idep/SetupScript/ubuntu/setup.sh)

2. Run setup script in root:
   
```
sudo sh setup.sh
```
Wait until the script shows 'iDep is ready.'

## Start system

```
sudo docker-compose up -d --scale webapp=15 
```

Note: `webapp=15` indecate the web application count. Based on your system capacity, you can increase or decrease this number.

## Stop system

You can bring everything down, removing the containers entirely, with the down command. Pass `--voluems` to also remove the data volume.

```
sudo docker-compose down --volumes
```


## Documentation
https://idepsite.wordpress.com/
http://docs.rstudio.com/shiny-server/


## Resources
### Docker-Compoer documentation
https://docs.docker.com/compose/reference/overview/
