# iDEP: Integrated Differential Expression and Pathway analysis


[iDEP](http://ge-lab.org/idep/) is a Shiny app for analyzing RNA-seq or other transcriptomic data. See [documentation](https://idepsite.wordpress.com/) and [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6). Based on annotation of 220 animal and plant genomes in Ensembl BioMart as of 6/4/2018. Additional data from KEGG, Reactome, MSigDB (human), GSKB (mouse) and araPath (arabidopsis). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.

## iDEP is a web application hosted at http://bioinformatics.sdstate.edu/idep/ 

## Local installation
Local installation of this software is possible through steps below. But it is not supported or updated freqently. 

## To run iDEP on your laptop you will need to download the database and follow these instructions:

1. Upgrade to the most recent version of R and Rstudio.
2. Install all the R packages by running [this script](https://github.com/iDEP-SDSU/idep/blob/master/classes/librarySetup.R) under Rstudio. This can take hours.
3. Download iDEP source code and example data files from Github. The best is to click the Clone or download button on this [page](https://github.com/iDEP-SDSU/idep). And unzip to a folder such as C:/IDEP.
4. Download database files. This can take several hours. Unzip this files to a folder under the IDEP folder such as C:/IDEP/data/data92. 
[pathwayDB](https://sdsu.box.com/shared/static/c24f792ojoikpzu0lkpng8uuf9ychwm7.gz), 
[Motif](https://sdsu.box.com/shared/static/9v1ao6mwhduvrcx793j3answph9gqnkt.gz), 
[geneInfo](https://sdsu.box.com/shared/static/mns0k1uvwtfnsohoc89b984ih36nmnz9.gz), 
[data_go](https://sdsu.box.com/shared/static/qwpdh36vcisgy1hcmadck8i8ezhvr2fh.gz), 
[convertIDs](https://sdsu.box.com/shared/static/sorewt7w6iypmhg2k2xhyi8myeit156o.gz)
5. Start Rstudio and load the ui.R and server.R script. And then click on Run app. 

## To install iDEP as a server, follow the following instructions:

Requirements
+ Storage should be more than 200GB
+ Memory should be more than 2GB

A [user](https://github.com/wresch) has contributed scripts to install a standalone version using [Singularity](https://www.sylabs.io/). Following the instruction in this [folder.](https://github.com/iDEP-SDSU/idep/tree/master/singularity_standalone)  

The following are instructions based on [Docker](https://www.docker.com/).
1. Download following script based on your system:
+ For Ubuntu: [Ubuntu](https://raw.githubusercontent.com/iDEP-SDSU/idep/master/docs/SetupScripts/ubuntu/setup.sh)
Note: We are working on the script for other system.

2. Run setup script in root:
```
sudo sh setup.sh
```
Wait until the script shows 'iDEP is ready.' It can take several hours, as the script installs dozens of R pacakges and also copies a large database automatically.

3. Start system
```
sudo docker-compose up -d --scale webapp=15 
```
Now the server is running. 
Note: `webapp=15` indecate the web application count. Based on your system capacity, you can increase or decrease this number.
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
