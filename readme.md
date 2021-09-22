# iDEP: Integrated Differential Expression and Pathway analysis


[iDEP](http://ge-lab.org/idep/) is a Shiny app for analyzing RNA-seq or other transcriptomic data. See [documentation](https://idepsite.wordpress.com/) and [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.

## iDEP is a web application hosted at http://bioinformatics.sdstate.edu/idep/ 

## Local installation
Local installation of this software is possible through steps below. But it is not supported or updated freqently. Local install is for non-profit organizations only. For-profit businesses please contact us to license the database.

## To run iDEP on your local machine (Windows, MacOS, Linux):

1. Upgrade to the most recent version of R and Rstudio.
2. Start RStudio and install all the R packages by using [this script](https://github.com/iDEP-SDSU/idep/blob/master/classes/librarySetup.R). From RStudio console window:
```
source https://raw.githubusercontent.com/iDEP-SDSU/idep/master/classes/librarySetup.R
```
As we need so many R packages, this may take several hours.

3. Download iDEP source code and example data files from Github. The best is to click the Clone or download button on this [page](https://github.com/iDEP-SDSU/idep). And unzip to a folder such as C:/IDEP.

4. Download the most recet [database files](https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz) and unzip to the same folder (C:/IDEP), so that your database should be at C:/IDEP/data/data104. 

5. Start Rstudio and load the ui.R and server.R scripts in the folder C:/IDEP/shinyapps/idep94. And then click on Run app. Similarily, the ShinyGO app could be started at the folder, C:/IDEP/shinyapps/go74/. 

## To install iDEP on a Linux server:

Requirements:
+ Storage should be more than 200GB
+ Memory should be more than 4GB
+ A Linux system with [Docker](https://docs.docker.com/get-docker/), [Docker-compose](https://docs.docker.com/compose/install/) and [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) installed

1. Download iDEP source code to a folder.
```
mkdir ~/idep
cd ~/idep
git clone https://github.com/iDEP-SDSU/idep.git
```
2. Download database and build docker image.
For Ubuntu, you can run this [script](https://raw.githubusercontent.com/iDEP-SDSU/idep/master/docs/SetupScripts/ubuntu/setup.sh)

```
sudo sh setup.sh
```
Wait until the script shows 'iDEP docker images and databases are ready!' It can take several hours, as the script installs dozens of R packages and also copies a large database automatically.

3. Start the Shiny server with Docker-compose.
```
sudo docker-compose up -d --scale webapp=15 
```
Now the server is running. You should be able to use iDEP from a web browser with http://12.12.12.12/idep94/, where 12.12.12.12 is the IP address of the server. The server's port 80 should be exposed.

You can bring the Shiny server down, Pass `--volumes` to also remove the data volume.
```
cd ~/IDEP/
sudo docker-compose down --volumes
```

A [user](https://github.com/wresch) has contributed scripts to install a standalone version using [Singularity](https://www.sylabs.io/). Following the instruction in this [folder.](https://github.com/iDEP-SDSU/idep/tree/master/singularity_standalone)

## Documentation
https://idepsite.wordpress.com/

