# iDEP and ShinyGO:   Graphical tools for analyzing genomics data  
## Let the data speak for itself!

[iDEP](http://ge-lab.org/idep/) is a Shiny app for analyzing RNA-seq or other transcriptomic data. See [documentation](https://idepsite.wordpress.com/) and [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6). For feedbacks or data contributions (genes and GO mapping of any species), please contact us, or visit our homepage. Send us suggestions or any error message to help improve iDEP.

As a younger sister, [ShinyGO](http://ge-lab.org/go/) focuses on the enrichment analysis of gene lists. The same database is used. They live in the same GitHub repository. ShinyGO is published [here](https://doi.org/10.1093/bioinformatics/btz931). 

## No Warranty
iDEP and ShinyGO are both under active development, and have not been fully tested. With the help of many users, we are constantly fixing bugs/errors in both source code and sometimes database. We recommend that you confirm your findings using other tools before publishing or further studies.

## Local installation
Local installation of this software is possible through steps below. But it is not supported or updated freqently. Local install is for non-profit organizations only. For-profit businesses please contact us. If you install locally, 
please check for updates regularily here and download the latest versions of both the source code and database.  For most users, we recommend to visit our website. 

## To run iDEP and ShinyGO on your local machine (Windows, MacOS, Linux):
Requirements:
+ More than 250GB available storage
+ More than 4GB memory
+ Most recent version of R and RStudio installed.

1. Upgrade to the most recent version of R and Rstudio.(Important!!)
2. Start RStudio and install all the R packages. As we need so many R packages, this may take several hours. You can let it run and get started on steps 3 and 4 below to save time. From RStudio console window:
```
source https://raw.githubusercontent.com/iDEP-SDSU/idep/master/classes/librarySetup.R
```

3. Download iDEP source code and example data files from GitHub. The best is to click the green "Code" button and select "Download ZIP" on this [page](https://github.com/iDEP-SDSU/idep). Unzip to a folder such as C:/IDEP, so that it contains all the subfolders such as config, classes, shinyapps, and so on.

4. Download all the database files from [here](http://18.235.92.206:8080/). Unzip all files to a folder (C:/IDEP/data/data104b), so that your database can be found by the most recent versions of iDEP and ShinyGO. For example, the convertIDs.db files should be at C:/IDEP/data/data104b/convertIDs.db, and the pathway information files should be at C:/IDEP/data/data104b/pathwayDB. 
Below is an example folder structure. 

<p align="center">
  <img width="313" height="299" src="docs/folders.png">
</p>


5. Start Rstudio and load the ui.R and server.R scripts in the folder C:/IDEP/shinyapps/idep95. And then click on Run app. Similarily, the ShinyGO app could be started at the folder, C:/IDEP/shinyapps/go75/. 

## To install iDEP and ShinyGO on a Linux server:

Requirements:
+ More than 250GB available storage
+ More than 4GB memory
+ A Linux system with port 80 open for web access. 

Cloud computing providers such as Amazon AWS provides temporary servers that can be used inexpensively, especially the spot instances.
See this [video](https://youtu.be/m-3vyGNYDOQ) for some step-by-step screen-shot.

1. Install [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), if it is not available.  For Ubuntu:
```
sudo apt install git-all
```

2. Install [Docker](https://docs.docker.com/get-docker/)
```
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
```

3. Install [Docker-compose](https://docs.docker.com/compose/install/)
```
sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
```

4. Copy iDEP source code from GitHub. Here we are installing iDEP in the home directory; change to other folder if desired.
```
cd ~
git clone https://github.com/iDEP-SDSU/idep.git
```
 
5. Download docker images from DockHub, and download database file from our FTP server. It can take as much as 2 hours, as the script downloads and unzips large files. Once finished total storage usage is about 161GB!
```
sudo sh idep/setup.sh
```
The setup.sh script was tested only on Ubuntu. Wait until the script shows 'iDEP docker images and databases are ready!' 

6. Start the Shiny server with Docker-compose from the idep folder.
```
cd ~/idep
sudo docker-compose up -d --scale webapp=15 
```
Now the server is running with 10 containers to serve many concurrent users. You should be able to use iDEP from a web browser with http://12.12.12.12/idep95/, where 12.12.12.12 is the IP address of the server. ShinyGO can used via http://12.12.12.12/go75/. The server's port 80 should be available and exposed.


A [user](https://github.com/wresch) has contributed scripts to install a standalone version using [Singularity](https://www.sylabs.io/). Following the instruction in this [folder.](https://github.com/iDEP-SDSU/idep/tree/master/singularity_standalone)

## Documentation
https://idepsite.wordpress.com/

