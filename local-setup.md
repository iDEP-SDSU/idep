Set up iDEP locally
================

1.  (Optional) Download and install Docker on your workstation

<https://docs.docker.com/get-docker/>

2.  Pull the Docker image

<!-- end list -->

  - Download the latest version of the image

<!-- end list -->

``` bash
docker pull villegar/idep:latest
```

Alternatively,

``` bash
docker pull villegar/idep:R-4.0.5
```

3.  Set up global variables and create local directories

<!-- end list -->

``` bash
IDEP_PATH=/path/to/idep/idep94
IDEP_DATA=/path/to/idep/data/
LOCAL_PORT=3838

mkdir -p $IDEP_DATA/data104
```

<!-- 4. Download the latest version of iDEP from Github -->

<!-- ```bash -->

<!-- cd $IDEP_PATH -->

<!-- git clone https://github.com/iDEP-SDSU/idep -->

<!-- ``` -->

4.  Download and extract the database

<!-- end list -->

``` bash
wget https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz -O data104.tar.gz
tar xvzf data104.tar.gz --directory $IDEP_DATA/../
rm data104.tar.gz
```

<!-- ```bash -->

<!-- cd $IDEP_DATA/data104 -->

<!-- wget https://sdsu.box.com/shared/static/c24f792ojoikpzu0lkpng8uuf9ychwm7.gz -O pathwayDB.tar.gz -->

<!-- tar xvzf pathwayDB.tar.gz -->

<!-- rm pathwayDB.tar.gz -->

<!-- wget https://sdsu.box.com/shared/static/9v1ao6mwhduvrcx793j3answph9gqnkt.gz -O motif.tar.gz -->

<!-- tar xvzf motif.tar.gz -->

<!-- rm motif.tar.gz -->

<!-- wget https://sdsu.box.com/shared/static/mns0k1uvwtfnsohoc89b984ih36nmnz9.gz -O geneInfo.tar.gz -->

<!-- tar xvzf geneInfo.tar.gz -->

<!-- rm geneInfo.tar.gz -->

<!-- wget https://sdsu.box.com/shared/static/qwpdh36vcisgy1hcmadck8i8ezhvr2fh.gz -O data.go.tar.gz -->

<!-- tar xvzf data.go.tar.gz -->

<!-- rm data.go.tar.gz -->

<!-- wget https://sdsu.box.com/shared/static/sorewt7w6iypmhg2k2xhyi8myeit156o.gz -O convertIDs.db.tar.gz -->

<!-- tar xvzf convertIDs.db.tar.gz -->

<!-- rm convertIDs.db.tar.gz -->

<!-- cd $IDEP_DATA/data_go -->

<!-- wget https://raw.githubusercontent.com/villegar/idep/master/data/STRING11_species.csv -->

<!-- ``` -->

5.  Download the latest version of iDEP from GitHub:

<!-- end list -->

``` bash
git clone https://github.com/iDEP-SDSU/idep $IDEP_PATH
```

6.  Deploy a local container

<!-- end list -->

``` bash
docker run -dp $LOCAL_PORT:3838 \
  -v $IDEP_DATA:/srv/data \
  -v $IDEP_PATH/shinyapps:/srv/shiny-server \
  villegar/idep:latest
```

Alternatively,

``` bash
docker run -dp $LOCAL_PORT:3838 \
  -v $IDEP_DATA:/srv/data \
  -v $IDEP_PATH/shinyapps:/srv/shiny-server \
  villegar/idep:R-4.0.5
```

<!-- ```bash -->

<!-- LOCAL_PORT=3838 -->

<!-- docker run -dp $LOCAL_PORT:3838 \ -->

<!--   -v $IDEP_DATA:/srv/shiny-server/data \ -->

<!--   -v $IDEP_PATH:/srv/shiny-server/ \ -->

<!--   villegar/idep:latest -->

<!-- ``` -->

7.  Open a web browser and navigate to the following URL

<http://localhost:3838/idep94>
