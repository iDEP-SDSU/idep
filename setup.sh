# Local install of iDEP and ShinyGO on a Linux server
# requires docker, docker-compose
#  Updated 9/23/2021 by Xijin.Ge@sdstate.edu

# idep folder should already be created by git clone.
cd idep

#Build nginx image (this is fast)
docker build ./nginx/. -t nginx  --pull

# Build iDEP main docker container; this can take 3 hours; This is now skip by using DockerHub
# First pull the most recent version of Shiny container; otherwise R won't upgrade
# docker build . -t webapp --pull

#Using the docker image hosted on DockerHub. Make sure this is up-to-date. If in doubt, build the image as above.
docker pull gexijin/idep:latest
sudo docker tag gexijin/idep webapp 


echo 'Docker images have been built. Start downloading data.'

# Download and unzip database files; this can take 2 hours
# use Amazon server for faster download


mkdir data
cd data
mkdir data104b
cd data104b


wget http://bioinformatics.sdstate.edu/data/current_version/convertIDs.db.tar.gz
tar xvzf convertIDs.db.tar.gz
rm convertIDs.db.tar.gz

wget http://bioinformatics.sdstate.edu/data/current_version/data_go.tar.gz
tar xvzf data_go.tar.gz
rm data_go.tar.gz

wget http://bioinformatics.sdstate.edu/data/current_version/geneInfo.tar.gz
tar xvzf geneInfo.tar.gz
rm geneInfo.tar.gz

wget http://bioinformatics.sdstate.edu/data/current_version/motif.tar.gz
tar xvzf motif.tar.gz
rm motif.tar.gz

wget http://bioinformatics.sdstate.edu/data/current_version/pathwayDB.tar.gz
tar xvzf pathwayDB.tar.gz
rm pathwayDB.tar.gz

#slower server
#wget --no-check-certificate https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz .


echo 'Data has been downloaded and unziped'

echo 'iDEP docker images and databases are ready!'

