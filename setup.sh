# Local install of iDEP and ShinyGO on a Linux server
# requires docker, docker-compose
#  Updated 9/23/2021 by Xijin.Ge@sdstate.edu

# idep folder should already be created by git clone.
cd idep

#Build nginx image (this is fast)
docker build ./nginx/. -t nginx  

# Build iDEP main docker container; this can take 3 hours; This is now skip by using DockerHub
# docker build . -t webapp 

#Using the docker image hosted on DockerHub. Make sure this is up-to-date. If in doubt, build the image as above.
docker pull gexijin/idep:latest
sudo docker tag gexijin/idep webapp


echo 'Docker images have been built. Start downloading data.'

# Download and unzip database files; this can take 2 hours
# use Amazon server for faster download
wget http://18.235.92.206:8080/data104.tar.gz 

#slower server
#wget --no-check-certificate https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz .
tar xvzf data104.tar.gz
rm data104.tar.gz

echo 'Data has been downloaded and unziped'

echo 'iDEP docker images and databases are ready!'

