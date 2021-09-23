# Local install of iDEP and ShinyGO on a Linux server
# requires docker, docker-compose
#  Updated 9/23/2021 by Xijin.Ge@sdstate.edu

# copy the source code of iDEP from GitHub
git clone https://github.com/iDEP-SDSU/idep.git


#Build nginx image (this is fast)
docker build ./nginx/. -t nginx  

# Build iDEP main docker container; this can take 3 hours
#docker pull villegar/idep:latest
docker build . -t webapp 

echo 'Docker images have been built. Start downloading data.'

# Download and unzip database files; this can take 2 hours
wget --no-check-certificate https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz .
tar xvzf data104.tar.gz
rm data104.tar.gz

echo 'Data has been downloaded and unziped'

echo 'iDEP docker images and databases are ready!'

