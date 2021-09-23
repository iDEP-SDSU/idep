# Based on digital ocean One click app: Docker-ce

docker build ./nginx/. -t nginx  #nginx image should be build very quick

#docker pull villegar/idep:latest

docker build . -t webapp #webapp image need hours to build

echo 'Docker images have been built. Start downloading data.'

mkdir data
cd data
wget https://mft.sdstate.edu/public/file/3Y66fppA0Eym0G41taPtRw/data104.tar.gz .
tar xvzf data104.tar.gz
rm data104.tar.gz

echo 'Data has been downloaded and unziped'

echo 'iDEP docker images and databases are ready!'

