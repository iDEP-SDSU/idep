# Ubuntu setup 

echo 'Checking required software.'
DEBIAN_FRONTEND=noninteractive

if ! [ -x "$(command -v R)" ]; then
	echo 'R is not installed. Installing r-base.'
	apt-get -qq -y update && -qq -y upgrade
	apt-get -y -qq install r-base
fi


if ! [ -x "$(command -v git)" ]; then
	echo 'Git is not installed. Installing git.'
	apt-get -qq -y update && -qq -y upgrade
	apt-get -y -qq install git
fi


if ! [ -x "$(command -v docker)" ]; then
	echo 'Docker is not installed. Installing docker.io and docker-compose.'
	apt-get -qq -y update && -qq -y upgrade
	apt-get -y -qq install docker.io
	apt-get -y -qq install docker-compose
fi
echo ''
echo 'Software setup done.'
echo ''
echo ''

echo 'Clone code from GitHub'
git clone https://github.com/iDEP-SDSU/idep.git -q 
echo ''
echo 'Code Clone Done.'
echo ''
echo ''


echo 'Start building docker image. This process may take hours.'
nohup docker daemon -H tcp://0.0.0.0:2375 -H unix:///var/run/docker.sock
cd idep
docker build -q ./nginx/. -t nginx  #nginx image should be build very quick
docker build -q . -t webapp #webapp image need hours to build
echo 'Docker images have been built. Start build data files.'
echo ''
echo ''

mkdir data
cd data
mkdir data92
cd data92

echo 'Downloading Data(1/5): pathwayDB'
wget -q --show-progress https://sdsu.box.com/shared/static/c24f792ojoikpzu0lkpng8uuf9ychwm7.gz -O pathwayDB.tar.gz
tar xvzf pathwayDB.tar.gz
rm pathwayDB.tar.gz
echo 'Download and Unzip Finished (1/5): pathwayDB'

echo 'Downloading Data(2/5): motif'
wget https://sdsu.box.com/shared/static/9v1ao6mwhduvrcx793j3answph9gqnkt.gz -O motif.tar.gz
tar xvzf motif.tar.gz
rm motif.tar.gz
echo 'Download and Unzip Finished (2/5): motif'

echo 'Downloading Data(3/5): geneInfo'
wget https://sdsu.box.com/shared/static/mns0k1uvwtfnsohoc89b984ih36nmnz9.gz -O geneInfo.tar.gz
tar xvzf geneInfo.tar.gz
rm geneInfo.tar.gz
echo 'Download and Unzip Finished (3/5): geneInfo'

echo 'Downloading Data(4/5): data.go'
wget https://sdsu.box.com/shared/static/qwpdh36vcisgy1hcmadck8i8ezhvr2fh.gz -O data.go.tar.gz
tar xvzf data.go.tar.gz
rm data.go.tar.gz
echo 'Download and Unzip Finished (4/5): data.go'

echo 'Downloading Data(5/5): convertIDs.db'
wget https://sdsu.box.com/shared/static/sorewt7w6iypmhg2k2xhyi8myeit156o.gz -O convertIDs.db.tar.gz
tar xvzf convertIDs.db.tar.gz
rm convertIDs.db.tar.gz
echo 'Download and Unzip Finished (5/5): convertIDs.db'
echo ''
echo 'Data has been downloaded and unziped'
echo ''
echo ''
echo 'iDep is ready.'
echo ''
echo ''
