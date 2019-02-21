# Based on digital ocean One click app: Docker-ce

docker build ./nginx/. -t nginx  #nginx image should be build very quick

docker build . -t webapp #webapp image need hours to build

echo 'Docker images have been built. Start downloading data.'

mkdir data
cd data

wget https://sdsu.box.com/shared/static/c24f792ojoikpzu0lkpng8uuf9ychwm7.gz -O pathwayDB.tar.gz
tar xvzf pathwayDB.tar.gz
rm pathwayDB.tar.gz
wget https://sdsu.box.com/shared/static/9v1ao6mwhduvrcx793j3answph9gqnkt.gz -O motif.tar.gz
tar xvzf motif.tar.gz
rm motif.tar.gz
wget https://sdsu.box.com/shared/static/mns0k1uvwtfnsohoc89b984ih36nmnz9.gz -O geneInfo.tar.gz
tar xvzf geneInfo.tar.gz
rm geneInfo.tar.gz
wget https://sdsu.box.com/shared/static/qwpdh36vcisgy1hcmadck8i8ezhvr2fh.gz -O data.go.tar.gz
tar xvzf data.go.tar.gz
rm data.go.tar.gz
wget https://sdsu.box.com/shared/static/sorewt7w6iypmhg2k2xhyi8myeit156o.gz -O convertIDs.db.tar.gz
tar xvzf convertIDs.db.tar.gz
rm convertIDs.db.tar.gz

echo 'Data has been downloaded and unziped'

echo 'All image are ready to run'

