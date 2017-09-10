# idep Shipy app 

wget -O geneInfo.zip https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb
curl -sL --retry 3 \
   \
  | unzip \
  | tar -x -C ./geneInfo

wget -qO- https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb | tar --transform './geneInfo' -xvz

ADD https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb /tmp/
RUN tar -xjf /srv/shiny-server/idep/data/geneInfo \
  && rm /tmp/*.zip


docker build -t image_name .
docker build --build-arg CHECKPPOINT=$(Copy) -t image_name .


docker run -d --name container_name image_name
## Remove All exit container
sudo docker ps -a | grep Exit | cut -d ' ' -f 1 | xargs sudo docker rm


## Run bash
docker run -it image_name /bin/bash



https://stackoverflow.com/questions/24482822/how-to-share-my-docker-image-without-using-the-docker-hub