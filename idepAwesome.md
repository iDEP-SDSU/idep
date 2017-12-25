# idep awesome

docker save -o <save image to path> <image name>
docker load -i <path to image tar file>

# load balance 
https://www.sep.com/sep-blog/2017/02/28/load-balancing-with-nginx-and-docker/

docker-compose up -d --scale webapp=4

# remove all exited containers
sudo docker ps -a | grep Exit | cut -d ' ' -f 1 | xargs sudo docker rm

# remove all images
docker images -a | awk '{print $3}' | xargs docker rmi

# delete all 
docker system prune -a