#git pull

#git submodule update --init --recursive

#cd shinyapps/go
#git checkout master
#git pull

#cd ../idepgolem1
#git checkout master
#git pull
#cd ../..

# Restart docker dameon; sometimes it causes problem
sudo systemctl restart docker
#sudo docker system prune --volumes
sudo docker-compose down 
#free buffer/cache  https://unix.stackexchange.com/questions/87908/how-do-you-empty-the-buffers-and-cache-on-a-linux-system
sudo sh -c 'echo 1 >/proc/sys/vm/drop_caches'
sudo sh -c 'echo 2 >/proc/sys/vm/drop_caches'
sudo sh -c 'echo 3 >/proc/sys/vm/drop_caches'
sudo docker-compose up -d --scale webapp=25
