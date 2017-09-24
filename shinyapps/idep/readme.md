# idep Shipy app 

## Build docker image
docker build -t image_name .


mv *.* ..

## Run bash
docker exec -it idep_web_1 bash

## Remove All exit container
sudo docker ps -a | grep Exit | cut -d ' ' -f 1 | xargs sudo docker rm
https://stackoverflow.com/questions/24482822/how-to-share-my-docker-image-without-using-the-docker-hub

### How to delete dangling volume 
https://stackoverflow.com/questions/30604846/docker-error-no-space-left-on-device

### known issue 
[issue1] rscript doesn't execute when docker compose up
[issue2] during install dplyr: [error] virtual memory exhausted: Cannot allocate memory

+ means your vm doesn't have sufficient memory to compile it. (Need more Memory)

R -e "shiny::runApp('./shinyapps/idep')"