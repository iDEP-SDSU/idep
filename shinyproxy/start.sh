#!/bin/bash
# Start the ShinyProxy test stack: ShinyProxy on :8080 + nginx front end on :80.
cd "$(dirname "$0")"
JAR=$(ls shinyproxy-*.jar | tail -1)

setsid nohup java -jar "$JAR" > startup.log 2>&1 < /dev/null &
disown
until curl -s -o /dev/null --max-time 2 http://127.0.0.1:8080/; do sleep 2; done
echo "ShinyProxy up on :8080 ($JAR, log: startup.log)"

docker rm -f sp-nginx >/dev/null 2>&1
docker run -d --name sp-nginx --network host --restart unless-stopped \
  -v "$PWD/nginx/nginx.conf:/etc/nginx/nginx.conf:ro" nginx:latest >/dev/null
echo "nginx up on :80"
echo "  http://$(hostname -I | awk '{print $1}')/idep/"
echo "  http://$(hostname -I | awk '{print $1}')/go/"
