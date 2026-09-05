#!/bin/bash
# Stop the ShinyProxy test stack and remove the app containers it left behind.
docker rm -f sp-nginx >/dev/null 2>&1 && echo "nginx stopped"
PID=$(ss -tlnpH 'sport = :8080' 2>/dev/null | grep -oP 'pid=\K[0-9]+' | head -1)
[ -n "$PID" ] && kill "$PID" && echo "ShinyProxy stopped (pid $PID)"
sleep 3
for c in $(docker ps -aq --filter "label=openanalytics.eu/sp-proxy-id"); do docker rm -f "$c" >/dev/null; done
echo "app containers removed"
