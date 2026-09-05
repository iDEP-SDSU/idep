#!/bin/bash
# Cold-start the app the way the current production stack does it:
# a webapp container running shiny-server, first request to /idep250.
PORT=3940
docker rm -f idep_ss_bench >/dev/null 2>&1
docker run -d --rm --name idep_ss_bench -p 127.0.0.1:$PORT:3838 \
  -v /home/exouser/idep/shinyapps/:/srv/shiny-server/ \
  -v /home/exouser/idep/data/:/srv/data/:ro \
  -v /home/exouser/idep/countsData/:/srv/countsData/ \
  -v /home/exouser/idep/shinylog/:/var/log/shiny-server/ \
  -v /home/exouser/idep/config/:/etc/shiny-server/ \
  -v /home/exouser/idep/classes/:/usr/local/src/myscripts/ \
  -v /home/exouser/idep/usage/:/srv/usage/ \
  -e IDEP_TELEMETRY_SQLITE=../../usage/idep_telemetry.sqlite \
  webapp:latest >/dev/null

# wait for shiny-server itself to accept connections (not the app)
until [ "$(curl -s -o /dev/null -w '%{http_code}' --max-time 2 http://127.0.0.1:$PORT/ 2>/dev/null)" != "000" ]; do sleep 0.2; done
echo "shiny-server listening"

# FIRST user hits /idep250 -> shiny-server must fork a fresh R process
START=$(date +%s.%N)
code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 300 http://127.0.0.1:$PORT/idep250/)
END=$(date +%s.%N)
printf 'cold  (first user, R process spawned): HTTP %s after %.2fs\n' "$code" "$(echo "$END-$START"|bc)"

# SECOND user reuses the already-running R process
S2=$(date +%s.%N)
code2=$(curl -s -o /dev/null -w '%{http_code}' --max-time 300 http://127.0.0.1:$PORT/idep250/)
E2=$(date +%s.%N)
printf 'warm  (second user, process reused) : HTTP %s after %.2fs\n' "$code2" "$(echo "$E2-$S2"|bc)"
docker rm -f idep_ss_bench >/dev/null 2>&1
