#!/bin/bash
# Times a fresh anonymous ShinyProxy session all the way to served app HTML.
# Uses SHINYPROXY_PUBLIC_PATH so it works with and without seat pre-init.
N=${1:-4}; GAP=${2:-20}; BASE=${BASE:-http://127.0.0.1:8080}; SPEC=${SPEC:-idep}
for run in $(seq 1 $N); do
  CJ=$(mktemp)
  curl -s -c $CJ -b $CJ -L -o /dev/null ${BASE:-http://127.0.0.1:8080}/
  START=$(date +%s.%N)
  RESP=$(curl -s -c $CJ -b $CJ -X POST -H 'Content-Type: application/json' -d '{}' ${BASE:-http://127.0.0.1:8080}/api/proxy/$SPEC)
  UP=$(date +%s.%N)
  ID=$(echo "$RESP"   | grep -oP '"id":"\K[^"]+' | head -1)
  PATH_=$(echo "$RESP"| grep -oP '"SHINYPROXY_PUBLIC_PATH":"\K[^"]+')
  if [ -z "$PATH_" ]; then echo "run $run FAILED: $RESP"; rm -f $CJ; continue; fi
  code=""
  for i in $(seq 1 1200); do
    code=$(curl -s -c $CJ -b $CJ -o /dev/null -w '%{http_code}' --max-time 5 "${BASE:-http://127.0.0.1:8080}${PATH_}")
    [ "$code" = "200" ] && break
    sleep 0.1
  done
  END=$(date +%s.%N)
  printf 'run %d: proxy ready %.2fs | app HTML (HTTP %s) %.2fs\n' "$run" "$(echo "$UP-$START"|bc)" "$code" "$(echo "$END-$START"|bc)"
  curl -s -c $CJ -b $CJ -X DELETE "${BASE:-http://127.0.0.1:8080}/api/proxy/$ID" >/dev/null
  rm -f $CJ
  sleep $GAP
done
