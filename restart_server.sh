#!/bin/bash
set -euo pipefail

cd /docker/idep

docker compose down
docker compose up -d --scale webapp=40
