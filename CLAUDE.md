# CLAUDE.md

## What this repo is

Deployment code for the Ge lab web server at **bioinformatics.sdstate.edu**,
which hosts the Shiny apps [iDEP](https://github.com/gexijin/idepGolem) (RNA-seq
analysis) and [ShinyGO](https://github.com/gexijin/shinygo) (gene set
enrichment), plus every older version of both so that URLs cited in published
papers keep working. The app source lives in the repos above; this repo holds
the server layout, the Docker image, the reverse proxy, and the per-version app
directories.

Two hosting stacks exist side by side:

- **`master` branch (production today):** nginx in front of a fixed number of
  `shiny-server` containers (`docker compose up --scale webapp=50`). Config in
  the repo root: `docker-compose.yml`, `nginx/`, `config/shiny-server.conf`,
  `restart_server.sh`.
- **`shinyproxy` branch (under test):** ShinyProxy launches one container per
  user session, nginx routes the same URLs to it. Everything lives in
  `shinyproxy/`; see `shinyproxy/README.md` for architecture and operations.
  Do not run both stacks at once, both bind :80/:443.

Both stacks use the same `webapp` image built from the root `Dockerfile`
(rocker/shiny + R packages from `classes/librarySetup.R`) and the same
`shinyapps/` and `data/` mounts.

## Folder structure

```
Dockerfile, docker-compose.yml   webapp image and the master-branch stack
setup.sh                         fresh-server install: pull image, download databases
restart_server.sh                master-branch restart (compose down/up, 50 replicas)
classes/                         librarySetup.R installs R packages into the image
config/shiny-server.conf         URL -> app directory map for shiny-server (master)
nginx/                           nginx image + conf for the master stack; TLS cert/key untracked
shinyapps/                       one directory per app version, mounted at /srv/shiny-server
  idep, idep11 ... idep250       iDEP versions (some are git submodules)
  go, go41 ... go86              ShinyGO versions
  dist/                          static landing site served at /
  RSet, combine, reads, go2gmt   smaller helper apps
shinyproxy/                      ShinyProxy stack (this branch)
  SETUP.md                       one-time install checklist for a new host
  application.yml                app specs: current tier (/idep/, /go/) and pinned legacy tier
  nginx/nginx.conf               routes /<app>/ -> ShinyProxy on 127.0.0.1:8080
  idep.sh                        start/stop/status/check/pin/update, the only entry point
  shinyproxy.service.in          systemd unit template, rendered by ./idep.sh unit
  templates/                     ShinyProxy page templates (frame-escape for iframe links)
  bench*.sh                      cold-start benchmarks for both stacks
data/                            species databases (gitignored, ~100s of GB), /srv/data
countsData/, shinylog/, usage/   runtime data: uploads, shiny logs, telemetry sqlite
docs/                            server setup notes, Docker/Postgres how-tos
singularity_standalone/          Singularity build for running iDEP without Docker
```

## Conventions

- App URLs are versioned (`/idep250/`, `/go86/`); `/idep/` and `/go/` are
  aliases to the current version. Never change what an existing versioned URL
  serves.
- Legacy app versions must stay on a pinned image (`webapp:2026`); only the
  current versions track `webapp:latest`.
- Fix link/iframe issues in the app itself (`target="_blank"`) rather than in
  the proxy layer where possible.
- `data/`, `shinylog/`, `countsData/`, TLS keys, the ShinyProxy JAR and logs are
  gitignored; do not commit them.

## Routine app updates (shinyproxy branch)

`/idep/` runs an R package baked into the image; `/go/` and all legacy apps are
bind-mounted from `shinyapps/`.

```bash
# ShinyGO / legacy: no rebuild, picked up by the next container
cd shinyapps/go86 && git pull

# iDEP: rebuild the image (tags the outgoing one webapp:pre-<date>, then restarts)
cd shinyproxy && ./idep.sh update          # --pull to take gexijin/idep:latest instead
```

Never plain `docker build` — it overwrites `webapp:latest` with no snapshot, and
without a restart old containers keep serving the old build for days.

```bash
curl -s 127.0.0.1:9090/actuator/recyclable   # activeConnections: sessions a restart drops
./idep.sh status && ./idep.sh check idep250  # after updating
```
