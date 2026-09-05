# ShinyProxy deployment for iDEP and ShinyGO

Alternative to the nginx + shiny-server stack in the repo root. ShinyProxy gives
each user session its own container started from the existing `webapp:latest`
image, running R directly — shiny-server is not involved.

| URL | spec | command in the container |
| --- | --- | --- |
| `/idep/` | `idep` | `idep250::run_app()` — R package installed in the image |
| `/go/` | `go` | `shiny::runApp("/srv/shiny-server/go86")` — mounted from `shinyapps/go86` |
| `/` | — | ShinyProxy's app list |

```
:80   nginx        container, --network host, config in nginx/nginx.conf
       |  /idep/ -> 127.0.0.1:8080/app/idep/
       |  /go/   -> 127.0.0.1:8080/app/go/
       |  /      -> 127.0.0.1:8080/            (JS, CSS, /api/, /app_proxy/)
       v
:8080 ShinyProxy   JAR on the host, config in application.yml,
       |           drives Docker via /var/run/docker.sock
       +-- webapp:latest containers, published on 127.0.0.1:2000x
```

## Operating it

```bash
./start.sh          # ShinyProxy on :8080 + nginx on :80
./stop.sh           # stops both, removes leftover app containers
./bench.sh 4        # time 4 cold sessions; BASE=http://127.0.0.1 SPEC=go ./bench.sh
./bench_shinyserver.sh   # same timing for the old shiny-server path, to compare
```

- **Nothing survives a reboot except nginx.** The nginx container has
  `--restart unless-stopped` but the ShinyProxy JAR does not, so after a reboot
  nginx answers :80 with 502 until you run `./start.sh`. Add a systemd unit if
  this becomes the production path.
- **Do not run the root `docker-compose.yml` at the same time** — both want :80.
- Logs: `startup.log` (console) and `shinyproxy.log` (structured). App container
  stdout via `docker logs sp-container-…`.

Health check:

```bash
docker ps --format '{{.Names}}' | grep -c '^sp-container'   # containers up (60 idle)
grep -c 'Created Seat' startup.log                          # seats created (108 idle)
curl -so /dev/null -w '%{http_code}\n' -L http://127.0.0.1/idep/
```

## After rebuilding `webapp:latest`

Running containers keep the old image, and ShinyProxy never cycles them on its
own. **`/admin/delegate-proxy` is not usable here** — `authentication: none`
means there are no admin users, so it returns 403. Restart instead:

```bash
docker pull gexijin/idep:latest && docker tag gexijin/idep webapp   # if pulling
./stop.sh && ./start.sh
```

The pool is back to 60 warm containers in well under a minute.

## Adding an app

Copy a block under `proxy.specs` in `application.yml`, then add two lines to
`nginx/nginx.conf`:

```yaml
    - id: reads
      display-name: ...
      container-image: webapp:latest
      container-cmd: [R, -e, 'shiny::runApp("/srv/shiny-server/reads", host = "0.0.0.0", port = 3838)']
      port: 3838
      container-volumes: [/home/exouser/idep/shinyapps/reads:/srv/shiny-server/reads:ro, ...]
      container-env: {IDEP_DATABASE: /srv/data}
      container-memory-limit: 8g
      minimum-seats-available: 4
```

```nginx
location = /reads { return 301 /reads/; }
location /reads/  { proxy_pass http://shinyproxy/app/reads/; }
```

Database paths must be passed explicitly — the apps' relative-path fallbacks
only resolve under shiny-server's working directory. iDEP reads `IDEP_DATA_DIR`,
ShinyGO reads `IDEP_DATABASE`; both point at `/srv/data`, mounted read-only
from `../data`.

## Configuration reference

Current values in `application.yml`, and what to change them for.

| | iDEP | ShinyGO | |
| --- | --- | --- | --- |
| `minimum-seats-available` | 48 | 60 | floor on **free** seats, not a total — the pool grows past it under load |
| `seats-per-container` | 1 | 5 | users sharing one R process |
| `allow-container-re-use` | `false` | *(unset)* | `false` = fresh R process per user; **only valid when `seats-per-container` is 1** |
| idle containers | 48 | 12 | **60 total, 108 seats, ~13.5 GB** |

Global:

- `max-total-instances: 108` — counts **seats (users), not containers**. Past
  it, new sessions get "not enough capacity"; existing ones are untouched.
- `container-memory-limit: 8g` — a *cap*, not a reservation: containers sit at
  224 MB and grow only as the user loads data. Exceeding it OOM-kills that
  container (exit 137) and nothing else.
- `heartbeat-timeout: 900000` — reclaim 15 min after the browser tab closes.
- `scale-down-delay: 5` — burst containers linger 5 min before the pool shrinks.
- `hide-navbar: true` — suppresses ShinyProxy's top bar.

Why iDEP is not shared: its sessions are memory-heavy and run long
single-threaded computations, so users on one R process would block each other
and share one OOM fate. ShinyGO is light — gene lists, a 5 MB upload cap,
sub-second start — so five per container is cheap. ShinyProxy spreads arrivals
across containers before packing them.

## Measured on this host (32 cores, 122 GB)

| | iDEP | ShinyGO |
| --- | --- | --- |
| warm seat available | 1.9 – 2.5 s | 0.3 – 0.4 s |
| cold container (pool empty) | 6.4 – 6.7 s | 3.8 s |
| shiny-server equivalent | 3.9 s cold / 1.9 s warm | — |

Idle container: 224 MB, 0.05 % CPU. 40 cold containers booted at once: all ready
in 14 s, peak load 8 of 32. The 8.9 GB database is a read-only bind mount, so it
sits in the host page cache once and is shared by every container.

**Capacity caveat.** 108 seats is an aggressive cap. Idle it is fine, but 108
*working* sessions at even 2 GB each would exceed 122 GB — the cap bounds the
count, not the sum of real usage. There is no measurement here of what a working
iDEP session costs, only idle (224 MB) and the 8 GB ceiling. Watch
`docker stats` under real traffic and set the cap to `~106 GB / observed p95`.

## Gotchas

- **Do not enable `track-app-url`.** `/idep/` is a real proxy pass, not a
  redirect, and works because every asset ShinyProxy emits is rooted at `/`.
  That one setting would rewrite the browser URL back to `/app/idep/`.
- ShinyProxy has no per-app URL setting (`target-path` is the path *inside* the
  container), which is why nginx does the mapping.
- `shinyproxy-3.2.4.jar` and `*.log` are gitignored; re-download the jar from
  <https://www.shinyproxy.io/downloads/>.
