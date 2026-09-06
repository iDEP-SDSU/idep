# ShinyProxy deployment for iDEP and ShinyGO

Replaces the nginx + shiny-server stack in the repo root. Every app that
`config/shiny-server.conf` used to publish is served here, at the same URL, so
existing bookmarks and links in published papers keep working.

Each user session gets its own container. shiny-server is not involved.

First-time installation on a new host: [SETUP.md](SETUP.md).

## The two tiers

| | current | legacy |
| --- | --- | --- |
| apps | `/idep/` `/go/` | the other 33 |
| image | `webapp:latest` | `webapp:2026` (pinned) |
| pool | pre-warmed, 20 + 20 free seats | none — cold start on first visit |
| idle cost | ~9 GB | **zero** |
| first load | 0.3 – 2.5 s | 1.6 – 9.6 s |

The legacy tier is pinned to a dated image tag **on purpose**. Rebuilding
`webapp:latest` every ~2 years moves R and every Bioconductor package under the
apps. Old iDEP and ShinyGO versions were never tested against those, and they
exist precisely so that published analyses stay reproducible — pinning is what
makes that true. `container-image` is a per-spec setting, which is the one
thing the old shiny-server layout could not express: there, all 30 apps shared
whatever R the single image happened to have.

```
:443  nginx        container, --network host, config in nginx/nginx.conf
:80    |           (:80 only redirects to https, except /data/)
       |  /            -> shinyapps/dist  (static, straight off disk)
       |  /data/       -> data/shared     (static, autoindex)
       |  /idep/ /go/  -> 127.0.0.1:8080/app/<id>/
       |  /idep73/ ... -> 127.0.0.1:8080/app/<id>/
       |  anything else with no matching file -> 127.0.0.1:8080
       v
:8080 ShinyProxy   JAR on the host, loopback only, config in application.yml,
       |           drives Docker via /var/run/docker.sock
       +-- webapp:latest / webapp:2026 containers on 127.0.0.1:2000x
```

Only :80 and :443 are reachable from outside. ShinyProxy (:8080) and its
actuator (:9090) bind 127.0.0.1, and every app container publishes its port on
127.0.0.1 too, so nginx is the only way in.

## Operating it

Everything goes through one script:

```bash
./idep.sh start           # ShinyProxy on :8080 + nginx on :80
./idep.sh stop            # stops both, sweeps leftover app containers
./idep.sh restart
./idep.sh status          # process / container / capacity summary
./idep.sh logs            # tail -f startup.log
./idep.sh check <id>      # really launch one app and time it to served HTML
./idep.sh check all       # smoke test every app (~2 min)
./idep.sh pin <tag>       # freeze webapp:latest as webapp:<tag>
./idep.sh update [--pull] # rebuild or pull the image, then restart
./idep.sh unit            # print the systemd unit for this checkout and user
```

`check` drives ShinyProxy's API the way a browser does and polls until R
actually serves HTML, then releases the seat. A plain `curl` of an app URL is
not a health check — ShinyProxy returns its loading page with HTTP 200 before
the container exists.

`parity` reads `config/shiny-server.conf` and asserts that every `location` it
publishes has an nginx route here, that the route resolves to a real spec, and
that the spec runs the same directory shiny-server ran. Run it after editing
either config, and on production before decommissioning shiny-server. Retire
the subcommand along with `shiny-server.conf` — once that file is gone,
`application.yml` is the only source of truth and there is nothing to compare
against.

For reboot survival and crash recovery install the unit — otherwise nginx
comes back by itself and answers :443 with 502 until someone runs
`./idep.sh start`:

```bash
./idep.sh unit | sudo tee /etc/systemd/system/shinyproxy.service >/dev/null
sudo systemctl daemon-reload
sudo systemctl enable --now shinyproxy
```

`unit` renders `shinyproxy.service.in` with the current user and the path of
this directory, because systemd cannot expand variables in `User=`,
`WorkingDirectory=` or `ExecStart=`. Re-run it if the checkout moves. The
unit runs the script via `/bin/bash` so that SELinux (enforcing on RHEL)
lets systemd execute it from a checkout that is not labeled `bin_t`.

The unit tracks the JVM through `shinyproxy.pid` and restarts it if it
crashes or is OOM-killed; `start` sweeps the containers a dead instance left
behind. Once installed, prefer `systemctl start|stop|restart shinyproxy` over
calling the script directly so systemd's view stays right.

## TLS

nginx terminates TLS with the certificate production already has. `start`
mounts it from `/etc/pki/nginx/server.pem` and
`/etc/pki/nginx/private/server.key` — the paths the previous nginx image
used internally — when **both** exist on the host, and otherwise from the
untracked `../nginx/idep_ssl.pem` and `../nginx/idep.key`, the pair the
previous image was built from (gitignored, key mode 600). `start` prints
which one it used. Either path can be overridden:

```bash
CERT_PEM=/path/to/fullchain.pem CERT_KEY=/path/to/privkey.key ./idep.sh start
# or Environment=CERT_PEM=... CERT_KEY=... in the [Service] section of shinyproxy.service.in
```

The `.pem` must contain the server certificate followed by any intermediate,
as nginx expects. `start` refuses to run without both files, and validates
`nginx.conf` with `nginx -t` before replacing the running nginx.

The previous production nginx sent `Strict-Transport-Security` with a one-year
max-age, so browsers that have visited the site refuse plain http. This nginx
sends the same header and redirects everything on :80 to https — except
`/data/`, which is served on both ports because the apps download missing
species databases from `http://<host>/data/` at runtime.

`server.forward-headers-strategy: native` makes ShinyProxy honour the
`X-Forwarded-Proto` nginx sets, so its redirects and cookies use https.

## Rate limiting app starts

There is no authentication, so a client that discards its cookie is a new user
on every request and could start containers until `max-total-instances` is
exhausted. nginx limits the three ways to start one — `POST /app_i/…` (the app
page), `POST /api/proxy/<spec>` (the REST API), and `GET /app_direct*/…` — to
20 per minute per client IP with a burst of 40: the first 20 go through at
once, the next 20 are queued, anything beyond that gets 429. A classroom
behind one NAT address starting together fits inside the burst; a script
hammering the API does not. Status polling and `/api/proxyspec` are GETs on
the same prefix and are not counted.

Benchmarks: `./bench.sh 4` times cold sessions (`SPEC=go74 ./bench.sh`);
`./bench_shinyserver.sh` does the same against the old stack, to compare.

**Do not run the root `docker-compose.yml` at the same time** — both want :80.

## Updating the image

Roughly every two years, with updated R and packages:

```bash
./idep.sh update            # build from the repo Dockerfile
./idep.sh update --pull     # or pull gexijin/idep:latest
```

`update` tags the outgoing image `webapp:pre-<date>` before replacing it, then
restarts — running containers keep the old image otherwise, and ShinyProxy
never cycles them on its own. When `shinyproxy.service` is active the restart
goes through `sudo systemctl restart shinyproxy`, so systemd keeps tracking
the new JVM. (`/admin/delegate-proxy` is not usable here:
`authentication: none` means there are no admin users, so it returns 403.)

Only `/idep/` and `/go/` move to the new image. The legacy specs stay on
`webapp:2026` and **should never be repointed at `:latest`** — that is the whole
point of the tier. When you eventually want to freeze a newer vintage for a new
generation of legacy apps:

```bash
./idep.sh pin 2028          # webapp:latest -> webapp:2028
```

then set `container-image: webapp:2028` on the specs that should hold there.
Images are ~30 GB on disk; keeping one per vintage is cheap next to the 1.4 TB
free here.

## Moving /idep/ and /go/ to a new version

Spec ids carry the version: `idep250`, `go86`, `idep210`, ... . `/idep/` and
`/go/` are nginx aliases onto the current spec, so a release never edits an
existing spec — the old version keeps its URL and its behaviour. To ship
iDEP 2.60 (ShinyGO is the same with `go`):

1. Build the new image (`./idep.sh update`) with the `idep260` package in it.
2. Copy the `idep250` block in `application.yml` to a new spec `idep260`,
   changing the id, display name and the `idep260::run_app` call. Give it
   the `minimum-seats-available` pool.
3. Demote `idep250`: drop `minimum-seats-available` and
   `allow-container-re-use`, set `seats-per-container: 5`, and point
   `container-image` at a tag that still has the `idep250` package
   (`./idep.sh pin <tag>` before the update creates one). It now costs
   nothing until someone visits.
4. In `nginx/nginx.conf`, add `location = /idep260` and `location /idep260/`
   lines like the ones for `/idep250/`, then repoint `location /idep/` at
   `/app/idep260/`.
5. `./idep.sh restart`, then `./idep.sh check idep260` and `check idep250`.

The `/idep/` route is a proxy alias, not a redirect, so the address bar keeps
showing `/idep/` and bookmarks keep following the latest version. Change the
line to `return 301 /idep250/;` if you would rather the version show.

## Adding an app

Three lines in `application.yml` — everything else is inherited from the
`&legacy` anchor by YAML merge key:

```yaml
    - <<: *legacy
      id: go90
      display-name: ShinyGO 0.90
      container-cmd: [R, -e, 'shiny::runApp("/srv/shiny-server/go90", host = "0.0.0.0", port = 3838)']
```

and two in `nginx/nginx.conf`:

```nginx
location = /go90 { return 301 /go90/; }
location /go90/  { proxy_pass http://shinyproxy/app/go90/; }
```

`container-cmd` is spelled out per app rather than inherited, so the file says
plainly which directory each spec runs.

## How the legacy apps find their databases

They hardcode relative paths — `"../../data/data104b/"`, `"../go/geneInfo/"`.
`shiny::runApp("/srv/shiny-server/idep73")` sets the working directory to the
app directory, exactly as shiny-server did, so with the whole `shinyapps` tree
mounted at `/srv/shiny-server` and the databases at `/srv/data`, every one of
those paths resolves to the same file it did before. That is why the legacy
anchor mounts the whole tree instead of a single app directory.

The golem-packaged versions read an environment variable instead and append
their own database version to it. Both names are set on every legacy spec —
`IDEP_DATABASE` (ShinyGO, `idepGolem*`) and `IDEP_DATA_DIR` (`idep250`) — which
is harmless where unused:

| app | package | database |
| --- | --- | --- |
| `/idep/` | `idep250` 2.5.0 | `data115` |
| `/idep210/` | `idepGolem` 2.4.4 | `data113` |
| `/idep20/` | `idepGolem201` 2.0 | `data107` |
| `/go/` | `go86` source | `data115` |

The tree is mounted **read-only**, which shiny-server's was not. An app that
writes into its own directory will fail — drop the `:ro` if one does.

## Status on this test host

`./idep.sh parity` passes: all 37 `location` entries in the old
`shiny-server.conf` are published here, on the same directories. That is 35
specs — `/idep` and `/idep250` both route to spec `idep250`, `/go` and `/go86`
to `go86`.

`./idep.sh check all` passes 22 of 35. Every failure is a missing input on
*this* machine, not configuration:

- `idep11` `idepg` `go80` `go77` — submodules not checked out here
  (`git submodule update --init`)
- `go85` `go82` `go81` `datamap` — directories not present on this host
- `go74` `go75` `go76` — need `data104b`
- `go65` — needs `data103`
- `goc` — needs `customDB`

Only `data96`, `data113` and `data115` are present here. Apps differ in when
they touch the database: `go74` opens it in `global.R` and so fails to start at
all, while `go41` and `idep73` start fine and only fail once a user picks a
species — so a green check is not proof the database underneath is right.

**On production, where all of these exist, `./idep.sh check all` should be 35 of
35.** Anything still failing there is a real problem.

## Links between apps

ShinyProxy shows every app inside an `<iframe>` on its own page. A link from
one app to another that has no `target` -- iDEP's "old versions are still
usable", the legacy apps' "try the new version" -- navigates that iframe, so
the second app loads nested inside the first one's page: the address bar
keeps the first URL, the first container keeps its seat because its heartbeat
never stops, and a reload brings the first app back.

Two fixes, both in use. The apps themselves should put `target="_blank"` on
such links; that is the real fix and needs nothing from this directory. As a
safety net for links that lack it, `templates/frame-escape.html` is a few
lines of JavaScript that make a framed app page take over the tab.
`./idep.sh start` inserts it into the stock `app.html` from the JAR and writes
the result to `templates/app.html`, which `proxy.template-path` points
ShinyProxy at. The generated file is gitignored and rebuilt on every start, so
a ShinyProxy upgrade cannot leave a stale copy behind; the snippet is the only
thing to maintain. If an upgrade renames the `<head lang="en">` anchor the
insertion targets, `start` refuses to run rather than silently serving the
stock page.

It works only when both apps are on the same origin. Chrome blocks a
cross-origin frame from navigating the top window unless that document itself
has a user gesture, and the click happened in the previous document, so the
attempt is refused with "Unsafe attempt to initiate navigation ... nor has it
received a user gesture". The links in the apps are absolute URLs to
production, so on this test host they leave the stack, show production inside
the frame, and the escape does nothing; it only takes effect once production
runs this stack and the links become same-origin. For the same reason it does
not stop a third-party site from embedding an app.

## Configuration reference

| | iDEP | ShinyGO | legacy |
| --- | --- | --- | --- |
| `minimum-seats-available` | 20 | 20 | unset — no pre-init |
| `seats-per-container` | 1 | 1 | 5 |
| `allow-container-re-use` | `false` | `false` | *(unset)* |
| idle containers | 20 | 20 | 0 |

`minimum-seats-available` is a floor on **free** seats, not a total — every
seat a user takes is replaced at once, so the pool grows past it under load
and 20 spares serve any number of users as long as fewer than 20 arrive
within one cold start (~7 s). Setting it at all is what enables
pre-initialization, which is why omitting it on the legacy tier makes those
apps cost nothing.

`allow-container-re-use: false` gives every user a fresh R process and is
**only valid when `seats-per-container` is 1**.

Global:

- `max-total-instances: 200` — counts **seats (users), not containers**,
  across every spec. Since the pools grow without bound under load, this is
  the only cap on concurrent sessions. Past it, new sessions get "not enough
  capacity"; existing ones are untouched.
- `container-memory-limit: 15g` — a *cap*, not a reservation: containers sit at
  224 MB and grow only as the user loads data. Exceeding it OOM-kills that
  container (exit 137) and nothing else.
- `heartbeat-timeout: 900000` — reclaim 15 min after the browser tab closes.
- `scale-down-delay: 5` — burst containers linger 5 min before the pool shrinks.
- `hide-navbar: true` — suppresses ShinyProxy's top bar.

Why the current apps are not shared: an iDEP session is memory-heavy and runs
long single-threaded computations, so users on one R process would block each
other and share one OOM fate. ShinyGO is light enough that a fresh R process
per user costs little, and it keeps the two specs identical. The legacy tier
still packs five users per container, as shiny-server did.

## Measured on this host (32 cores, 122 GB)

| | iDEP | ShinyGO | legacy (cold) |
| --- | --- | --- | --- |
| warm seat available | 1.9 – 2.5 s | 0.3 – 0.4 s | — |
| cold container (pool empty) | 6.4 – 6.7 s | 3.8 s | 1.6 – 9.6 s |
| shiny-server equivalent | 3.9 s cold / 1.9 s warm | — | ~3 s |

Idle container: 224 MB, 0.05 % CPU. 40 cold containers booted at once: all ready
in 14 s, peak load 8 of 32. The database is a read-only bind mount, so it sits
in the host page cache once and is shared by every container.

**Capacity caveat.** 200 seats is an aggressive cap. Idle it is fine, but 200
*working* sessions at even 1 GB each would exceed 122 GB — the cap bounds the
count, not the sum of real usage. There is no measurement here of what a working
iDEP session costs, only idle (224 MB) and the 15 GB ceiling. The total is
bounded by `memory-cap.sh` instead, below; watch `docker stats` under real
traffic to see how far from it a normal day sits.

## Capping total container memory

ShinyProxy limits memory per container only (`container-memory-limit: 15g`)
and has no setting for the total, so the seat cap alone leaves the host
exposed to `200 × 15 GB`. `memory-cap.sh` closes that gap one level down:

```bash
sudo ./memory-cap.sh install          # once per host; budget is MEMORY_MAX at the top of the script
sudo ./memory-cap.sh status
sudo ./memory-cap.sh remove
```

`install` writes a systemd slice, `containers.slice`, with that `MemoryMax`,
and sets `cgroup-parent` in `/etc/docker/daemon.json` so Docker starts every
container inside it. When the containers' combined usage reaches the budget
the kernel first drops the page cache charged to the slice (the database
bind mount, mostly), then OOM-kills the largest process in it — one heavy R
session, which its user sees as "This app has crashed" — and everything
outside the slice (the OS, Docker, the ShinyProxy JVM, this shell) is never
touched. Verified here: a slice held to 200 MB kills an R process the moment
it allocates 320 MB, exit 137, `oom_kill` counted in the slice's
`memory.events`.

`install` runs **once per host**. Both files are ordinary persistent config:
Docker re-reads `daemon.json` on every start, and systemd applies the slice's
`MemoryMax` whenever the slice comes up, so the cap survives `systemctl
restart docker` and reboots without any help. `install` is idempotent and only
restarts Docker — which restarts every container, then ShinyProxy to rebuild
the pool — when `daemon.json` actually changed. The budget applies to *all*
containers on the host, ShinyProxy's or not; nginx (~70 MB) and the Guacamole
pair here share it. Containers created before the cap stay outside it until
recreated; `status` lists them.

The budget is `MEMORY_MAX` at the top of the script, 140G, sized for the
production server; keep at least 16 GiB of host RAM outside it for the OS, the
page cache and the ShinyProxy JVM (2 GB heap cap, ~1.3 GB resident after a
busy hour). On a host with less RAM than the budget, such as this 122 GB test
server, the slice never binds and only the per-container limit protects the
host. To change it, edit that line and run
`install` again; the new value is applied to the running slice in place, no
Docker restart.

## Gotchas

- **Do not enable `track-app-url`.** `/idep/` is a real proxy pass, not a
  redirect, and works because every asset ShinyProxy emits is rooted at `/`.
  That one setting would rewrite the browser URL back to `/app/idep250/`.
- `templates/app.html` is generated; edit `templates/frame-escape.html` (see
  "Links between apps").
- ShinyProxy has no per-app URL setting (`target-path` is the path *inside* the
  container), which is why nginx does the mapping.
- nginx's `try_files ... @shinyproxy` fallback is what lets the static site and
  ShinyProxy share the root. ShinyProxy's own `/js/`, `/css/`, `/webjars/` and
  `/app_proxy/` live inside the JAR, never on disk, so they never collide with
  a file in `dist/` and never need to be listed.
- Nothing is tied to a checkout path or user. `application.yml` writes its
  bind mounts as `${IDEP_ROOT}/...`, which Spring resolves from the
  environment `idep.sh` exports, and the systemd unit is rendered by
  `./idep.sh unit`. `idep.sh` itself locates the repo relative to its own
  path.
- The nginx image is pinned (`NGINX_IMAGE` in `idep.sh`); `nginx:latest` is
  mainline and would move under us on every start. Pre-pull the tag on a new
  host so a reboot does not depend on Docker Hub.
- The JVM runs with `-Xmx2g` (`JAVA_OPTS` in `idep.sh`); without a cap it
  would grow the heap to a quarter of host RAM before collecting.
- `shinyproxy-3.2.4.jar`, `*.log`, `shinyproxy.pid` and the TLS files in
  `../nginx/` are gitignored; re-download the jar from
  <https://www.shinyproxy.io/downloads/>.
