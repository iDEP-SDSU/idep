# Initial setup on a new host

One-time steps to bring the ShinyProxy stack up on a server that currently
runs the master-branch stack (nginx + shiny-server via `docker-compose.yml`).
Day-to-day operation is in [README.md](README.md).

Nothing in the stack is tied to a checkout path or login user: `idep.sh`
finds the repo relative to its own location, `application.yml` mounts
`${IDEP_ROOT}/...` which the script exports, and the systemd unit is rendered
for the current user and directory. The examples below use `/docker/idep` and
the user `gex`; substitute your own.

## 1. Host prerequisites

The script starts ShinyProxy as a plain JAR on the host, so the host needs:

| what | why | check |
| --- | --- | --- |
| Java 17 or newer | ShinyProxy 3.2.x is a Spring Boot 3 JAR and does not start on older runtimes | `java -version` (RHEL 9's default is often 11; see below) |
| `unzip`, `bc`, `curl`, `ss`, `python3` | `start` extracts `app.html` from the JAR with unzip; `check` and `memory-cap.sh status` use bc; every probe uses curl; `ss` finds a ShinyProxy whose pidfile is gone; `memory-cap.sh` edits `daemon.json` with python3 | the loop below prints nothing missing |
| Docker on cgroup v2 with the systemd driver | `memory-cap.sh` puts every container in one systemd slice and refuses otherwise | `docker info -f '{{.CgroupDriver}} {{.CgroupVersion}}'` prints `systemd 2` |
| the login user in the `docker` group | ShinyProxy talks to `/var/run/docker.sock`; the old stack was run with `sudo`, so this may not be set up yet | `docker ps` without sudo |
| `nginx:1.31-alpine` image pulled | pinned in `idep.sh`; the same base the old nginx image was built from, so production has it already. Pre-pull anyway so a reboot does not depend on Docker Hub | `docker pull nginx:1.31-alpine` |

Production is RHEL 9, which ships cgroup v2 and Docker CE defaults to the
systemd driver there, so the third row should already hold:

```bash
sudo dnf install -y java-21-openjdk-headless unzip bc iproute python3 curl
for t in java unzip bc curl ss python3; do command -v $t >/dev/null || echo "MISSING $t"; done
sudo usermod -aG docker "$USER"   # then log out and back in
docker ps                          # must work without sudo
docker info -f '{{.CgroupDriver}} {{.CgroupVersion}}'
docker pull nginx:1.31-alpine
getenforce                         # Enforcing on a stock RHEL host; see step 9
```

(On Debian/Ubuntu: `apt-get install -y openjdk-21-jre-headless unzip bc iproute2 python3 curl`.)

If `java -version` still reports 11 after the install, the system default
was not changed (`dnf` leaves an existing alternative in place). Either
leave it, in case something else on the host wants Java 11, and point the
stack at 21 explicitly, or switch the default:

```bash
ls /usr/lib/jvm                                   # jre-21, jre-21-openjdk, ...
export JAVA=/usr/lib/jvm/jre-21/bin/java          # for every ./idep.sh call in this shell
$JAVA -version                                    # 21.x
# or instead: sudo alternatives --config java      # makes 21 the default for everyone
```

With `JAVA` exported, `./idep.sh unit` (step 9) writes it into the systemd
unit, so run `unit` from the same shell; the rehearsal (step 8) uses it too. `start` checks the version and
refuses anything older than 17.

`idep.sh start` refuses to run without `java` or `unzip` on the PATH.

## 2. Check out the branch and its submodules

Production has files that were edited in place rather than committed (the
root `nginx/nginx.conf` almost certainly, since master commits 40 upstreams
and runs 50; possibly `config/shiny-server.conf`). `git checkout` refuses to
switch branches over a modified tracked file, and a `git stash` would
conflict on `nginx/nginx.conf`, which this branch rewrites. Commit them on
`master` instead; that also means a rollback (`git checkout master`) brings
production's own files back:

```bash
cd /docker/idep
git status                        # tracked files with local edits?
git commit -am "Production in-place edits"   # if any
git fetch origin
git checkout shinyproxy
git submodule update --init       # idep11, go80, idepgolem1, go
```

The new stack does not read the root `nginx/` or `docker-compose.yml`, but
`./idep.sh parity` (step 8) does read `config/shiny-server.conf`. If
production's copy differs from the committed one, put it back for the
comparison so parity checks what production really served:

```bash
git diff master -- config/shiny-server.conf          # anything?
git checkout master -- config/shiny-server.conf      # then: git checkout shinyproxy -- config/shiny-server.conf
```

Confirm every app directory that `config/shiny-server.conf` publishes is
present under `shinyapps/` (some, such as `go86`, were copied onto the server
rather than committed):

```bash
grep -oP 'app_dir /srv/shiny-server/\K\S+(?=;)' config/shiny-server.conf | sort -u | while read d; do
  [ -d "shinyapps/$d" ] || echo "MISSING shinyapps/$d"
done
```

## 3. Download the ShinyProxy JAR

It is gitignored. Get **3.2.4**, the version the frame-escape template was
built against, from <https://www.shinyproxy.io/downloads/> and put it in
this directory:

```bash
cd /docker/idep/shinyproxy
curl -LO https://www.shinyproxy.io/downloads/shinyproxy-3.2.4.jar
ls -l shinyproxy-*.jar           # ~150 MB
```

`idep.sh` picks the highest-versioned `shinyproxy-*.jar` it finds here.

## 4. Run as the user who owns the checkout

The old stack was driven with `sudo docker compose`, so everything it did
was done by root. The new stack is different: ShinyProxy is not a container
but a Java process that `idep.sh` starts, and the script writes into this
directory (`shinyproxy.pid`, `startup.log`, `shinyproxy.log`,
`templates/app.html`). Whoever runs it owns those files.

So pick one ordinary login account (`gex` in these examples) and use it for
**every** `./idep.sh` command, never `sudo ./idep.sh`. Three things depend
on that being the same account throughout:

- `./idep.sh unit` (step 9) writes `User=<whoever ran it>` into the
  systemd unit. After that the service runs as that account and has to be
  able to overwrite the files above; a `shinyproxy.pid` left behind by a
  root run makes the unit fail to start.
- The account must be able to read the checkout, including the private key
  under `nginx/` (mode 600), and to run `git checkout` in it (step 2), which
  means it should own the checkout.
- It must be in the `docker` group (step 1). That is all the privilege the
  stack needs: Docker's daemon, which is root, binds :80 and :443 for the
  nginx container, and ShinyProxy itself listens on the unprivileged :8080.

The account does **not** need write access to `usage/`, `countsData/` or
`data/`: the app containers run as root inside and write through the bind
mounts as root, exactly as before.

Check the ownership before going on. If the checkout was cloned with sudo
and is root-owned, hand it to the account (`data/` is large but `chown` only
touches inodes, so this takes seconds to a minute):

```bash
ls -ld /docker/idep /docker/idep/.git         # owner should be gex, not root
sudo chown -R gex:gex /docker/idep             # only if it is root-owned
```

`sudo` is still used for the host-level steps: installing packages, the
memory cap, and the systemd unit.

## 5. Pin the legacy image

Every legacy spec runs `webapp:2026`. Create that tag from the current image
**before the first start and before any `./idep.sh update`**, because `update`
replaces `webapp:latest`:

```bash
docker image ls webapp            # confirm webapp:latest exists
./idep.sh pin 2026                # webapp:latest -> webapp:2026
```

Without it `start` stops with `image webapp:2026 not found`.

## 6. Make sure the mounted directories exist

The specs bind-mount these from the repo root. `data/`, `countsData/` and
`shinyapps/` already exist on a server that ran the old stack. `usage/` holds
the telemetry SQLite database and is not in git:

```bash
mkdir -p /docker/idep/usage
```

If it is missing Docker creates it as root; that also works, since the app
containers run as root, but the directory then cannot be removed without sudo.

## 7. TLS

`start` looks for the certificate in two places and takes the host pair
when **both** of its files exist, otherwise the checkout pair:

1. `/etc/pki/nginx/server.pem` + `/etc/pki/nginx/private/server.key` on the
   host (the RHEL convention).
2. `nginx/idep_ssl.pem` + `nginx/idep.key` in the checkout, untracked.

The old stack **copied** the pair from `nginx/` into its nginx image at the
`/etc/pki/nginx` paths (see `nginx/Dockerfile`), so the checkout copy is the
one production definitely has; whether the host also has a copy under
`/etc/pki` is unknown. Check both, and if both exist make sure the host one
is not an older certificate, because it wins silently:

```bash
cd /docker/idep
ls -l /etc/pki/nginx/server.pem /etc/pki/nginx/private/server.key nginx/idep_ssl.pem nginx/idep.key
openssl x509 -enddate -subject -noout -in /etc/pki/nginx/server.pem   # if it exists
openssl x509 -enddate -subject -noout -in nginx/idep_ssl.pem
```

`start` prints which pair it used. To force a particular one set `CERT_PEM`
and `CERT_KEY` (see README, "TLS"). `start` refuses to run without both
files. The `.pem` must be the
server certificate followed by any intermediate, as before.

## 8. Rehearse ShinyProxy alone, with the old stack still serving

ShinyProxy listens on 127.0.0.1:8080, which the old stack does not use, so
it can run on production before anything is stopped. This starts the JAR
without its nginx, runs every app against the real directories and
databases, and stops it again. It is the only step that exercises all 35
apps before the site depends on them; ten minutes, no disruption. The
pre-warmed pools add about 9 GB alongside the old stack while it runs.

```bash
cd /docker/idep/shinyproxy
IDEP_ROOT=/docker/idep $JAVA -Xmx2g -jar shinyproxy-3.2.4.jar > startup.log 2>&1 &
until curl -s -o /dev/null http://127.0.0.1:8080/; do sleep 2; done
./idep.sh parity                  # every shiny-server location has a route and a spec
./idep.sh check all               # ~2 min; must be 35 of 35 on production
./idep.sh stop                    # stops the JAR and sweeps its containers ("nginx stopped" does not print)
```

`parity` compares against `config/shiny-server.conf`; anything it reports
means a URL from the old server would now 404. Any `check` failure on
production is a real problem, since all data directories and app versions
exist there (the test host fails 13 only because they do not). Fix it now,
while the old stack is still serving.

## 9. Install and enable the systemd unit, without starting it

Installing it before the cutover means the first start already goes
through systemd, instead of a manual start followed by a stop and a second
start. `enable` only registers the service for boot; nothing runs until
`systemctl start` in step 11.

```bash
./idep.sh unit                    # review: User=, WorkingDirectory=, ExecStart=, Environment=JAVA= must match this host
./idep.sh unit | sudo tee /etc/systemd/system/shinyproxy.service >/dev/null
sudo systemctl daemon-reload
sudo systemctl enable shinyproxy
```

The unit runs `idep.sh` through `/bin/bash` on purpose. With SELinux
enforcing (the RHEL default), systemd may only execute files whose label is
an executable type, and a checkout under `/docker` or a home directory is
`default_t` or `user_home_t`; executing the script directly fails with
`status=203/EXEC` and "Permission denied" in the journal even though
`./idep.sh start` works from a shell. Executing `/bin/bash` and handing it
the script sidesteps that.

## 10. Disable whatever brings the old stack back

Both stacks bind :80 and :443, so nothing may restart the old one after the
cutover: a cron entry or a systemd unit that runs `restart_server.sh` or
`docker compose up`, or containers with `restart: always` left from earlier
compose versions.

```bash
sudo crontab -l; crontab -l
systemctl list-units --type=service | grep -i 'idep\|compose\|shiny'
docker ps -a --filter name=idep --format '{{.Names}} {{.Status}}'
```

## 11. Cutover

Pick a quiet hour; sessions on the old stack end when it stops. The site is
down from `compose down` until `systemctl start` returns, about a minute.
The memory cap goes in between because it restarts Docker, which is cheap
while nothing is running.

```bash
docker stats --no-stream | sort -k3 -h | tail    # any busy sessions?
cd /docker/idep && sudo docker compose down      # site down from here
docker ps -a --filter name=idep                  # should list nothing
cd /docker/idep/shinyproxy
sudo ./memory-cap.sh install                     # budget is MEMORY_MAX at the top of the script (140G)
sudo systemctl start shinyproxy                  # site up when this returns
./idep.sh status
```

Keep at least 16 GiB of host RAM outside the memory budget. On a host with
less RAM than `MEMORY_MAX` the slice never binds; edit the value first.

If the start fails, get the site up first and debug second. The unit is not
active after a failed start, so the script allows a manual start:

```bash
./idep.sh start
sudo journalctl -u shinyproxy -n 30
sudo ausearch -m avc -ts recent           # SELinux denials, if any
```

Once the unit is in charge use `systemctl start|stop|restart shinyproxy`
rather than `./idep.sh` directly, so systemd's view stays right;
`./idep.sh update` does this by itself when the unit is active.

## 12. Verify, live

```bash
./idep.sh check all               # 35 of 35 again, now through the unit
sudo ./memory-cap.sh status
```

Then from a browser, over https:

- `/` landing page, `/data/` listing (also over plain http).
- `/idep/` and `/go/`, and a couple of old versions such as `/idep73/`,
  `/go41/`.
- Upload the example data in `/idep/` and run through to a plot, so the
  database and `countsData` mounts are exercised, not only the app start.
- Follow the "old versions" link inside iDEP; it should open in a new tab or
  take over the tab, never nest inside the frame.

Reboot once at a convenient time and confirm the site answers without
manual intervention.

## 13. Afterwards

- Watch `docker stats` and `./idep.sh status` during the first busy day to
  see how far real traffic sits from the memory budget and the seat cap.
- Once the old stack is gone for good, `config/shiny-server.conf` and the
  `parity` subcommand can be retired; `application.yml` becomes the only
  source of truth for which apps exist.

## Rolling back

If the new stack has to come down, the old one goes back up from the
`master` branch. Do not run `docker compose` from the `shinyproxy` branch: it
carries test-host commits to the root `nginx/nginx.conf` (sticky cookie,
least_time) that use NGINX Plus-only directives and list only 20 upstreams,
so the nginx image built from that branch does not start.

```bash
cd /docker/idep/shinyproxy
sudo systemctl disable --now shinyproxy 2>/dev/null || ./idep.sh stop   # unit not installed yet? stop by hand
./idep.sh status                          # both DOWN, 0 containers
cd /docker/idep
git checkout master                       # brings back production's committed files (step 2)
docker tag webapp:2026 webapp:latest      # only if ./idep.sh update has ever run: compose runs everything on :latest
sudo docker compose up -d --no-build --scale webapp=50
```

The container memory cap can stay; it simply bounds the 50 compose
containers instead, and `webapp:2026` is just a tag. To lift the cap anyway,
do it **before** `git checkout master` (the script is not on master) and
before `compose up` (it needs a Docker restart, which would bounce every
container):

```bash
sudo ./memory-cap.sh remove && sudo systemctl restart docker
```
