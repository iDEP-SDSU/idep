# Initial setup on a new host

One-time steps to bring the ShinyProxy stack up on a server that currently
runs the master-branch stack (nginx + shiny-server via `docker-compose.yml`).
Day-to-day operation is in [README.md](README.md).

Everything below assumes the repo is checked out at `/home/exouser/idep` and
that the login user is `exouser`. See step 4 if that is not the case.

## 1. Host prerequisites

The script starts ShinyProxy as a plain JAR on the host, so the host needs:

| what | why | check |
| --- | --- | --- |
| Java 17 or newer | ShinyProxy 3.2.x is a Spring Boot JAR | `java -version` |
| `unzip` | `idep.sh start` extracts `app.html` from the JAR | `command -v unzip` |
| Docker, with the login user in the `docker` group | ShinyProxy talks to `/var/run/docker.sock` | `docker ps` without sudo |
| `nginx:1.30` image pulled | pinned in `idep.sh`; pre-pull so a reboot does not depend on Docker Hub | `docker pull nginx:1.30` |

```bash
sudo apt-get install -y openjdk-21-jre-headless unzip
docker pull nginx:1.30
```

`idep.sh` does not check for Java. Without it, `start` reports ShinyProxy as
failed and `startup.log` contains only `java: command not found`.

## 2. Check out the branch and its submodules

```bash
cd /home/exouser/idep
git fetch origin
git checkout shinyproxy
git submodule update --init      # idep11, go80, idepgolem1, go
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
cd /home/exouser/idep/shinyproxy
wget https://www.shinyproxy.io/downloads/shinyproxy-3.2.4.jar
ls shinyproxy-*.jar
```

`idep.sh` picks the highest-versioned `shinyproxy-*.jar` it finds here.

## 4. Adjust paths if the checkout or user differs

Two files hardcode `/home/exouser/idep` and the user `exouser`:

- `application.yml`: the `container-volumes` bind mounts (8 lines).
- `shinyproxy.service`: `User=`, `PIDFile=`, `WorkingDirectory=`,
  `ExecStart=`, `ExecStop=`.

If production uses a different path or login, replace them before starting:

```bash
sed -i 's#/home/exouser/idep#/actual/path/to/idep#g' application.yml shinyproxy.service
sed -i 's#^User=exouser#User=actualuser#' shinyproxy.service
```

`idep.sh` and `nginx/nginx.conf` derive paths from where the script lives and
need no change.

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
mkdir -p /home/exouser/idep/usage
```

If it is missing Docker creates it as root; that also works, since the app
containers run as root, but the directory then cannot be removed without sudo.

## 7. TLS

`start` uses `/etc/pki/nginx/server.pem` and
`/etc/pki/nginx/private/server.key` when they exist, the same files the old
nginx used, so production needs nothing here. On a host without them set
`CERT_PEM` and `CERT_KEY` (see README, "TLS"). `start` refuses to run without
both files.

## 8. Stop the old stack

Both stacks bind :80 and :443.

```bash
cd /home/exouser/idep
sudo docker compose down
```

Also disable anything that would bring it back: a cron entry or systemd unit
that runs `restart_server.sh` or `docker compose up`, and remove
`restart: always` containers left from earlier compose versions
(`docker ps -a --filter name=idep`).

## 9. Install the container memory cap

Do this **before** starting the stack. It edits `/etc/docker/daemon.json` and
restarts Docker, which would restart every running container.

```bash
cd /home/exouser/idep/shinyproxy
sudo ./memory-cap.sh install      # budget is MEMORY_MAX at the top of the script (140G)
sudo ./memory-cap.sh status
```

Keep at least 16 GiB of host RAM outside the budget. On a host with less RAM
than `MEMORY_MAX` the slice never binds; edit the value first.

## 10. First start and verification

```bash
./idep.sh start
./idep.sh status
./idep.sh parity                  # every shiny-server location has a route and a spec
./idep.sh check all               # ~2 min; must be 35 of 35 on production
```

`parity` compares against `config/shiny-server.conf`; anything it reports
means a URL from the old server would now 404. Any `check` failure on
production is a real problem, since all data directories and app versions
exist there (the test host fails 13 only because they do not).

Then from a browser, over https:

- `/` landing page, `/data/` listing (also over plain http).
- `/idep/` and `/go/`, and a couple of old versions such as `/idep73/`,
  `/go41/`.
- Upload the example data in `/idep/` and run through to a plot, so the
  database and `countsData` mounts are exercised, not only the app start.
- Follow the "old versions" link inside iDEP; it should open in a new tab or
  take over the tab, never nest inside the frame.

## 11. Install the systemd unit

```bash
sudo cp shinyproxy.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable shinyproxy
./idep.sh stop && sudo systemctl start shinyproxy
systemctl status shinyproxy
```

From now on use `systemctl start|stop|restart shinyproxy` rather than
`./idep.sh` directly, so systemd's view stays right. Reboot once and confirm
the site answers without manual intervention.

## 12. Afterwards

- Watch `docker stats` and `./idep.sh status` during the first busy day to
  see how far real traffic sits from the memory budget and the seat cap.
- Once the old stack is gone for good, `config/shiny-server.conf` and the
  `parity` subcommand can be retired; `application.yml` becomes the only
  source of truth for which apps exist.
