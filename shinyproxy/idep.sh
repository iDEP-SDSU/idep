#!/bin/bash
#
# Single entry point for the iDEP / ShinyGO ShinyProxy deployment.
#
#   ./idep.sh start | stop | restart | status | logs
#   ./idep.sh check <spec-id>      launch one app and time it
#   ./idep.sh parity               check nothing shiny-server published was dropped
#   ./idep.sh pin <tag>            freeze the current webapp:latest as webapp:<tag>
#   ./idep.sh update [--pull]      rebuild (or pull) webapp:latest, then restart
#
# The stack is a ShinyProxy JAR on the host plus one nginx container on :80.
# App containers are started and reaped by ShinyProxy itself.
set -u

cd "$(dirname "$(readlink -f "$0")")"
ROOT=$(cd .. && pwd)          # the idep repo checkout
PIDFILE=shinyproxy.pid
SP_PORT=8080
NGINX_NAME=sp-nginx

die() { echo "error: $*" >&2; exit 1; }

jar_path() { ls shinyproxy-*.jar 2>/dev/null | sort -V | tail -1; }

# PID of the running ShinyProxy: trust the pidfile, fall back to whoever holds
# the port (covers a JAR started before this script existed).
sp_pid() {
    local pid
    if [ -f "$PIDFILE" ]; then
        pid=$(cat "$PIDFILE" 2>/dev/null)
        if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then echo "$pid"; return 0; fi
    fi
    pid=$(ss -tlnpH "sport = :$SP_PORT" 2>/dev/null | grep -oP 'pid=\K[0-9]+' | head -1)
    [ -n "$pid" ] && { echo "$pid"; return 0; }
    return 1
}

# Every container ShinyProxy started, whether or not ShinyProxy is still alive.
app_containers() { docker ps -aq --filter "label=openanalytics.eu/sp-proxy-id"; }

cmd_start() {
    local jar; jar=$(jar_path)
    [ -n "$jar" ] || die "no shinyproxy-*.jar here; download it from https://www.shinyproxy.io/downloads/"

    # Fail loudly now rather than on every app launch: a missing pinned image
    # means every legacy app 500s while /idep/ and /go/ look fine.
    for img in $(grep -oP '^\s*container-image:\s*\K\S+' application.yml | sort -u); do
        docker image inspect "$img" >/dev/null 2>&1 || die "image $img not found (./idep.sh pin <tag> creates the pinned tag)"
    done

    if sp_pid >/dev/null; then
        echo "ShinyProxy already running (pid $(sp_pid))"
    else
        setsid nohup java -jar "$jar" > startup.log 2>&1 < /dev/null &
        echo $! > "$PIDFILE"
        disown
        printf 'starting ShinyProxy (%s)' "$jar"
        for _ in $(seq 1 60); do
            curl -s -o /dev/null --max-time 2 "http://127.0.0.1:$SP_PORT/" && break
            kill -0 "$(cat $PIDFILE)" 2>/dev/null || { echo; die "ShinyProxy exited; see startup.log"; }
            printf '.'; sleep 2
        done
        echo
        curl -s -o /dev/null --max-time 2 "http://127.0.0.1:$SP_PORT/" \
            || die "ShinyProxy did not answer on :$SP_PORT; see startup.log"
        echo "ShinyProxy up on :$SP_PORT (pid $(cat $PIDFILE), log: startup.log)"
    fi

    docker rm -f "$NGINX_NAME" >/dev/null 2>&1
    docker run -d --name "$NGINX_NAME" --network host --restart unless-stopped \
        -v "$PWD/nginx/nginx.conf:/etc/nginx/nginx.conf:ro" \
        -v "$ROOT/shinyapps/dist:/srv/dist:ro" \
        -v "$ROOT/data/shared:/srv/shared:ro" \
        nginx:latest >/dev/null || die "nginx failed to start"
    echo "nginx up on :80"
    echo "  http://$(hostname -I | awk '{print $1}')/"
}

cmd_stop() {
    # nginx first, so no request can land on a half-stopped ShinyProxy.
    docker rm -f "$NGINX_NAME" >/dev/null 2>&1 && echo "nginx stopped"

    local pid
    if pid=$(sp_pid); then
        kill "$pid" 2>/dev/null
        for _ in $(seq 1 15); do kill -0 "$pid" 2>/dev/null || break; sleep 1; done
        kill -0 "$pid" 2>/dev/null && kill -9 "$pid" 2>/dev/null
        echo "ShinyProxy stopped (pid $pid)"
    else
        echo "ShinyProxy not running"
    fi
    rm -f "$PIDFILE"

    # ShinyProxy reaps its containers on a clean shutdown, but not if it was
    # killed; sweep whatever is left so the next start begins from zero.
    # Errors are ignored: ShinyProxy's own reaping runs concurrently, so some
    # of the ids listed a moment ago are already gone by the time we get here.
    local before after
    before=$(app_containers | wc -l)
    if [ "$before" -gt 0 ]; then
        app_containers | xargs -r docker rm -f >/dev/null 2>&1
        after=$(app_containers | wc -l)
        echo "app containers: $before -> $after"
    fi
}

cmd_status() {
    local pid rc=0
    if pid=$(sp_pid); then
        echo "ShinyProxy   running (pid $pid, :$SP_PORT)"
    else
        echo "ShinyProxy   DOWN"; rc=1
    fi

    if [ -n "$(docker ps -q -f "name=^${NGINX_NAME}$")" ]; then
        echo "nginx        running (:80)"
    else
        echo "nginx        DOWN"; rc=1
    fi

    local specs
    specs=$(curl -s --max-time 5 "http://127.0.0.1:$SP_PORT/api/proxyspec" 2>/dev/null \
            | grep -o '"id":' | wc -l)
    echo "specs        ${specs:-0} configured"

    local n mem
    n=$(docker ps -q --filter "label=openanalytics.eu/sp-proxy-id" | wc -l)
    if [ "$n" -gt 0 ]; then
        mem=$(docker stats --no-stream --format '{{.MemUsage}}' \
              $(docker ps -q --filter "label=openanalytics.eu/sp-proxy-id") 2>/dev/null \
              | awk '{s+=$1} END {printf "%.1f", s/1024}')
        echo "containers   $n up (~${mem:-?} GB)"
    else
        echo "containers   0 up"
    fi

    # Static site only. Deliberately does not touch the app paths: a GET to a
    # legacy path would boot a container just to answer the health check.
    local code
    code=$(curl -so /dev/null -w '%{http_code}' --max-time 5 http://127.0.0.1/ 2>/dev/null)
    echo "static site  HTTP ${code:-000}"
    [ "$code" = "200" ] || rc=1

    return $rc
}

cmd_check() {
    local id=${1:-}
    # Declared before the "all" branch: the RETURN trap set by the recursive
    # calls below is global, and fires on this frame's return too.
    local cj=""
    [ -n "$id" ] || die "usage: ./idep.sh check <spec-id>       (or: check all)"

    if [ "$id" = "all" ]; then
        local rc=0
        for spec in $(curl -s --max-time 5 "http://127.0.0.1:$SP_PORT/api/proxyspec" \
                      | grep -oP '"id":"\K[^"]+'); do
            cmd_check "$spec" || rc=1
        done
        return $rc
    fi

    # A plain GET only fetches ShinyProxy's loading page, which returns 200
    # before the container exists. Drive the API the way the browser does:
    # POST to start the session, then poll the app's own URL until R answers.
    local resp proxy_id app_path start code
    cj=$(mktemp); trap 'rm -f "$cj"' RETURN
    curl -s -c "$cj" -b "$cj" -L -o /dev/null "http://127.0.0.1:$SP_PORT/"

    start=$(date +%s.%N)
    resp=$(curl -s -c "$cj" -b "$cj" -X POST -H 'Content-Type: application/json' \
           -d '{}' --max-time 200 "http://127.0.0.1:$SP_PORT/api/proxy/$id")
    proxy_id=$(echo "$resp" | grep -oP '"id":"\K[^"]+' | head -1)
    app_path=$(echo "$resp" | grep -oP '"SHINYPROXY_PUBLIC_PATH":"\K[^"]+')

    if [ -z "$app_path" ]; then
        printf '%-10s FAILED to start: %s\n' "$id" "$(echo "$resp" | head -c 200)"
        return 1
    fi

    code=000
    for _ in $(seq 1 1200); do
        code=$(curl -s -c "$cj" -b "$cj" -o /dev/null -w '%{http_code}' --max-time 5 \
               "http://127.0.0.1:$SP_PORT$app_path")
        [ "$code" = "200" ] && break
        sleep 0.1
    done

    printf '%-10s HTTP %s  %.1fs\n' "$id" "$code" "$(echo "$(date +%s.%N) - $start" | bc)"

    # Release the seat so a check does not leave a container idling for 15 min.
    [ -n "$proxy_id" ] && curl -s -c "$cj" -b "$cj" -X DELETE \
        "http://127.0.0.1:$SP_PORT/api/proxy/$proxy_id" >/dev/null

    [ "$code" = "200" ]
}

cmd_pin() {
    local tag=${1:-}
    [ -n "$tag" ] || die "usage: ./idep.sh pin <tag>    e.g. ./idep.sh pin 2026"
    docker image inspect webapp:latest >/dev/null 2>&1 || die "webapp:latest not found"
    docker tag webapp:latest "webapp:$tag" || die "tag failed"
    echo "webapp:latest frozen as webapp:$tag"
    echo "point legacy specs in application.yml at webapp:$tag to hold them there."
}

cmd_update() {
    # Freeze the outgoing image before it is replaced. Legacy specs reference
    # a dated tag, so this snapshot is what keeps them reproducible; without it
    # the environment they were tested against is gone for good.
    local stamp; stamp=$(date +%Y%m%d)
    if docker image inspect webapp:latest >/dev/null 2>&1; then
        docker tag webapp:latest "webapp:pre-$stamp"
        echo "outgoing image saved as webapp:pre-$stamp"
    fi

    if [ "${1:-}" = "--pull" ]; then
        docker pull gexijin/idep:latest || die "pull failed"
        docker tag gexijin/idep:latest webapp:latest
    else
        ( cd "$ROOT" && docker build -t webapp:latest . ) || die "build failed"
    fi

    echo "restarting so running containers pick up the new image..."
    cmd_stop
    cmd_start
}

cmd_parity() {
    # Proves every app shiny-server published is still published here, and that
    # each one still points at the same directory. Run it on production before
    # decommissioning shiny-server, and after editing either config.
    #
    # Delete this subcommand once config/shiny-server.conf is gone -- at that
    # point application.yml is the only source of truth and there is nothing
    # left to compare against.
    local ss="$ROOT/config/shiny-server.conf"
    [ -f "$ss" ] || die "$ss not found (already retired? then so is this check)"

    # spec id -> the shinyapps/ directory it runs, read out of application.yml
    local specs
    specs=$(awk '
        function flush() { if (id != "") print id "\t" (dir != "" ? dir : pkg); id=""; dir=""; pkg="" }
        /^    - /       { flush() }
        /^    - id: /   { id=$3 }
        /^      id: /   { id=$2 }
        # container-cmd and the go86 volume both name the app directory; the
        # legacy anchor mounts bare /srv/shiny-server, which this will not match.
        match($0, /\/srv\/shiny-server\/[A-Za-z0-9_.-]+/) {
            dir = substr($0, RSTART + 18, RLENGTH - 18)
        }
        # /idep/ runs the installed package, not a directory under shinyapps.
        /idep250::run_app/ { pkg = "idep250" }
        END { flush() }
    ' application.yml)

    # url -> spec id, read out of nginx.conf
    local routes
    routes=$(sed -nE 's|^ *location (/[A-Za-z0-9_-]+)/ *\{ *proxy_pass http://shinyproxy/app/([A-Za-z0-9_-]+)/; *\}.*|\1 \2|p' nginx/nginx.conf)

    local rc=0 n=0
    while read -r url want; do
        [ -n "$url" ] || continue
        n=$((n + 1))
        local id got
        id=$(echo "$routes" | awk -v u="$url" '$1 == u {print $2; exit}')
        if [ -z "$id" ]; then
            printf 'MISSING  %-10s no nginx route\n' "$url"; rc=1; continue
        fi
        got=$(echo "$specs" | awk -F'\t' -v i="$id" '$1 == i {print $2; exit}')
        if [ -z "$got" ]; then
            printf 'MISSING  %-10s nginx routes to spec "%s", which is not in application.yml\n' "$url" "$id"; rc=1
        elif [ "$got" != "$want" ]; then
            printf 'MISMATCH %-10s shiny-server ran %s, spec "%s" runs %s\n' "$url" "$want" "$id" "$got"; rc=1
        fi
    # shiny-server.conf has CRLF line endings; strip them or every directory
    # comes back with a trailing \r and nothing ever compares equal.
    done <<< "$(awk '
        { sub(/\r$/, "") }
        /^[[:space:]]*location /  { loc = $2 }
        /app_dir/ {
            gsub(/;/, "", $2); sub(/.*\/srv\/shiny-server\//, "", $2)
            print loc, $2
        }' "$ss")"

    if [ $rc -eq 0 ]; then
        echo "parity ok: all $n shiny-server locations are published, on the same directories"
    fi
    return $rc
}

cmd_logs() { tail -n "${2:-50}" -f startup.log; }

case "${1:-}" in
    start)   cmd_start ;;
    stop)    cmd_stop ;;
    restart) cmd_stop; echo; cmd_start ;;
    status)  cmd_status ;;
    check)   shift; cmd_check "$@" ;;
    pin)     shift; cmd_pin "$@" ;;
    update)  shift; cmd_update "$@" ;;
    parity)  cmd_parity ;;
    logs)    cmd_logs "$@" ;;
    *)       sed -n '3,12p' "$0" | sed 's/^# \?//'; exit 1 ;;
esac
