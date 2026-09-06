#!/bin/bash
#
# Cap the memory of every Docker container on this host *together*.
#
#   sudo ./memory-cap.sh install          put all containers under one budget
#   sudo ./memory-cap.sh status           show the budget, usage and OOM kills
#   sudo ./memory-cap.sh remove           lift the cap again
#
# ShinyProxy can only limit memory per container (container-memory-limit,
# 8 GB here) and has no notion of a total, so 200 seats x 8 GB is nominally
# 1.6 TB on a 122 GB host. This puts every container Docker starts into one
# systemd slice, containers.slice, with a MemoryMax. When the containers'
# combined usage reaches it the kernel first drops reclaimable page cache
# inside the slice, then OOM-kills the largest process in it -- one heavy R
# session, which ShinyProxy reports to that user as a crashed app -- while
# the host keeps the remainder for the OS, Docker and the ShinyProxy JVM,
# which live outside the slice.
#
# Run install ONCE per host. Two ordinary config files do the work and both
# survive Docker restarts and reboots, so there is nothing to re-run:
#   /etc/systemd/system/containers.slice     the budget
#   /etc/docker/daemon.json  "cgroup-parent"  makes Docker use it
# install is idempotent and only restarts Docker (and with it every
# container, then ShinyProxy) when daemon.json actually changed. Running it
# again is harmless; it finds nothing to do, or applies an edited budget.
#
# The slice covers *all* Docker containers, not only ShinyProxy's: nginx
# (~70 MB) and anything else on the host (Guacamole here) share the budget.
set -eu

# ---------------------------------------------------------------- the budget
# Combined memory for every container on the host. What is left over goes to
# the OS, the page cache and the ShinyProxy JVM (2 GB heap cap): 16 GiB on a
# 122 GB host is comfortable. Any systemd size works ("90G", "80%").
# To change it later, edit this line and run `sudo ./memory-cap.sh install`
# again: the new value is applied in place, no Docker restart.
MEMORY_MAX=106G
# ----------------------------------------------------------------------------

SLICE=containers.slice
UNIT=/etc/systemd/system/$SLICE
DAEMON_JSON=/etc/docker/daemon.json
CG=/sys/fs/cgroup/$SLICE

die() { echo "error: $*" >&2; exit 1; }
need_root() { [ "$(id -u)" = 0 ] || die "run with sudo"; }

# Which slice a running container is actually in, read from the kernel.
container_slice() {
    local pid; pid=$(docker inspect -f '{{.State.Pid}}' "$1" 2>/dev/null) || return
    [ "${pid:-0}" -gt 0 ] || return
    sed -nE 's|^0::/([^/]+\.slice)/.*|\1|p' "/proc/$pid/cgroup" 2>/dev/null
}

cmd_install() {
    need_root
    local size=$MEMORY_MAX
    [ "$(docker info -f '{{.CgroupDriver}} {{.CgroupVersion}}' 2>/dev/null)" = "systemd 2" ] \
        || die "needs Docker on cgroup v2 with the systemd cgroup driver (docker info)"
    command -v python3 >/dev/null || die "python3 is required to edit $DAEMON_JSON"

    # 1. the slice with its budget
    local unit; unit=$(printf '[Unit]\nDescription=Memory budget shared by every Docker container (see shinyproxy/memory-cap.sh)\n\n[Slice]\nMemoryMax=%s\n' "$size")
    if [ "$(cat $UNIT 2>/dev/null)" != "$unit" ]; then
        printf '%s' "$unit" > $UNIT
        systemctl daemon-reload
        echo "$UNIT: MemoryMax=$size"
    fi
    systemctl start $SLICE
    [ -r $CG/memory.max ] || die "$CG not created; is systemd managing cgroups?"

    # 2. make Docker start every container under it
    local changed
    changed=$(python3 - "$DAEMON_JSON" "$SLICE" <<'PY'
import json, sys
path, slice_ = sys.argv[1], sys.argv[2]
try:
    cfg = json.load(open(path))
except FileNotFoundError:
    cfg = {}
if cfg.get("cgroup-parent") == slice_:
    print("no"); sys.exit()
cfg["cgroup-parent"] = slice_
json.dump(cfg, open(path, "w"), indent=4); open(path, "a").write("\n")
print("yes")
PY
    )
    if [ "$changed" = yes ]; then
        echo "$DAEMON_JSON: cgroup-parent=$SLICE -- restarting Docker (every container restarts)"
        systemctl restart docker
        # ShinyProxy's containers do not survive a Docker restart and it does
        # not adopt what it finds; a restart rebuilds the pool under the slice.
        if systemctl is-active --quiet shinyproxy; then
            systemctl restart shinyproxy
            echo "shinyproxy restarted"
        fi
    fi

    cmd_status
}

cmd_remove() {
    need_root
    python3 - "$DAEMON_JSON" <<'PY'
import json, sys
path = sys.argv[1]
try:
    cfg = json.load(open(path))
except FileNotFoundError:
    sys.exit()
if cfg.pop("cgroup-parent", None) is not None:
    json.dump(cfg, open(path, "w"), indent=4); open(path, "a").write("\n")
PY
    rm -f $UNIT
    systemctl daemon-reload
    echo "cap removed from config; restart Docker (and then shinyproxy) for running containers to leave the slice"
}

cmd_status() {
    local rc=0
    if [ -r $CG/memory.max ]; then
        local max cur
        max=$(cat $CG/memory.max); cur=$(cat $CG/memory.current)
        if [ "$max" = max ]; then
            echo "budget       none (slice exists, no MemoryMax)"; rc=1
        else
            printf 'budget       %.1f GiB, in use %.1f GiB (%d%%)\n' \
                "$(echo "$max/1073741824" | bc -l)" "$(echo "$cur/1073741824" | bc -l)" $((cur * 100 / max))
        fi
        echo "oom kills    $(awk '/^oom_kill /{print $2}' $CG/memory.events) since boot"
    else
        echo "budget       none ($SLICE not present)"; rc=1
    fi
    grep -q "\"cgroup-parent\": \"$SLICE\"" $DAEMON_JSON 2>/dev/null \
        && echo "docker       new containers go into $SLICE" \
        || { echo "docker       cgroup-parent not set in $DAEMON_JSON"; rc=1; }

    local inside=0 outside=0 names=""
    for c in $(docker ps -q 2>/dev/null); do
        if [ "$(container_slice "$c")" = "$SLICE" ]; then inside=$((inside + 1))
        else outside=$((outside + 1)); names="$names $(docker inspect -f '{{.Name}}' "$c" | tr -d /)"; fi
    done
    echo "containers   $inside inside the budget, $outside outside"
    if [ $outside -gt 0 ]; then
        # Containers created before the cap stay where they were until recreated.
        echo "  outside: $(echo $names | tr ' ' '\n' | head -3 | cut -c1-40 | paste -sd ' ')$([ $outside -gt 3 ] && echo " ... ($outside total)") -- started before the cap; restart them"
        rc=1
    fi
    return $rc
}

case "${1:-}" in
    install) cmd_install ;;
    remove)  cmd_remove ;;
    status)  cmd_status ;;
    *) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
