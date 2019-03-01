#! /bin/bash

# Build script to set up a single instance of the idep shiny server running in
# a singularity container. No kind of load balancing is done - only the
# built-in shiny web server is running **as the user**

################################################################################
#                                  Functions                                   #
################################################################################
function fail() {
    printf "### ERROR ### %s\n" "$@" >&2
    exit 1
}

function info() {
    printf "### INFO  ### %s\n" "$@" >&2
}

function download() {
    local name="${1}"
    local url="${2}"
    if [[ -e "${name}" ]]; then
        info "${name} already exists - skipping download"
    else
        info "downloading ${name}"
        tmp="${name}.tgz"
        (
          wget -O "${tmp}" "${url}" \
          && tar -xzf "${tmp}" \
          && rm -f "${tmp}" \
        ) || fail "download failed"
    fi
}

################################################################################
#                                   Settings                                   #
################################################################################
wd="$(pwd)"
def_file="${wd}/idep.def"
img_file="${wd}/idep.sif"
wrapper="${wd}/idep-server"

################################################################################
#                                     Main                                     #
################################################################################
prefix="${1:-none}"

if [[ "${prefix}" == "none" ]]; then
    printf "USAGE: build_idep <install prefix>\n"
    exit 1
fi


###
### set up the folder structure for the application
###
cd "${prefix}" || fail "prefix for install does not exist"
mkdir -p libexec share/idep/data bin || fail "could not create directories in ${prefix}"

cp "${wrapper}" bin

###
### build the container - this takes a long time
###

if command -v singularity &> /dev/null; then
    info "Found singularity version $(singularity --version)"
else
    fail "singularity not found"
fi


if [[ ! -f libexec/idep.sif ]]; then
    if [[ -f "${img_file}" ]]; then
        info "using existing singularity image in build script directory"
        cp "${img_file}" libexec/idep.sif
    else
        [[ -f "${def_file}" ]] || fail "container definition file does not exist"
        info "Building container"
        sudo singularity build libexec/idep.sif share/idep/idep.def || fail "singularity build failed"
    fi
    cp "${def_file}" share/idep
else
    info "Container already exists"
fi

###
### fetch the data
###
cd share/idep/data || fail "could not cd to ${prefix}/share/idep/data"
download pathwayDB https://sdsu.box.com/shared/static/c24f792ojoikpzu0lkpng8uuf9ychwm7.gz
download motif https://sdsu.box.com/shared/static/9v1ao6mwhduvrcx793j3answph9gqnkt.gz
download geneInfo https://sdsu.box.com/shared/static/mns0k1uvwtfnsohoc89b984ih36nmnz9.gz
download data_go https://sdsu.box.com/shared/static/qwpdh36vcisgy1hcmadck8i8ezhvr2fh.gz
download convertIDs.db https://sdsu.box.com/shared/static/sorewt7w6iypmhg2k2xhyi8myeit156o.gz
