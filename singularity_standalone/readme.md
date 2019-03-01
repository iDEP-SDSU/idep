# Build script and wrapper to run iDEP in standalone mode using singularity

[Singularity](https://www.sylabs.io/) containers are an alternative to docker
containers often found at HPC sites. Tools in this directory allow creation of
an iDEP installation that uses singularity.

### Requirements

- singularity > 3
- the container must be built on a system where you have sudo privileges
  with access to the internet to download data files

## Set up the application directory and build the container

```ShellSession
$ ./build_idep.sh /path/to/install/prefix
```

Notes:

- This would best be installed to an isolated prefix, not `/usr/local`.
- Currently this uses the head of the iDEP repo

### Running the shiny server

```ShellSession
$ /path/to/install/prefix/bin/idep-server
```

Notes:

- this will pick a random port between 49000 and 50000
- if running this on a compute node, you may need to set up a tunnel to access the server

