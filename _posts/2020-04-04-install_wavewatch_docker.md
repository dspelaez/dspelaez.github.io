---
layout: post
title: Install WAVEWATCH III wave model in a Docker container
date: 2020-04-04
author: D.S. Peláez-Zapata
categories: Posts
tags: wave-modelling, docker, containers	
thumbnail: "img/thumbnails/wavewatch.jpg"
---

![docker]({{ site.url }}/img/thumbnails/docker-art.png){: .center-image }

[Docker](https://www.docker.com/) is an open source project that make extremely
easy deployment of software and applications thought the use of containers.  A
container, as indicated in the docker project web site, is a "standard unit of
software that packages up code and all its dependencies so the application runs
quickly and reliably from one computing environment to another."

It is not a secret that one of the most difficult steps when running a numerical
or computational model is the compilation/installation, not only the model
itself, but the huge amount of dependencies that sometimes are required. An
undergrad or grad student can spend more time compiling the model than actually
running it in real or test cases. In that sense, it is a really good idea to use
a docker container to install and deploy the model along with all of its
dependencies, that way one can distribute it, as a package, without worrying
about to compile and install confusing libraries. Everything is in the docker
image and it works out-of-the-box.

OK, let's see how to create an image containing the [WAVEWATCH III
5.16](https://polar.ncep.noaa.gov/waves/wavewatch/) using docker and how to use
the model afterwards.


### Creating the Dockerfile

The first part is to create a file called `Dockerfile` that contains all the
steps to create the image. The header of out docker file is:

```
FROM ubuntu:18.04
MAINTAINER Daniel Santiago <dspelaez@gmail.com>
```

which states that we want to use Ubuntu 18.04 as a base OS. Then we need to
prepare the update the system, install the dependencies and prepare the folder
to copy the model source code. In that case WAVEWATCH III depends on `gfortran`
and `gcc` as Fortran and C/C++ compilers:

```
# set home variable
ENV HOME /root
ENV TERM xterm

RUN apt-get update
RUN apt-get -yq install build-essential m4 \
    gcc gfortran g++ \
    curl libcurl4-gnutls-dev wget \
    git tar

# set environmental variables
ENV NCDIR /usr/local/netcdf
ENV WWATCH3_NETCDF NC4
ENV NETCDF_CONFIG /usr/local/netcdf/bin/nf-config

# create folders to put libs and ww3 source in
RUN mkdir -p /root/libs /root/ww3
```

This is the most tricky part. WAVEWATCH III uses netCDF4 to handle the input and
output files (can be plain text files as well, but for large datasets it becomes
unmanageable). To install those libraries in the image we can use the following
[script](https://raw.githubusercontent.com/dspelaez/install-netcdf-fortran/master/install_netcdf.sh):

```
# get model source code
ADD https://raw.githubusercontent.com/dspelaez/install-netcdf-fortran/master/install_netcdf.sh /root/libs
WORKDIR /root/libs
RUN chmod 755 ./install_netcdf.sh && \
    CC=gcc FC=gfortran MAINDIR=/usr/local/netcdf ./install_netcdf.sh && \
    export LD_LIBRARY_PATH=/usr/local/netcdf/lib:$LD_LIBRARY_PATH && \
    echo export LD_LIBRARY_PATH=/usr/local/netcdf/lib:'$LD_LIBRARY_PATH' >> /root/.bashrc
```

If this work well, we are almost done. The next step is to install the model
itself. Since the 5.16 version the NOAA started to use Github to host the source
code, it makes everything easier. So the only thing we need to do is clone the
repo and run the installation script using some tricks to avoid the prompting of
the questions. Note this will work unless the model developers change the order
of the questions in future versions.

```
# get tar file from repository
WORKDIR /root/ww3
RUN wget https://github.com/NOAA-EMC/WW3/releases/download/6.07/wwatch3.v6.07.tar.gz

# untar and run instalation script
RUN tar -xvf wwatch3.v6.07.tar.gz && \
    chmod 755 ./install_ww3_tar && \
    printf "y\nL\ny\ny\nprinter\ngfortran\ngcc\n/root/ww3/model/tmp\nyes\nyes\ny\ny" > answers && \
    ./install_ww3_tar < answers
```

The next step is configure and compile the model. This part is well explained in
the user manual: we need to have two files `comp` and `link`. Fortunately there
are some preset versions distributed with the source so the only thing we need
to do is copy them to the `bin` folder and compile the entire code:

```
# copy switch file
COPY switch /root/ww3/model/bin
RUN cp /root/ww3/model/bin/link.Gnu /root/ww3/model/bin/link
RUN cp /root/ww3/model/bin/comp.Gnu /root/ww3/model/bin/comp

# compile the entire code
RUN /root/ww3/model/bin/w3_clean && \
    /root/ww3/model/bin/w3_new   && \
    /root/ww3/model/bin/w3_make
```

Finally, let's delete some unnecessary files, just to have a lighter image, and
then set the working directory:

```
# delete unnecessary files
RUN rm -rf /root/ww3/README_manual_reference /root/ww3/manual.*.pdf \
           /root/ww3/arc /root/ww3/regtests /root/ww3/smc_docs \
           /root/libs

# set working directory
WORKDIR /root/ww3/model/work
```

The `Dockerfile` is ready and now we can create the image.

### Building the image

A Docker image is just a file containing a set of multiple layers, according to
the instructions specified in the `Dockerfile`, that is used to execute code in
a Docker container. More information of the difference between container,
images, client, kernel etc., can be found
[here](https://docs.docker.com/engine/docker-overview/).

To build the image we just type (within the folder containing the `Dockerfile`):

```
docker build -t ww3-docker .
```

where `ww3-docker` is as we will call our new image, it can be anything else.
Now the image containing the WAVEWATCH III 5.16 complied model and all its
dependences has been created. It is time to run the mode.


### Running the model

Once the image have been built, now we can run the model interactively with:

```
docker run -it --rm --name=ww3 ww3-docker /bin/bash
```

Of course we will want to run the model using our own data, that is to say, our
own bottom, wind, ice, currents and other forcings, and our own config files. To
do so the logic is simple, the only thing we have to do is to mirror the
container work directory with our own work directory, in other words, we have to
mount a volume that shares the host data with the container instance. To do that
just type:

```
docker run -it --rm --name=ww3 -v $(pwd):/root/ww3/work ww3-docker /bin/bash
```

in this case, `$(pwd)` represents the current directory, i.e., where our input
data is stored.
