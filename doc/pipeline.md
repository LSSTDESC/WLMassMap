# WL mass-mapping pipeline

This document details how to setup and run the mass-mapping pipeline

## Basic local setup

This implementation of the pipeline relies on [descpipe](https://github.com/joezuntz/descpipe)
which is fairly easy to install:
```
pip3 install --user pyyaml py-dag descpipe
```
and make sure that the descpipe executable is in your path.

Each step of the pipeline is implemented as a Docker container, which first need
to be built. In the root folder of the project, run:
```
descpipe build pipeline/pipeline.yaml
```

The next step is to generate the bash script that runs the pipeline:
```
descpipe local pipeline/pipeline.yaml local.sh
```

The pipeline can then be executed by running the script:
```
bash local.sh
```

## Deploying at NERSC

The first step is to push local docker images to NERSC's private registry, here
is how to do it. You first need to open an SSH tunnel to NERSC on a random port:
```
ssh -L localhost:6000:registry.services.nersc.gov:443 yourusername@edison.nersc.gov
```
while keeping that connection open, in another terminal push the pipeline
containers to NERSC:
```
descpipe push pipeline/pipeline.yaml
```

Once  the transfer is completed, you need to pull the docker images at nersc.
On a terminal running at nersc:
```
cd pipeline/images
make pull
```
This takes a little while but makes the docker images available to shifter.

The last step is to generate an execution script at NERSC:
```
descpipe nersc pipeline/pipeline.yaml nersc.sh
```

Running the pipeline:
```
salloc -N 1 -p debug -t 00:10:00
bash nersc.sh
```
