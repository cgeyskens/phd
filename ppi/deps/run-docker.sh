#!/bin/bash

## to build this docker image:
# docker build -t ppi-con  

# to run interactively the image interactively:
# docker container run -it  ppi-con /bin/bash

# adding variables
# the conda env name
CONDA_ENV="ppi"
# the python script
PYTHON_SCRIPT="biogrid-string-intact-apid.py"
# the folder in local host for binding between local host and container
BINDFOLDERS_HOST="/Users/cgeyskens/Documents/code/phd/ppi/deps/container-in-out"
# the folder in container for binding between local host and container
BINDFODLER_CONTAINER="/opt/data"
# protein of interests
PROTEIN_OF_INTERST="CD6_MOUSE"

# run the container
docker container run \
    -v $BINDFOLDERS_HOST:$BINDFODLER_CONTAINER ppi-con \
    /bin/bash -c "source activate $CONDA_ENV \
                 && cd opt/data \
                 && python3 $BINDFODLER_CONTAINER/$PYTHON_SCRIPT $PROTEIN_OF_INTERST"