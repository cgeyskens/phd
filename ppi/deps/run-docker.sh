#!/bin/bash

# adding variables
CONDA_ENV="ppi"
PYTHON_SCRIPT="biogrid-string-intact-apid.py"
BINDFOLDERS_HOST="/Users/cgeyskens/Documents/code/phd/ppi/deps/container-in-out"
BINDFODLER_CONTAINER="/opt/data"
PROTEIN_OF_INTERST="CD6_MOUSE"

# run the container
docker container run -it \
    -v $BINDFOLDERS_HOST:$BINDFODLER_CONTAINER ppi-con \
    /bin/bash -c "source activate $CONDA_ENV \
                 && cd opt/data \
                 && python3 $BINDFODLER_CONTAINER/$PYTHON_SCRIPT $PROTEIN_OF_INTERST"