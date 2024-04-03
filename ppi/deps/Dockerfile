FROM continuumio/miniconda3:24.1.2-0

# copy the conda environment file
COPY ppi-env.yml /opt

# make mounting directory
RUN mkdir /opt/data

# initialize Conda
RUN conda init bash

# create the environment and install required packages
RUN conda env create -f /opt/ppi-env.yml && \
    echo "conda activate ppi" >> ~/.bashrc && \
    /bin/bash -c "source ~/.bashrc"

# debugging statement
RUN echo "Contents of ~/.bashrc:" && cat ~/.bashrc

# start bash shell with activated Conda environment
CMD ["conda", "run", "-n", "ppi", "/bin/bash"]

# to build this container: docker build -t ppi-con .
# to run this container, see run-docker.sh