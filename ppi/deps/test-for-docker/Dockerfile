# this dockerfile lets me setup a conda env inside a docker where the scripts and the data are local but the analysis happens in the contianer

# # run this container interactively 
# # docker container run -it myimage4 /bin/bash 

# Mount via data directory
# docker container run -v /Users/cgeyskens/Documents/code/phd/ppi/conda_apptainer_docker/test/data:/opt/data myimage4 /bin/bash -c "source activate test_env && python3 opt/data/test.py"

FROM continuumio/miniconda3:24.1.2-0

# Copy the Python file & environment file
COPY test_env.yml /opt

# make mounting directory
RUN mkdir /opt/data

# Initialize Conda
RUN conda init bash

# Create the environment and install required packages
RUN conda env create -f /opt/test_env.yml && \
    echo "conda activate test_env" >> ~/.bashrc && \
    /bin/bash -c "source ~/.bashrc && conda install pandas"

# Debugging statement
RUN echo "Contents of ~/.bashrc:" && cat ~/.bashrc

# Start bash shell with activated Conda environment
CMD ["conda", "run", "-n", "test_env", "/bin/bash"]
