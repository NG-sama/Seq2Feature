#Run this docker file to automatically create a container and run the Seq2Feature python 
# Use the Miniconda3 base image
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /seqannotate_v1

# Copy the environment.yml file into the container
COPY environment.yml .

# Create the Conda environment using the environment.yml file
RUN conda env create -f environment.yml \
    && conda clean --all -f -y

# Activate the Conda environment
SHELL ["conda", "run", "-n", "seqannotate", "/bin/bash", "-c"]

# Copy the rest of your application code into the container
COPY . .
# Ensure the Conda environment is activated for subsequent commands
RUN source activate seqannotate && pip install .

# Set the default command to run when the container starts
CMD ["conda", "run", "-n", "seqannotate", "python", "-m", "seqannotate.main"]