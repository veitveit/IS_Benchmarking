FROM nfcore/base
LABEL authors="Niclas Kildegaard Nielsen" \
      description="Docker image containing all requirements for maxquant pipeline"

COPY environment.yml /
# Create the environment:
RUN conda env create -f /environment.yml && conda clean -a
# Make RUN commands use the new environment:
# Copies the SDRF parser
COPY sdrf/ ./sdrf/

# Activates the env for the program to run in
RUN echo "source activate nf-core-maxquant" > ~/.bashrc
ENV PATH /opt/conda/envs/nf-core-maxquant/bin:$PATH

# Changes the Work dir to the SDRF file area
WORKDIR /sdrf
RUN python setup.py install
WORKDIR /
# ENTRYPOINT ["/bin/bash"]

