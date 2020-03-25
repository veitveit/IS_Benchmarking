FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-moff-ms1"


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a;

RUN conda install python=3.6; \
	conda install -c conda-forge -c bioconda moff;


ENV PATH /opt/conda/envs/nf-core-moff-ms1/bin:$PATH


SHELL ["conda", "run", "-n", "nf-core-moff-ms1", "/bin/bash", "-c"]

