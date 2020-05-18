FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-moff-ms1"


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a;


#SHELL ["conda", "run", "-n", "moFF", "/bin/bash", "-c"]


ENV PATH /opt/conda/envs/moFF/bin:$PATH


#ENTRYPOINT ["conda", "run", "-n", "moFF", "/bin/bash", "-c"]
