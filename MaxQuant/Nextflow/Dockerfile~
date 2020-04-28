FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for tpp pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
COPY MergeOutput.R  .
ENV PATH /opt/conda/envs/nf-core-tpp/bin:$PATH
