FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for tpp pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# copied from TPP docker file
WORKDIR /usr/local/

RUN apt-get update && apt-get install -y libfindbin-libs-perl unzip && apt-get clean all
RUN apt-get update && apt-get install -y libxml-sax-expat-perl && apt-get clean all

RUN wget https://github.com/BioContainers/software-archive/releases/download/TPP/tpp-5.0.zip && \
  unzip tpp-5.0.zip
RUN chmod a+x /usr/local/tpp/bin/*

ENV PATH /usr/local/tpp/bin/:$PATH

COPY MergeOutput.R  .
ENV PATH /opt/conda/envs/nf-core-tpp/bin:$PATH
