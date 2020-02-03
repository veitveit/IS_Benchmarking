FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-searchgui-labelfree"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-searchgui-labelfree/bin:$PATH
# files will be changed by searchgui and thus need to be readable by any user
RUN chmod -R a+rwx /opt/conda/envs/nf-core-searchgui-labelfree/share/searchgui-3.3.10-0/resources
RUN chmod -R a+rwx /opt/conda/envs/nf-core-searchgui-labelfree/share/peptide-shaker-1.16.36-0/resources
