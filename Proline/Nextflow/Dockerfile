FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-proline-labelfree"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-proline-labelfree/bin:$PATH
# files will be changed by searchgui and thus need to be readable by any user
RUN apt-get update; apt-get install -y unzip  ; apt-get clean 
RUN chmod -R a+rwx /opt/conda/envs/nf-core-proline-labelfree/share/searchgui-*/resources
RUN chmod -R a+rwx /opt/conda/envs/nf-core-proline-labelfree/share/peptide-shaker-*/resources
RUN update-ca-certificates
RUN mkdir /usr/share/man/man1
RUN mkdir  mzdbconvert; cd mzdbconvert; wget https://bitbucket.org/veitveit/prolineworkflow/downloads/mzML2mzDB.zip; \
unzip mzML2mzDB.zip; cd ..
RUN chmod a+x /mzdbconvert/mzML2mzDB; ln -s /mzdbconvert/mzML2mzDB /usr/bin/mzML2mzDB
RUN wget https://github.com/profiproteomics/proline-cli/releases/download/0.2.0-SNAPSHOT-2019-10-04/proline-cli-0.2.0-SNAPSHOT-bin.zip
RUN apt-get update; apt-get install -y unzip openjdk-11-jdk; apt-get clean 
# Very dirty solution to get java 8 in as it is needed for Proline
RUN apt-get update; apt-get install -y apt-transport-https ca-certificates wget dirmngr gnupg software-properties-common
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add -
RUN add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/
RUN apt update; apt install -y adoptopenjdk-8-hotspot
###
RUN unzip proline-cli-0.2.0-SNAPSHOT-bin.zip; mv proline-cli-0.2.0-SNAPSHOT proline
RUN chmod -R a+rwx proline/
RUN mkdir proline_results
RUN chmod -R a+rwx proline_results/

