FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for tpp pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# copied from TPP docker file

RUN apt -y update && apt -y upgrade && apt -y install gnuplot && apt -y install unzip && apt -y install libexpat1 libexpat1-dev && apt -y install apache2 && apt -y install build-essential && apt -y install xsltproc && apt clean all 

RUN export PERL_MM_USE_DEFAULT=1 && cpan install CGI && cpan install XML::Parser && cpan install FindBin::libs && cpan install JSON 

WORKDIR /usr/local/
RUN wget "http://www.tppms.org/sw/TPP5.2/TPP5.2_ubuntu.tgz" && tar xvfz TPP5.2_ubuntu.tgz && chmod -R a+rwx tpp && rm -f ${tarfile} 

RUN a2enmod rewrite && a2enmod cgid
WORKDIR /usr/local/tpp/conf/
RUN cat httpd-2.4-tpp.conf | sed 's/#Require all granted/Require all granted/' > httpd-2.4-tpp-1.conf && perl -pi -e 's/_TPP_PORT_/10401/g' httpd-2.4-tpp-1.conf && cp -p /usr/local/tpp/conf/httpd-2.4-tpp-1.conf /etc/apache2/conf-available/httpd-tpp.conf && ln -s /etc/apache2/conf-available/httpd-tpp.conf /etc/apache2/conf-enabled/httpd-tpp.conf && perl -pi -e 's/www-data/biodocker/g' /etc/apache2/envvars

RUN apt -y update &&  apt -y upgrade && apt -y install nano && apt -y install vim && apt clean all 

ENV PERL5LIB="/usr/local/tpp/lib/perl"
ENV PATH="/usr/local/tpp/bin/:/usr/local/tpp/cgi-bin/:${PATH}"
ENV PATH="/opt/conda/envs/nf-core-tpp/bin:${PATH}"

