FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/esga-1.3/bin:/opt/pasa/2.5.1/bin:$PATH
ENV PASAHOME /opt/pasa/2.5.1

RUN apt-get -y update && apt-get -y install wget build-essential zlib1g-dev

RUN mkdir -p /opt/pasa && cd /opt/pasa && wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.5.1/PASApipeline-v2.5.1.FULL.tar.gz && tar -xvf PASApipeline-v2.5.1.FULL.tar.gz && mv PASApipeline-v2.5.1 2.5.1 && rm *.tar.gz
RUN cd /opt/pasa/2.5.1 && make -j2

RUN cp /opt/pasa/2.5.1/pasa_conf/pasa.CONFIG.template /opt/pasa/2.5.1/pasa_conf/conf.txt
RUN cpan -i URI::Encode


