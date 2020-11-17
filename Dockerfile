FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/esga-1.1/bin:/opt/bin:/opt/spaln/2.4.1/bin:$PATH
ENV PASAHOME /opt/conda/envs/esga-1.1/opt/pasa-2.4.1
ENV EVM_HOME /opt/conda/envs/esga-1.1/opt/evidencemodeler-1.1.1
RUN cp /opt/conda/envs/esga-1.1/opt/pasa-2.4.1/pasa_conf/pasa.CONFIG.template /opt/conda/envs/esga-1.1/opt/pasa-2.4.1/pasa_conf/conf.txt
RUN mkdir -p /opt/bin && cd /opt/bin && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords && chmod +x faSomeRecords
RUN apt-get -y update && apt-get -y install make gcc
RUN mkdir -p /opt/spaln && cd /opt/spaln && wget https://github.com/ogotoh/spaln/archive/ver.2.4.10.tar.gz && tar -xvf ver.2.4.10.tar.gz \
	&& mv spaln-ver.2.4.10 2.4.1 && rm ver.2.4.10.tar.gz \
	&& cd 2.4.1/src && ./configure --use_zlib=0 && make && make install
