FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/esga-1.2/bin:/opt/conda/envs/esga-1.2/opt/pasa-2.4.1:/opt/bin:/opt/spaln/2.4.5/bin:$PATH
ENV PASAHOME /opt/conda/envs/esga-1.2/opt/pasa-2.4.1
ENV EVM_HOME /opt/conda/envs/esga-1.2/opt/evidencemodeler-1.1.1
RUN cp /opt/conda/envs/esga-1.2/opt/pasa-2.4.1/pasa_conf/pasa.CONFIG.template /opt/conda/envs/esga-1.2/opt/pasa-2.4.1/pasa_conf/conf.txt
RUN mkdir -p /opt/bin && cd /opt/bin && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords && chmod +x faSomeRecords
RUN apt-get -y update && apt-get -y install make gcc g++ zlib1g-dev zlib1g
RUN mkdir -p /opt/spaln && cd /opt/spaln && wget https://github.com/ogotoh/spaln/archive/refs/tags/ver2.4.5.tar.gz && tar -xvf ver2.4.5.tar.gz \
	&& mv spaln-ver2.4.5 2.4.5 && rm ver2.4.5.tar.gz \
	&& cd 2.4.5/src && ./configure --use_zlib=1 && make && make install
RUN cpan -i URI::Encode
