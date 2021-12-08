FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/esga-1.3/bin:/opt/bin:/opt/spaln/2.4.6/bin:/opt/conda/envs/esga-1.3/opt/pasa-2.4.1:$PATH
ENV EVM_HOME /opt/conda/envs/esga-1.3/opt/evidencemodeler-1.1.1
ENV PASAHOME /opt/conda/envs/esga-1.3/opt/pasa-2.4.1
RUN mkdir -p /opt/bin && cd /opt/bin && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords && chmod +x faSomeRecords
RUN apt-get -y update && apt-get -y install make gcc g++ zlib1g-dev zlib1g
RUN mkdir -p /opt/spaln && cd /opt/spaln && wget https://github.com/ogotoh/spaln/archive/refs/tags/ver.2.4.6.tar.gz && tar -xvf ver.2.4.6.tar.gz \
	&& mv spaln-ver.2.4.6 2.4.6 && rm ver.2.4.6.tar.gz \
	&& cd 2.4.6/src && ./configure --use_zlib=1 && make && make install
RUN cpan -i URI::Encode

