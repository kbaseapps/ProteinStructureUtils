FROM kbase/sdkbase2:python
MAINTAINER KBase Developer

# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN pip install --upgrade pip \
    && pip install biopython --upgrade \
    && pip install --upgrade requests \
    && pip install Cython \
    && pip install wheel \
    && pip install pandas \
    && pip install mock==4.0.3\
    && pip install openpyxl \
    && pip install python-graphql-client

RUN apt-get update \
    && apt-get -y install wget \
    && apt-get -y install libgomp1 \
    && apt-get install -y gcc
#RUN rm -rf /miniconda/lib/python3.6/site-packages/numpy
#RUN rm -rf /miniconda/lib/python3.6/site-packages/ruamel*
#RUN pip install networkx
#RUN pip install --use-deprecated=legacy-resolver git+https://github.com/ModelSEED/ModelSEEDpy.git@28de76a315c607944c992ce21984f303b2cf4e4b


ENV BLAST_VERSION='2.13.0'

RUN cd /opt \
    && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && tar zxvpf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

ENV PATH $PATH:/opt/ncbi-blast-${BLAST_VERSION}+/bin
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
