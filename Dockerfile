FROM centos:centos7
MAINTAINER Michael Panciera

ENV PYTHON_VERSION 3.6

RUN yum -y update && \
    yum -y install \
      cmake patch git curl

RUN yum install -y bzip2

COPY july-22-latest-3.sh /tmp/miniconda.sh

RUN bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh \
    && conda update conda \
    && conda install -y python=$PYTHON_VERSION \
    && conda clean --all --yes \
    && conda install -c conda-requirements.txt \
    && pip install pip-requirements.txt \
    && python setup.py develop
