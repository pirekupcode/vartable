FROM centos:centos7
MAINTAINER Michael Panciera

ENV PYTHON_VERSION 3.6

RUN yum -y update && \
    yum -y install curl bzip2

ADD . /vartable

WORKDIR /vartable

RUN bash install.sh /vartable/miniconda

ENV PATH=/vartable/miniconda/bin/:$PATH

RUN conda clean --all --yes && \ 
    rm miniconda3.sh && \
    rmp -e --nodeps curl bzip2 && \ 
    yum clean all # this inherited image should `yum clean all` automatically
