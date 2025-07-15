FROM /.../
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    && apt-get clean

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Anaconda в PATH
ENV PATH=$CONDA_DIR/bin:$PATH
COPY genozip_license/.genozip_license.v71 ~/.genozip_license.v71
RUN conda init
RUN conda create -n genozip -c conda-forge genozip
ENV PATH /opt/conda/envs/genozip/bin:$PATH
RUN sed -i ~/.profile -e 's/mesg n || true/tty -s \&\& mesg n/g'
RUN echo "conda activate genozip" >> ~/.bashrc

CMD ["source", "~/.profile"]
