FROM continuumio/miniconda3
RUN conda update -n base -c defaults conda

RUN mkdir -m777 /app
WORKDIR /app

ADD radiomics_process.py ./
ADD docker-entrypoint.sh ./
RUN chmod 777 ./docker-entrypoint.sh

RUN conda create -n env python==3.11
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
RUN chmod -R 777 /opt/conda/envs
RUN apt-get update && apt-get install --no-install-recommends --no-install-suggests -y git build-essential cmake pigz
RUN apt-get update && apt-get install --no-install-recommends --no-install-suggests -y libsm6 libxrender-dev libxext6 ffmpeg

ADD environment.yml ./
RUN conda env create -f ./environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate $(head -1 ./environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 ./environment.yml | cut -d' ' -f2)/bin:$PATH

RUN python -m pip uninstall -y opencv-python
RUN python -m pip install opencv-python==4.5.5.64
#workaround for issues with numpy v2
RUN conda uninstall -n $(head -1 ./environment.yml | cut -d' ' -f2) numpy -y
# Install NumPy version 1.24.0
RUN conda install -n $(head -1 ./environment.yml | cut -d' ' -f2) numpy=1.24.0 -y
RUN conda install -n $(head -1 ./environment.yml | cut -d' ' -f2) pywavelets -y
RUN conda install -n $(head -1 ./environment.yml | cut -d' ' -f2) pandas -y
RUN conda install -n $(head -1 ./environment.yml | cut -d' ' -f2) scipy -y


RUN chmod -R 777 /app
WORKDIR /app

CMD ["./docker-entrypoint.sh"]