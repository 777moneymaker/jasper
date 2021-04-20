FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive 
WORKDIR /app

COPY requirements.txt .
COPY install_dependencies.sh .
COPY jasper-vh .

RUN apt-get update \
 && bash install_dependencies.sh \
 && apt-get install -y python3-pip python3-dev \
 && cd /usr/local/bin \ 
 && ln -s /usr/bin/python3 python \
 && pip install --upgrade pip \
 && apt remove -y cmake autoconf g++ wget git \
 && apt-get clean && apt autoremove -y


RUN pip install -r requirements.txt \
 && pip install jasper-vh \
 && rm requirements.txt 

ENTRYPOINT ["./jasper-vh"]