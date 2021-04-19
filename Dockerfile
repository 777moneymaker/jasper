FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive 
WORKDIR /app

COPY requirements.txt .
COPY install_dependencies.sh .
# COPY jasper/ .
# COPY tests/ .
COPY jasper-vh .

RUN apt-get update \
 && apt-get install -y python3-pip python3-dev \
 && cd /usr/local/bin \
 && ln -s /usr/bin/python3 python \
 && pip3 install --upgrade pip

RUN python3 -m pip install -r requirements.txt \
 && python3 -m pip install jasper-vh \
 && bash install_dependencies.sh

RUN rm requirements.txt install_dependencies.sh \
 && sudo apt autoremove -y

ENTRYPOINT ["./jasper-vh"]