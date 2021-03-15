FROM ubuntu:bionic

RUN apt update && apt install -y cmake libnetcdf-dev libcgal-dev libboost-dev
##WORKDIR /home/pbarletta/labo/ANA2
##COPY . .
#RUN 

