FROM ubuntu

RUN apt-get update && apt-get install -y \
#libhdf5-serial-dev \
#g++ \
#libcurl4 \
#libcurl4-openssl-dev \
#mercurial \
#cmake \
#git \
libhdf5-dev \
libcfitsio-dev \
libccfits-dev \
libgomp1

# include the code
WORKDIR /home
COPY build/echellesimulator build/echellesimulator
COPY CMakeLists.txt CMakeLists.txt
#COPY src src
#COPY include include
COPY data/spectrographs data/spectrographs
#RUN mkdir build
WORKDIR /home/build
#RUN cmake -DCMAKE_BUILD_TYPE=Release ..
#RUN make

ENTRYPOINT ["/home/build/echellesimulator"]
CMD ["-h"]