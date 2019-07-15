FROM ubuntu

LABEL maintainer="Mathieu Fourment"

RUN apt-get update && apt-get install -y --no-install-recommends \
	build-essential \
	ca-certificates \
	cmake \
	default-jre \
	libgsl0-dev \
	libhmsbeagle1v5 \
	libhmsbeagle-java \
	python3 \
	python3-dev \
	python3-pip \
	unzip \
	wget

RUN rm /usr/lib/x86_64-linux-gnu/libhmsbeagle-opencl.so.*

# BEAST2
RUN wget https://github.com/CompEvol/beast2/releases/download/v2.5.2/BEAST.v2.5.2.Linux.tgz && tar zxvf BEAST.v2.5.2.Linux.tgz && ln -s /beast/bin/beast /usr/local/bin/beast2
# initialize local directory
RUN beast2 -help

# BWEST
RUN wget https://github.com/4ment/bwest/releases/download/v0.0.1/bwest.v0.0.1.zip && unzip bwest.v0.0.1.zip -d /root/.beast/2.5/BWEST/
RUN echo "package.path=\:/root/.beast/2.5/BEAST/lib/beast.src.jar\:/root/.beast/2.5/BEAST/lib/beast.jar\:/root/.beast/2.5/BWEST/lib/bwest.v0.0.1.jar" > /root/.beast/2.5/beauti.properties

# BEAST1
RUN wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz && tar zxvf BEASTv1.10.4.tgz && ln -s /BEASTv1.10.4/bin/beast /usr/local/bin/beast

# phylostan
RUN pip3 install --upgrade setuptools 
RUN pip3 install scons==3.0.5 phylostan

# physher
RUN wget https://github.com/4ment/physher/archive/vb-v1.0.zip && unzip vb-v1.0.zip
WORKDIR /physher-vb-v1.0/Release
RUN cmake -DBUILD_SHARED_LIBS=OFF .. && make && make install

ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/jni:/usr/lib/x86_64-linux-gnu
RUN ln -s /usr/bin/python3 /usr/bin/python 

WORKDIR /data

ENTRYPOINT ["scons"]

