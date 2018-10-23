FROM ubuntu:18.04 as compile

RUN apt update \
  && apt install -y libmpich-dev libgdal-dev libgeos++-dev libshp-dev \
    libboost-program-options-dev libboost-iostreams-dev
WORKDIR /src
COPY apps ./apps/
COPY src ./src/
COPY include ./include/
COPY makefile .
RUN cd /src; mkdir obj
RUN make

FROM ubuntu:18.04
RUN apt update \
  && apt install -y --no-install-recommends mpich libmpich12 \
  libgdal20 libshp2 libboost-program-options1.65.1 libboost-iostreams1.65.1
WORKDIR /app
COPY --from=compile /src/p_points2grid .

CMD /app/p_points2grid --help
