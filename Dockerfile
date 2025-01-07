FROM ubuntu:24.04

# Install build dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential cmake

# Add a non-root user
RUN groupadd -r hyphyuser && useradd -r -g hyphyuser hyphyuser

# Create a directory for the project
WORKDIR /hyphy

# Copy project files
COPY ./cmake /hyphy/cmake
COPY ./src /hyphy/src
COPY ./contrib /hyphy/contrib
COPY ./res /hyphy/res
COPY CMakeLists.txt /hyphy
COPY ./tests /hyphy/tests

RUN chown -R hyphyuser:hyphyuser .

# Install project
RUN cmake . && make -j install

USER hyphyuser
