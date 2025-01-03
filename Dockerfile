FROM alpine:3.21.0

# Install build dependencies
RUN apk add --no-cache build-base cmake

# Create a directory for the project
WORKDIR /hyphy

# Copy project files
COPY . /hyphy

# Install project
RUN cmake . && make install
