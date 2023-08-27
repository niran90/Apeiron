FROM ubuntu:22.04

WORKDIR /project

# Install dependencies/packages.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get -y --no-install-recommends install \
        build-essential \
        git \
        clang && \
    rm -rf /var/lib/apt/lists/*

# Copy all source files into container.
COPY . .

# Compile C++ source code.
RUN clang++ -std=c++20 -o hello src/main.cpp && \
    chmod a+x hello  

CMD ["./hello"]