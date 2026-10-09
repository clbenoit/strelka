FROM python:3.13-slim

# install system dependencies
# GNU parallel 20161222: newer releases truncate the --header/--pipe header
# for jobs after the first one, which breaks --parallel mode on VCFs with
# large headers
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget make perl && \
    wget -q https://ftp.gnu.org/gnu/parallel/parallel-20161222.tar.bz2 && \
    tar xjf parallel-20161222.tar.bz2 && \
    cd parallel-20161222 && \
    ./configure --prefix=/usr/local && make && make install && \
    cd / && rm -rf parallel-20161222 parallel-20161222.tar.bz2 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# update python package manager
RUN pip install -U pip setuptools

WORKDIR /opt/denovo

# install runtime dependencies (avoid cache invalidation)
COPY requirements.txt  /opt/denovo/
RUN pip install -U -r requirements.txt

# copy package content
COPY . .

# install test dependencies and run the tests
RUN pip install -U -r test/test_requirements.txt
RUN python -m pytest

# define bash as entry point
CMD ["/usr/bin/env", "bash"]
