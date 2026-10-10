#
# minimum environment required to run Strelka EVS training procedures
# (somatic + germline)
#
# Build with the repository root as the docker build context so the pinned
# requirements file is available:
#
#     docker build -f docs/userGuide/trainingSomaticEmpiricalScore/Dockerfile .
#

FROM ubuntu:22.04

#
# major required packages (versions pinned in
# src/python/scoringModelTraining/requirements-training.txt):
#
# - python 3 (>= 3.10)
# - numpy
# - pandas
# - scikit-learn
#
COPY src/python/scoringModelTraining/requirements-training.txt /tmp/requirements-training.txt
RUN apt-get update -qq && \
    apt-get install -qq python3 python3-pip python3-dev build-essential zlib1g-dev && \
    pip3 install -r /tmp/requirements-training.txt

#
# Optional packages for visualization discussed
# in user guide but not required for training:
#
# - R
# - R-ggplot2
#
RUN apt-get install -qq r-base r-cran-ggplot2
