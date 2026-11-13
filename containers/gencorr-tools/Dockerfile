# =============================================================================
# gencorr-tools: the helper scripts' runtime and LAVA. LDSC uses the official
# biocontainer quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2; LDAK (SumHer) and
# Popcorn have their own images (containers/gencorr-ldak, gencorr-popcorn).
#
#   * Python 3 + pandas/numpy/scipy  -> bin/*.py helper scripts
#   * R + LAVA (pinned commit)       -> bin/run_lava.R
#
# Every version is pinned so an image built today and in two years contains the
# same software:
#   - rocker/r-ver:4.4.2 pins R and points CRAN at a dated Posit Package Manager
#     snapshot, so data.table / matrixsampling / cpp11 versions are frozen.
#   - LAVA is installed from an exact Git commit (it is not on CRAN/conda).
#   - Python packages are pinned with ==. Python itself is Ubuntu 24.04's
#     python3 (3.12.3); the modules' environment.yml files pin the same
#     python=3.12.3 so conda and container runs report identical versions.
#
# Build and push (from the repository root):
#   docker build -t ghcr.io/mauricioroza/gencorr-tools:0.1.0 containers/gencorr-tools
#   docker push ghcr.io/mauricioroza/gencorr-tools:0.1.0
# Apptainer on a cluster without Docker:
#   apptainer build gencorr-tools_0.1.0.sif docker://ghcr.io/mauricioroza/gencorr-tools:0.1.0
#
# LICENCE NOTE: LAVA's DESCRIPTION says "All rights reserved". Keep this image
# private (or ask the LAVA authors) until redistribution is clarified.
# =============================================================================
FROM rocker/r-ver:4.4.2

LABEL org.opencontainers.image.title="gencorr-tools" \
      org.opencontainers.image.description="Python + R/LAVA runtime for the mauricioroza/gencorr pipeline" \
      org.opencontainers.image.source="https://github.com/mauricioroza/gencorr"

ARG LAVA_COMMIT=e729a245f7b6923967a96804fbf5246eadf2d6c6

# procps: Nextflow needs `ps` for task metrics
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-pip procps git \
 && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --break-system-packages \
        pandas==2.2.3 numpy==2.2.2 scipy==1.15.1

# CRAN packages come from the snapshot rocker configured for R 4.4.2
RUN Rscript -e 'install.packages(c("data.table", "matrixsampling", "cpp11", "remotes"))' \
 && Rscript -e 'remotes::install_github("josefin-werme/LAVA", ref = Sys.getenv("LAVA_COMMIT"), upgrade = "never", dependencies = FALSE)' \
 && Rscript -e 'stopifnot(packageVersion("LAVA") == "0.1.5"); library(LAVA)'

# Record exact versions inside the image
RUN { python3 -c 'import pandas, numpy, scipy; print("pandas", pandas.__version__); print("numpy", numpy.__version__); print("scipy", scipy.__version__)'; \
      Rscript -e 'for (p in c("LAVA","data.table","matrixsampling")) cat(p, as.character(packageVersion(p)), "\n")'; \
      echo "LAVA_commit ${LAVA_COMMIT}"; } > /opt/gencorr-tools.versions
