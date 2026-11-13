# =============================================================================
# gencorr-ldsc: LDSC for Python 3 (CBIIT/ldsc, branch ldsc313, pinned commit).
#
# CBIIT/ldsc is a port of LDSC 1.0.1 (Bulik-Sullivan et al.) to Python 3 with
# current numpy/pandas/scipy; it reports "Version 3.0.1" in its logs. The
# statistics are the original ones: the port changes Python 2 idioms and
# pandas/numpy API calls, not the regressions (see docs/ldsc_py3.md).
#
#   /opt/ldsc            CBIIT/ldsc checkout (ldsc.py, munge_sumstats.py, ldscore/)
#   /opt/ldsc/test       upstream unit tests; run with pytest at build time
#                        (run_upstream_tests.sh)
#   ldsc.py, munge_sumstats.py on PATH
#
# Every version is pinned:
#   - python:3.13.7-slim-bookworm
#   - LDSC from an exact Git commit (not on PyPI/bioconda)
#   - numpy/pandas/scipy as in the ldsc313 branch's requirements.txt;
#     bitarray 2.9.3 instead of 2.6.0 (no Python 3.13 wheel for 2.6.0; 2.x API
#     is unchanged, used only for computing LD scores from PLINK files)
#   - pytest + pynose: the upstream tests import nose.tools, and the original
#     nose does not import on Python >= 3.12; pynose is its maintained fork.
#
# The build fails if any upstream unit test fails.
#
# Build (from the repository root):
#   docker build -t ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996 containers/gencorr-ldsc
# Run the upstream tests again in the built image:
#   docker run --rm ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996 run_upstream_tests.sh
#
# LICENCE: LDSC is GPL-3.0, so this image may be public.
# =============================================================================
FROM python:3.13.7-slim-bookworm AS ldsc

LABEL org.opencontainers.image.title="gencorr-ldsc" \
      org.opencontainers.image.description="LDSC (CBIIT Python 3 port) for the mauricioroza/gencorr pipeline" \
      org.opencontainers.image.source="https://github.com/mauricioroza/gencorr" \
      org.opencontainers.image.licenses="GPL-3.0"

ARG LDSC_REPO=https://github.com/CBIIT/ldsc
ARG LDSC_COMMIT=128f99674d5cfb13eaaf60bc889362b4d3b585b4

# procps: Nextflow needs `ps` for task metrics; git only to fetch LDSC
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps git \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        numpy==2.4.3 pandas==2.3.3 scipy==1.15.0 bitarray==2.9.3 \
        pytest==8.4.2 pynose==1.5.5

RUN git clone --quiet "${LDSC_REPO}" /opt/ldsc \
 && cd /opt/ldsc \
 && git checkout --quiet "${LDSC_COMMIT}" \
 && rm -rf .git \
 && chmod +x ldsc.py munge_sumstats.py \
 && ln -s /opt/ldsc/ldsc.py /usr/local/bin/ldsc.py \
 && ln -s /opt/ldsc/munge_sumstats.py /usr/local/bin/munge_sumstats.py

ENV PYTHONPATH=/opt/ldsc \
    PYTHONDONTWRITEBYTECODE=1

COPY run_upstream_tests.sh /usr/local/bin/run_upstream_tests.sh
RUN chmod 755 /usr/local/bin/run_upstream_tests.sh

# `--target ldsc` stops here (untested image, for debugging)
FROM ldsc AS tested

# Upstream unit tests (harness adaptations and deselected tests are explained
# in run_upstream_tests.sh): the image is only built if they pass
RUN run_upstream_tests.sh /opt/ldsc

RUN ldsc.py -h > /dev/null && munge_sumstats.py -h > /dev/null \
 && { python -V; \
      python -c 'import numpy, pandas, scipy, bitarray; print("numpy", numpy.__version__); print("pandas", pandas.__version__); print("scipy", scipy.__version__); print("bitarray", bitarray.__version__)'; \
      echo "ldsc_commit ${LDSC_COMMIT}"; } > /opt/gencorr-ldsc.versions
