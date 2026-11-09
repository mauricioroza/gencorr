# =============================================================================
# gencorr-popcorn: Popcorn (Brown et al. 2016) for cross-ancestry genetic
# correlation (popcorn fit), Python 3 port pinned to an exact commit.
#
# Popcorn states "GNU General Public License v3" (program banner and README;
# the repository has no separate LICENSE file), so this image may be published.
# Its requirements.txt is unpinned; the versions below are pinned so the image
# rebuilds identically. six is needed by pysnptools but not declared by it.
# numpy < 1.24: Popcorn's jackknife still uses the removed np.bool alias.
#
# Build (from the repository root):
#   docker build -t ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8 containers/gencorr-popcorn
# Published: docker pull ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8
# =============================================================================
FROM python:3.11.11-slim-bookworm

ARG POPCORN_COMMIT=20a29e8cb9c43fe88de8e903537eb874d52f38b3

LABEL org.opencontainers.image.title="gencorr-popcorn" \
      org.opencontainers.image.description="Popcorn ${POPCORN_COMMIT} for the mauricioroza/gencorr pipeline" \
      org.opencontainers.image.source="https://github.com/mauricioroza/gencorr" \
      org.opencontainers.image.licenses="GPL-3.0"

# procps: Nextflow needs `ps` for task metrics; git: install from the pinned commit
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps git \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        numpy==1.23.5 scipy==1.11.4 pandas==2.2.3 pysnptools==0.5.13 \
        bottleneck==1.4.2 statsmodels==0.14.4 six==1.16.0 \
 && pip install --no-cache-dir --no-deps "git+https://github.com/brielin/popcorn.git@${POPCORN_COMMIT}" \
 && echo "${POPCORN_COMMIT}" > /opt/popcorn.commit \
 && popcorn -h > /dev/null
