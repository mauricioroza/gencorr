# =============================================================================
# gencorr-ldak: LDAK 6 for the SumHer method (ldak6 --sum-cors).
#
# LDAK is GPL-3.0 (https://github.com/dougspeed/LDAK/blob/main/LICENCE.txt), so
# this image may be published (unlike gencorr-tools, which contains LAVA).
# The version is pinned to the same package as the module's environment.yml
# (genomedk::ldak6=6.1), so conda and container runs use the same binary.
#
# Build (from the repository root):
#   docker build -t ghcr.io/mauricioroza/gencorr-ldak:6.1 containers/gencorr-ldak
# Published: docker pull ghcr.io/mauricioroza/gencorr-ldak:6.1
# =============================================================================
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="gencorr-ldak" \
      org.opencontainers.image.description="LDAK 6.1 (SumHer) for the mauricioroza/gencorr pipeline" \
      org.opencontainers.image.source="https://github.com/mauricioroza/gencorr" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

# procps: Nextflow needs `ps` for task metrics
RUN micromamba install -y -n base -c genomedk -c conda-forge ldak6=6.1 procps-ng \
 && micromamba clean -a -y

ENV PATH="/opt/conda/bin:${PATH}"
