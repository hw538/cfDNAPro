# Reuse the existing R/Bioconductor dependency environment at an immutable digest.
# cfDNAPro itself is always replaced with the source in this checkout.
FROM zetianjia/cfdnapro:1.7.3@sha256:d5c0bbbc6ad411b52e7afd11bd2e6ff70ec4b1f04bd60672aecf8495a2e53e1f

USER root
ARG VERSION
ARG REVISION
LABEL org.opencontainers.image.title="cfDNAPro" \
      org.opencontainers.image.description="cfDNA fragmentomic feature extraction" \
      org.opencontainers.image.source="https://github.com/hw538/cfDNAPro" \
      org.opencontainers.image.url="https://github.com/hw538/cfDNAPro" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${REVISION}"

COPY . /tmp/cfDNAPro
RUN micromamba run -n base Rscript --vanilla -e \
      'if (!requireNamespace("testthat", quietly=TRUE) || packageVersion("testthat") < "3.2.0") install.packages("testthat", repos="https://cloud.r-project.org"); stopifnot(packageVersion("testthat") >= "3.2.0")' \
    && micromamba run -n base R CMD INSTALL --install-tests /tmp/cfDNAPro \
    && micromamba run -n base Rscript --vanilla -e \
      'expected <- commandArgs(TRUE)[1]; stopifnot(nzchar(expected), as.character(packageVersion("cfDNAPro")) == expected); cat("Installed cfDNAPro", expected, "\n")' "${VERSION}" \
    && rm -rf /tmp/cfDNAPro

USER mambauser
WORKDIR /tmp
CMD ["R", "--no-save"]
