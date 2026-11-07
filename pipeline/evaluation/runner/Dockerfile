# Agent sandbox for the skills eval. One image, shared by all three arms.
# r2u gives prebuilt CRAN binaries via apt, so R packages install in seconds.
FROM rocker/r2u:jammy

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      python3 python3-pip python3-numpy python3-pandas python3-scipy python3-statsmodels \
      python3-sklearn python3-matplotlib python3-h5py \
      plink1.9 bcftools samtools tabix less jq ripgrep \
    && ln -sf /usr/bin/plink1.9 /usr/local/bin/plink \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages the causal-genomics Skills call for (binary installs).
RUN install.r data.table dplyr tidyr readr ggplot2 Matrix remotes \
      MendelianRandomization susieR coloc mediation

# Not on CRAN: TwoSampleMR / ieugwasr (MRC-IEU r-universe), MRPRESSO (GitHub).
RUN Rscript -e 'options(repos=c(mrcieu="https://mrcieu.r-universe.dev", CRAN="https://cloud.r-project.org")); install.packages(c("ieugwasr","TwoSampleMR"))' \
 && Rscript -e 'remotes::install_github("rondolab/MR-PRESSO", upgrade="never")' || true

WORKDIR /work

# ---- second layer: tools the wider task set's Skills call for ----
RUN install.r BiocManager MicrobiomeStat || true
RUN Rscript -e 'BiocManager::install(c("limma","edgeR","DESeq2","DEqMS","ALDEx2","phyloseq","clusterProfiler"), ask=FALSE, update=FALSE)' || true

# pip last so R layers stay cached. No --break-system-packages: jammy pip is 22.x and rejects the flag.
RUN pip3 install --no-cache-dir rdkit isocor scanpy anndata
# MAGeCK is not on PyPI; try the SourceForge release, keep going if it is unreachable.
RUN pip3 install --no-cache-dir "https://downloads.sourceforge.net/project/mageck/0.5/mageck-0.5.9.5.tar.gz" || echo "MAGECK NOT INSTALLED"

# A lab machine answers to `python` and `pip`; without these every script-first Skill route costs a failed step.
RUN ln -sf /usr/bin/python3 /usr/local/bin/python && ln -sf /usr/bin/pip3 /usr/local/bin/pip
