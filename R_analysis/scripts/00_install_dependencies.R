#!/usr/bin/env Rscript
# =============================================================================
# Dependency Installation Script
# Run this script once before starting the analysis
# =============================================================================

cat("=============================================================================\n")
cat("Installing R Dependencies for scRNA-seq Analysis\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Install BiocManager (if not present)
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) {
  cat("[1/5] Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
} else {
  cat("[1/5] BiocManager already installed\n")
}

# -----------------------------------------------------------------------------
# 2. Install CRAN packages
# -----------------------------------------------------------------------------

cat("[2/5] Installing CRAN packages...\n")

cran_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "tidyr",
  "tibble",
  "stringr",
  "patchwork",
  "ggpubr",
  "ggsci",
  "harmony",
  "remotes",
  "renv"
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# -----------------------------------------------------------------------------
# 3. Install Bioconductor packages
# -----------------------------------------------------------------------------

cat("[3/5] Installing Bioconductor packages...\n")

bioc_packages <- c(
  "SingleCellExperiment",
  "scater",
  "scran"
)

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# -----------------------------------------------------------------------------
# 4. Install SCP from GitHub
# -----------------------------------------------------------------------------

cat("[4/5] Installing SCP from GitHub...\n")

if (!require("SCP", quietly = TRUE)) {
  cat("  Installing SCP (this may take a few minutes)...\n")
  remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE, quiet = TRUE)
} else {
  cat("  SCP already installed\n")
}

# -----------------------------------------------------------------------------
# 5. Verify installation
# -----------------------------------------------------------------------------

cat("[5/5] Verifying installation...\n\n")

required_packages <- c("Seurat", "SCP", "harmony", "dplyr", "ggplot2", 
                       "tidyr", "patchwork", "ggpubr", "ggsci")

all_installed <- TRUE
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("  ✓ %s (%s)\n", pkg, version))
  } else {
    cat(sprintf("  ✗ %s - FAILED\n", pkg))
    all_installed <- FALSE
  }
}

cat("\n=============================================================================\n")
if (all_installed) {
  cat("SUCCESS: All dependencies installed correctly!\n")
  cat("You can now run: source('scripts/FINAL_ANALYSIS_CODE.R')\n")
} else {
  cat("WARNING: Some packages failed to install. Please check the errors above.\n")
}
cat("=============================================================================\n")

