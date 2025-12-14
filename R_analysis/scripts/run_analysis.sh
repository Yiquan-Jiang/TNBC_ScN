#!/bin/bash
# =============================================================================
# Quick Start Script for scRNA-seq Analysis
# Usage: bash run_analysis.sh
# =============================================================================

echo "============================================================================="
echo "Single-Cell RNA-seq Analysis Pipeline"
echo "TNBC: Chemotherapy-Induced Senescence Study"
echo "============================================================================="
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: R is not installed or not in PATH"
    echo "Please install R from https://www.r-project.org/"
    exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Project directory: $PROJECT_DIR"
echo ""

# Check if dependencies are installed
echo "[Step 1] Checking dependencies..."
Rscript -e 'if(!require("Seurat", quietly=TRUE)) quit(status=1)' 2>/dev/null
if [ $? -ne 0 ]; then
    echo "  Dependencies not found. Installing..."
    Rscript "$SCRIPT_DIR/00_install_dependencies.R"
else
    echo "  Dependencies OK"
fi
echo ""

# Run the main analysis
echo "[Step 2] Running main analysis..."
echo "  This may take 30-60 minutes depending on data size."
echo ""

cd "$PROJECT_DIR"
Rscript "$SCRIPT_DIR/FINAL_ANALYSIS_CODE.R"

echo ""
echo "============================================================================="
echo "Analysis complete!"
echo "Output files:"
echo "  - Figures: $PROJECT_DIR/output/figures/"
echo "  - Statistics: $PROJECT_DIR/output/statistics/"
echo "============================================================================="

