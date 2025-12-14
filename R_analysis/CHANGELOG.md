# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-12-14

### Added
- Initial release of the scRNA-seq analysis pipeline
- Complete workflow from raw data to publication figures
- Support for 10X Genomics data format
- Harmony batch correction integration
- Cell type annotation using canonical markers
- Hypothesis testing for:
  - Myeloid cell senescence
  - Fibroblast proportion changes
  - CXCL13+ T cell proportion
- CD8+ T cell functional analysis (cytotoxicity/exhaustion)
- Fibroblast CAF activation analysis
- Macrophage SASP analysis
- Publication-quality visualization (12 panels)
- Statistical summary export

### Technical Details
- **DotPlot Fix**: Resolved color mapping issues by using `group.by` instead of `split.by` and setting `scale = FALSE`
- **Color Palette**: Implemented Morandi color scheme for consistent aesthetics
- **Reproducibility**: Added `renv.lock` for exact package version tracking

## [Unreleased]

### Planned
- Integration with CellChat for cell-cell communication analysis
- Support for spatial transcriptomics data
- Interactive Shiny dashboard

