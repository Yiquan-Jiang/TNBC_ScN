"""
BRCA Spatial Transcriptomics Analysis Package
==============================================

This package provides tools for spatial transcriptomics analysis of 
breast cancer data, focusing on:

1. IL1R1+ CAF and CXCL13+ T cell colocalization
2. CDKN1A+ Macrophage and CAF spatial relationships
3. Tumor microenvironment cell-cell interactions

Modules:
--------
- spatial_colocalization_base: Base class for spatial analysis
- il1r1_cxcl13_extended_analysis: Extended IL1R1/CXCL13 analysis
- comprehensive_analysis_with_tumor_and_cdkn1a: Full analysis pipeline

Usage:
------
>>> from publication_code import ComprehensiveSpatialAnalyzer
>>> analyzer = ComprehensiveSpatialAnalyzer(data_dir=".")
>>> results = analyzer.run_comprehensive_analysis()

Author: Research Team
Date: 2024
"""

try:
    from .spatial_colocalization_base import SpatialColocalizationAnalyzer
    from .il1r1_cxcl13_extended_analysis import IL1R1_CXCL13_ColocalizationAnalyzer
    from .comprehensive_analysis_with_tumor_and_cdkn1a import ComprehensiveSpatialAnalyzer
except ImportError:
    # For direct execution
    pass

__all__ = [
    'SpatialColocalizationAnalyzer',
    'IL1R1_CXCL13_ColocalizationAnalyzer', 
    'ComprehensiveSpatialAnalyzer'
]

__version__ = '1.0.0'
__author__ = 'Research Team'

