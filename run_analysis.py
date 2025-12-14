#!/usr/bin/env python3
"""
====================================================================
BRCA Spatial Transcriptomics Analysis - Main Runner Script
====================================================================
This script runs the complete spatial transcriptomics analysis pipeline.

Usage:
    python run_analysis.py [--mode MODE]
    
    Modes:
        full    - Run complete analysis (default)
        basic   - Run basic colocalization only
        il1r1   - Run IL1R1/CXCL13 analysis only
        
Author: Research Team
Date: 2024
====================================================================
"""

import argparse
import sys
import os

# Add parent directory to path for data access
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_basic_analysis(data_dir):
    """Run basic colocalization analysis."""
    try:
        from spatial_colocalization_base import SpatialColocalizationAnalyzer
    except ImportError:
        from publication_code.spatial_colocalization_base import SpatialColocalizationAnalyzer
    
    print("=" * 70)
    print("Running Basic Spatial Colocalization Analysis")
    print("=" * 70)
    
    analyzer = SpatialColocalizationAnalyzer(data_dir)
    
    analyzer.load_spatial_data()
    if not analyzer.spatial_data:
        print("No data loaded. Exiting.")
        return None
    
    analyzer.annotate_cell_types()
    
    # Colocalization analysis
    coloc_results = analyzer.calculate_colocalization('IL1R1_pos_Fibroblasts', 'T_cells')
    analyzer.results['IL1R1_Fib_vs_T_cells_coloc'] = coloc_results
    
    # Correlation analysis
    corr_results = analyzer.calculate_correlation('IL1R1_pos_Fibroblasts', 'T_cells')
    analyzer.results['IL1R1_Fib_vs_T_cells_corr'] = corr_results
    
    # Visualizations
    if coloc_results:
        analyzer.plot_colocalization_results(coloc_results, 'IL1R1+ Fibroblasts vs T cells')
    
    if corr_results:
        analyzer.plot_correlation_results(corr_results, 'IL1R1+ Fibroblasts vs T cells')
    
    analyzer.save_results_summary()
    
    print("\nBasic analysis complete!")
    return analyzer.results


def run_il1r1_cxcl13_analysis(data_dir):
    """Run extended IL1R1/CXCL13 analysis."""
    try:
        from il1r1_cxcl13_extended_analysis import IL1R1_CXCL13_ColocalizationAnalyzer
    except ImportError:
        from publication_code.il1r1_cxcl13_extended_analysis import IL1R1_CXCL13_ColocalizationAnalyzer
    
    print("=" * 70)
    print("Running IL1R1+ CAF vs CXCL13+ T Cell Extended Analysis")
    print("=" * 70)
    
    analyzer = IL1R1_CXCL13_ColocalizationAnalyzer(data_dir)
    results = analyzer.run_extended_analysis()
    
    return results


def run_comprehensive_analysis(data_dir):
    """Run complete comprehensive analysis."""
    try:
        from comprehensive_analysis_with_tumor_and_cdkn1a import ComprehensiveSpatialAnalyzer
    except ImportError:
        from publication_code.comprehensive_analysis_with_tumor_and_cdkn1a import ComprehensiveSpatialAnalyzer
    
    print("=" * 70)
    print("Running Comprehensive Spatial Transcriptomics Analysis")
    print("=" * 70)
    
    analyzer = ComprehensiveSpatialAnalyzer(data_dir)
    results = analyzer.run_comprehensive_analysis()
    
    return results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='BRCA Spatial Transcriptomics Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python run_analysis.py                    # Run full analysis
    python run_analysis.py --mode basic       # Run basic analysis only
    python run_analysis.py --mode il1r1       # Run IL1R1/CXCL13 analysis only
    python run_analysis.py --data /path/to/data  # Specify data directory
        """
    )
    
    parser.add_argument(
        '--mode', 
        choices=['full', 'basic', 'il1r1'],
        default='full',
        help='Analysis mode (default: full)'
    )
    
    parser.add_argument(
        '--data',
        default='..',
        help='Path to data directory (default: parent directory)'
    )
    
    args = parser.parse_args()
    
    print("\n" + "=" * 70)
    print("BRCA Spatial Transcriptomics Analysis Pipeline")
    print("=" * 70)
    print(f"\nMode: {args.mode}")
    print(f"Data directory: {args.data}")
    print("")
    
    # Run appropriate analysis
    if args.mode == 'basic':
        results = run_basic_analysis(args.data)
    elif args.mode == 'il1r1':
        results = run_il1r1_cxcl13_analysis(args.data)
    else:  # full
        results = run_comprehensive_analysis(args.data)
    
    if results:
        print("\n" + "=" * 70)
        print("Analysis completed successfully!")
        print("=" * 70)
        print("\nCheck the output files for results.")
    else:
        print("\nAnalysis failed or no results generated.")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

