#!/usr/bin/env python3
"""
====================================================================
Breast Cancer Spatial Transcriptomics Analysis - Base Module
====================================================================
This module contains the base analyzer class for spatial colocalization
analysis of IL1R1+ Fibroblasts and T cells in breast cancer.

Author: Research Team
Date: 2024
Version: 1.0

Description:
    This script provides foundational functionality for:
    1. Loading 10X Visium spatial transcriptomics data
    2. Cell type annotation using marker gene signatures
    3. Basic spatial colocalization analysis
    4. Correlation analysis between cell types
====================================================================
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import warnings
from pathlib import Path
import glob
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
import matplotlib.patches as patches
from matplotlib.patches import Circle
import matplotlib.colors as mcolors

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

# Configure matplotlib for publication-quality PDF output
plt.rcParams.update({
    'pdf.fonttype': 42,      # Embed fonts as Type 42 (editable text)
    'ps.fonttype': 42,       # Ensure PostScript compatibility
    'font.family': 'Arial',
    'font.size': 12,
    'axes.titlesize': 25,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20,
    'figure.titlesize': 25
})


class SpatialColocalizationAnalyzer:
    """
    Base class for spatial transcriptomics colocalization analysis.
    
    This analyzer provides methods for:
    - Loading spatial transcriptomics data
    - Cell type annotation using marker genes
    - Spatial colocalization analysis
    - Correlation analysis between cell types
    
    Attributes:
        data_dir (Path): Directory containing spatial data
        spatial_data (dict): Dictionary of sample AnnData objects
        sc_reference (any): Single-cell reference data (optional)
        results (dict): Analysis results
    """
    
    def __init__(self, data_dir="."):
        """
        Initialize the SpatialColocalizationAnalyzer.
        
        Args:
            data_dir (str): Path to directory containing spatial data
        """
        self.data_dir = Path(data_dir)
        self.spatial_data = {}
        self.sc_reference = None
        self.results = {}
        
    def _create_cell_type_markers(self):
        """
        Create cell type marker gene dictionary for annotation.
        
        Returns:
            dict: Dictionary mapping cell types to marker genes
        """
        cell_type_markers = {
            'Fibroblasts': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA', 'IL1R1'],
            'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CXCL13'],
            'B_cells': ['CD19', 'CD20', 'MS4A1', 'CD79A'],
            'Macrophages': ['CD68', 'CD163', 'MARCO', 'CDKN1A'],
            'Epithelial': ['EPCAM', 'KRT18', 'KRT19'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5']
        }
        return cell_type_markers
    
    def load_spatial_data(self):
        """
        Load all spatial transcriptomics data from the data directory.
        
        Returns:
            dict: Dictionary of AnnData objects keyed by sample ID
        """
        print("Loading spatial transcriptomics data...")
        
        # Find all spatial data directories
        spatial_dirs = glob.glob(str(self.data_dir / "download_page_*_spatial"))
        
        for spatial_dir in spatial_dirs:
            sample_id = Path(spatial_dir).name.replace("download_page_", "").replace("_spatial", "")
            print(f"Loading sample: {sample_id}")
            
            try:
                # Find h5 file
                h5_files = glob.glob(f"{spatial_dir}/*.h5")
                if h5_files:
                    h5_file = os.path.basename(h5_files[0])
                    # Load 10x Visium data
                    adata = sc.read_visium(
                        path=spatial_dir,
                        count_file=h5_file,
                        load_images=True
                    )
                    
                    # Basic preprocessing
                    adata.var_names_make_unique()
                    adata.var['mt'] = adata.var_names.str.startswith('MT-')
                    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
                    
                    # Add sample information
                    adata.obs['sample_id'] = sample_id
                    
                    self.spatial_data[sample_id] = adata
                    print(f"  Sample {sample_id} loaded: {adata.shape}")
                    
            except Exception as e:
                print(f"  Failed to load sample {sample_id}: {e}")
                
        print(f"Total samples loaded: {len(self.spatial_data)}")
        return self.spatial_data
    
    def annotate_cell_types(self):
        """
        Annotate cell types in spatial transcriptomics data using marker genes.
        """
        print("Annotating cell types...")
        
        cell_type_markers = self._create_cell_type_markers()
        
        for sample_id, adata in self.spatial_data.items():
            print(f"Annotating sample: {sample_id}")
            
            # Normalize and log-transform
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # Calculate score for each cell type
            for cell_type, markers in cell_type_markers.items():
                # Use only available markers
                available_markers = [gene for gene in markers if gene in adata.var_names]
                
                if available_markers:
                    # Calculate cell type score (mean expression of markers)
                    adata.obs[f'{cell_type}_score'] = adata.X[:, [adata.var_names.get_loc(gene) 
                                                                for gene in available_markers]].mean(axis=1)
                else:
                    adata.obs[f'{cell_type}_score'] = 0
            
            # Assign primary cell type
            cell_types = ['Fibroblasts', 'T_cells', 'B_cells', 'Macrophages', 'Epithelial', 'Endothelial']
            score_columns = [f'{ct}_score' for ct in cell_types]
            
            # Find highest scoring cell type for each spot
            max_scores = adata.obs[score_columns].max(axis=1)
            adata.obs['cell_type'] = adata.obs[score_columns].idxmax(axis=1).str.replace('_score', '')
            
            # Set threshold for unknown cells
            threshold = 0.1
            adata.obs.loc[max_scores < threshold, 'cell_type'] = 'Unknown'
            
            # Define special phenotypes
            # IL1R1+ Fibroblasts
            fibroblast_mask = adata.obs['cell_type'] == 'Fibroblasts'
            if 'IL1R1' in adata.var_names:
                il1r1_high = adata.X[:, adata.var_names.get_loc('IL1R1')] > np.percentile(
                    adata.X[:, adata.var_names.get_loc('IL1R1')], 75)
                adata.obs['IL1R1_pos_Fibroblasts'] = fibroblast_mask & il1r1_high
            else:
                adata.obs['IL1R1_pos_Fibroblasts'] = False
            
            # CXCL13+ T cells
            t_cell_mask = adata.obs['cell_type'] == 'T_cells'
            if 'CXCL13' in adata.var_names:
                cxcl13_high = adata.X[:, adata.var_names.get_loc('CXCL13')] > np.percentile(
                    adata.X[:, adata.var_names.get_loc('CXCL13')], 75)
                adata.obs['CXCL13_pos_T_cells'] = t_cell_mask & cxcl13_high
            else:
                adata.obs['CXCL13_pos_T_cells'] = False
            
            print(f"  Sample {sample_id} annotation complete")
            print(f"  Cell type distribution:\n{adata.obs['cell_type'].value_counts()}")
    
    def calculate_colocalization(self, cell_type1, cell_type2, radius=50):
        """
        Calculate colocalization index between two cell types.
        
        Args:
            cell_type1 (str): First cell type identifier
            cell_type2 (str): Second cell type identifier
            radius (float): Distance threshold for neighbor detection
            
        Returns:
            dict: Colocalization results per sample
        """
        print(f"Calculating colocalization: {cell_type1} vs {cell_type2}...")
        
        colocalization_results = {}
        
        for sample_id, adata in self.spatial_data.items():
            print(f"  Analyzing sample: {sample_id}")
            
            # Get spatial coordinates
            coords = adata.obsm['spatial']
            
            # Get cell type masks
            if cell_type1 in adata.obs.columns:
                cells1_mask = adata.obs[cell_type1].astype(bool)
            else:
                cells1_mask = adata.obs['cell_type'] == cell_type1
                
            if cell_type2 in adata.obs.columns:
                cells2_mask = adata.obs[cell_type2].astype(bool)
            else:
                cells2_mask = adata.obs['cell_type'] == cell_type2
            
            cells1_coords = coords[cells1_mask]
            cells2_coords = coords[cells2_mask]
            
            if len(cells1_coords) == 0 or len(cells2_coords) == 0:
                print(f"    Sample {sample_id}: missing target cell types")
                continue
            
            # Calculate distance matrix
            distances = pairwise_distances(cells1_coords, cells2_coords)
            
            # Count neighbors within radius
            neighbors_within_radius = (distances <= radius).sum(axis=1)
            
            # Calculate colocalization index
            colocalization_index = neighbors_within_radius.mean()
            
            # Permutation test for random expectation
            n_permutations = 1000
            random_indices = []
            
            for _ in range(n_permutations):
                shuffled_coords2 = cells2_coords[np.random.permutation(len(cells2_coords))]
                random_distances = pairwise_distances(cells1_coords, shuffled_coords2)
                random_neighbors = (random_distances <= radius).sum(axis=1)
                random_indices.append(random_neighbors.mean())
            
            random_mean = np.mean(random_indices)
            p_value = np.sum(np.array(random_indices) >= colocalization_index) / n_permutations
            
            # Normalized colocalization index
            normalized_colocalization = colocalization_index / random_mean if random_mean > 0 else 0
            
            colocalization_results[sample_id] = {
                'colocalization_index': colocalization_index,
                'normalized_colocalization': normalized_colocalization,
                'random_mean': random_mean,
                'p_value': p_value,
                'n_cells1': len(cells1_coords),
                'n_cells2': len(cells2_coords)
            }
        
        return colocalization_results
    
    def calculate_correlation(self, cell_type1, cell_type2):
        """
        Calculate correlation between two cell type counts across samples.
        
        Args:
            cell_type1 (str): First cell type identifier
            cell_type2 (str): Second cell type identifier
            
        Returns:
            dict: Correlation results
        """
        print(f"Calculating correlation: {cell_type1} vs {cell_type2}...")
        
        counts1 = []
        counts2 = []
        sample_ids = []
        
        for sample_id, adata in self.spatial_data.items():
            if cell_type1 in adata.obs.columns:
                count1 = adata.obs[cell_type1].sum()
            else:
                count1 = (adata.obs['cell_type'] == cell_type1).sum()
                
            if cell_type2 in adata.obs.columns:
                count2 = adata.obs[cell_type2].sum()
            else:
                count2 = (adata.obs['cell_type'] == cell_type2).sum()
            
            counts1.append(count1)
            counts2.append(count2)
            sample_ids.append(sample_id)
        
        # Calculate correlations
        if len(counts1) > 2:
            pearson_r, pearson_p = pearsonr(counts1, counts2)
            spearman_r, spearman_p = spearmanr(counts1, counts2)
        else:
            pearson_r, pearson_p = np.nan, np.nan
            spearman_r, spearman_p = np.nan, np.nan
        
        return {
            'counts1': counts1,
            'counts2': counts2,
            'sample_ids': sample_ids,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p
        }
    
    def visualize_spatial_distribution(self, sample_id, cell_types_to_show):
        """
        Visualize spatial distribution of specified cell types.
        
        Args:
            sample_id (str): Sample identifier
            cell_types_to_show (list): List of cell type identifiers
        """
        if sample_id not in self.spatial_data:
            print(f"Sample {sample_id} not found")
            return
        
        adata = self.spatial_data[sample_id]
        
        # Create subplots
        n_types = len(cell_types_to_show)
        fig, axes = plt.subplots(1, n_types + 1, figsize=(6 * (n_types + 1), 6))
        if n_types == 0:
            axes = [axes]
        
        # Show tissue section
        if hasattr(adata, 'uns') and 'spatial' in adata.uns:
            img = adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['images']['hires']
            axes[0].imshow(img)
            axes[0].set_title('Tissue Section', fontsize=23)
            axes[0].axis('off')
        
        # Visualize each cell type
        for i, cell_type in enumerate(cell_types_to_show):
            ax = axes[i + 1] if n_types > 0 else axes[0]
            
            coords = adata.obsm['spatial']
            
            if cell_type in adata.obs.columns:
                mask = adata.obs[cell_type].astype(bool)
                title = cell_type.replace('_', ' ')
            else:
                mask = adata.obs['cell_type'] == cell_type
                title = cell_type
            
            # Plot all spots (gray)
            ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=10, alpha=0.5)
            
            # Plot target cell type (colored)
            if mask.sum() > 0:
                target_coords = coords[mask]
                ax.scatter(target_coords[:, 0], target_coords[:, 1], 
                          c='red', s=20, alpha=0.8, edgecolors='darkred', linewidth=0.5)
            
            ax.set_title(f'{title}\n(n={mask.sum()})', fontsize=23)
            ax.set_xlabel('Spatial X', fontsize=20)
            ax.set_ylabel('Spatial Y', fontsize=20)
            ax.tick_params(labelsize=16)
        
        plt.tight_layout()
        plt.savefig(f'spatial_distribution_{sample_id}.pdf', dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig(f'spatial_distribution_{sample_id}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def plot_colocalization_results(self, results, title):
        """
        Visualize colocalization analysis results.
        
        Args:
            results (dict): Colocalization results dictionary
            title (str): Plot title
        """
        samples = list(results.keys())
        colocalization_indices = [results[s]['normalized_colocalization'] for s in samples]
        p_values = [results[s]['p_value'] for s in samples]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Colocalization index bar chart
        colors = ['red' if p < 0.05 else 'gray' for p in p_values]
        bars = ax1.bar(range(len(samples)), colocalization_indices, color=colors, alpha=0.7)
        ax1.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Random expectation')
        ax1.set_xlabel('Sample ID', fontsize=20)
        ax1.set_ylabel('Normalized Colocalization Index', fontsize=20)
        ax1.set_title(f'{title}\nColocalization Analysis', fontsize=25)
        ax1.set_xticks(range(len(samples)))
        ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=16)
        ax1.tick_params(labelsize=16)
        ax1.legend(fontsize=20)
        
        # Add significance annotations
        for i, (bar, p_val) in enumerate(zip(bars, p_values)):
            height = bar.get_height()
            significance = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                    significance, ha='center', va='bottom', fontsize=16)
        
        # P-value distribution
        ax2.hist(p_values, bins=10, alpha=0.7, color='steelblue', edgecolor='black')
        ax2.axvline(x=0.05, color='red', linestyle='--', alpha=0.7, label='p = 0.05')
        ax2.set_xlabel('P-value', fontsize=20)
        ax2.set_ylabel('Frequency', fontsize=20)
        ax2.set_title('P-value Distribution', fontsize=23)
        ax2.tick_params(labelsize=16)
        ax2.legend(fontsize=20)
        
        plt.tight_layout()
        filename = title.replace(' ', '_').replace('+', 'pos').lower()
        plt.savefig(f'colocalization_{filename}.pdf', dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig(f'colocalization_{filename}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def plot_correlation_results(self, results, title):
        """
        Visualize correlation analysis results.
        
        Args:
            results (dict): Correlation results dictionary
            title (str): Plot title
        """
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        
        counts1 = results['counts1']
        counts2 = results['counts2']
        sample_ids = results['sample_ids']
        
        # Scatter plot
        ax.scatter(counts1, counts2, s=100, alpha=0.7, color='steelblue', edgecolors='darkblue')
        
        # Add sample labels
        for i, sample_id in enumerate(sample_ids):
            ax.annotate(sample_id, (counts1[i], counts2[i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=12)
        
        # Regression line
        if len(counts1) > 1:
            z = np.polyfit(counts1, counts2, 1)
            p = np.poly1d(z)
            ax.plot(counts1, p(counts1), "r--", alpha=0.8, linewidth=2)
        
        # Statistics text
        stats_text = f'Pearson r = {results["pearson_r"]:.3f} (p = {results["pearson_p"]:.3f})\n'
        stats_text += f'Spearman ρ = {results["spearman_r"]:.3f} (p = {results["spearman_p"]:.3f})'
        
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                verticalalignment='top', fontsize=16)
        
        ax.set_xlabel(title.split(' vs ')[0] + ' Count', fontsize=20)
        ax.set_ylabel(title.split(' vs ')[1] + ' Count', fontsize=20)
        ax.set_title(f'{title}\nCorrelation Analysis', fontsize=25)
        ax.tick_params(labelsize=16)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = title.replace(' ', '_').replace('+', 'pos').lower()
        plt.savefig(f'correlation_{filename}.pdf', dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig(f'correlation_{filename}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def save_results_summary(self):
        """Save analysis results summary to CSV."""
        summary = {
            'analysis_type': 'Spatial Colocalization and Correlation Analysis',
            'samples_analyzed': list(self.spatial_data.keys()),
            'n_samples': len(self.spatial_data)
        }
        
        # Add colocalization results summary
        for key, results in self.results.items():
            if 'coloc' in key:
                summary[f'{key}_significant_samples'] = sum(
                    1 for r in results.values() if r['p_value'] < 0.05
                )
                summary[f'{key}_mean_colocalization'] = np.mean([
                    r['normalized_colocalization'] for r in results.values()
                ])
            elif 'corr' in key:
                summary[f'{key}_pearson_r'] = results['pearson_r']
                summary[f'{key}_pearson_p'] = results['pearson_p']
                summary[f'{key}_spearman_r'] = results['spearman_r']
                summary[f'{key}_spearman_p'] = results['spearman_p']
        
        pd.DataFrame([summary]).to_csv('analysis_summary.csv', index=False)
        print("Results summary saved to analysis_summary.csv")


if __name__ == "__main__":
    print("=" * 70)
    print("BRCA Spatial Transcriptomics - Base Colocalization Analysis")
    print("=" * 70)
    
    analyzer = SpatialColocalizationAnalyzer()
    
    # Load data
    analyzer.load_spatial_data()
    
    if not analyzer.spatial_data:
        print("No spatial data loaded. Exiting.")
        exit(1)
    
    # Annotate cell types
    analyzer.annotate_cell_types()
    
    # Run colocalization analysis
    print("\n=== Colocalization Analysis ===")
    coloc_results = analyzer.calculate_colocalization('IL1R1_pos_Fibroblasts', 'T_cells')
    analyzer.results['IL1R1_Fib_vs_T_cells_coloc'] = coloc_results
    
    # Run correlation analysis
    print("\n=== Correlation Analysis ===")
    corr_results = analyzer.calculate_correlation('IL1R1_pos_Fibroblasts', 'T_cells')
    analyzer.results['IL1R1_Fib_vs_T_cells_corr'] = corr_results
    
    # Generate visualizations
    print("\n=== Generating Visualizations ===")
    if coloc_results:
        analyzer.plot_colocalization_results(coloc_results, 'IL1R1+ Fibroblasts vs T cells')
    if corr_results:
        analyzer.plot_correlation_results(corr_results, 'IL1R1+ Fibroblasts vs T cells')
    
    # Save results
    analyzer.save_results_summary()
    
    print("\nBase analysis complete!")

