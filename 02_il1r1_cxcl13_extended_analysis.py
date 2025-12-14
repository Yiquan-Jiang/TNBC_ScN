#!/usr/bin/env python3
"""
====================================================================
Breast Cancer Spatial Transcriptomics Analysis - Extended IL1R1/CXCL13 Module
====================================================================
Extended analysis module for IL1R1+ CAF and CXCL13+ T cell colocalization.

This module extends the base analyzer with:
1. Gene signature-based cell type scoring
2. IL1R1+ CAF (Cancer-Associated Fibroblast) identification
3. CXCL13+ T cell identification
4. Detailed spatial colocalization analysis
5. Spatial proximity analysis with permutation tests

Author: Research Team
Date: 2024
Version: 1.0
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
from sklearn.preprocessing import StandardScaler
from scipy import stats

# Import base module - adjust import based on how script is run
try:
    from spatial_colocalization_base import SpatialColocalizationAnalyzer
except ImportError:
    from publication_code.spatial_colocalization_base import SpatialColocalizationAnalyzer

# Suppress warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

# Configure matplotlib for publication-quality output
plt.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.family': 'Arial',
    'font.size': 12,
    'axes.titlesize': 25,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20,
    'figure.titlesize': 25
})


class IL1R1_CXCL13_ColocalizationAnalyzer(SpatialColocalizationAnalyzer):
    """
    Extended analyzer for IL1R1+ CAF and CXCL13+ T cell colocalization analysis.
    
    This class extends the base SpatialColocalizationAnalyzer with:
    - Refined gene signatures for cell type identification
    - Composite scoring for IL1R1+ CAF and CXCL13+ T cells
    - Detailed spatial proximity analysis
    - Enhanced visualization methods
    
    Attributes:
        gene_signatures (dict): Gene signatures for cell types
        combined_adata (AnnData): Combined data from all samples
        il1r1_caf_scores (dict): IL1R1+ CAF scores per sample
        cxcl13_t_scores (dict): CXCL13+ T cell scores per sample
        colocalization_results (dict): Detailed colocalization results
    """
    
    def __init__(self, data_dir="."):
        """Initialize the extended analyzer."""
        super().__init__(data_dir)
        self.gene_signatures = {}
        self.combined_adata = None
        self.il1r1_caf_scores = {}
        self.cxcl13_t_scores = {}
        self.colocalization_results = {}
        
    def define_gene_signatures(self):
        """
        Define gene signatures for cell type identification.
        
        These signatures are based on published literature and
        single-cell RNA-seq reference data.
        """
        print("Defining gene signatures...")
        
        # CAF (Cancer-Associated Fibroblast) signature
        self.gene_signatures['CAF'] = [
            'COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A1',
            'DCN', 'LUM', 'POSTN', 'FAP', 'PDGFRA', 'PDGFRB',
            'ACTA2', 'FN1', 'VIM', 'THY1', 'IL1R1'
        ]
        
        # IL1R1 specific signature
        self.gene_signatures['IL1R1_high'] = ['IL1R1']
        
        # T cell signature
        self.gene_signatures['T_cells'] = [
            'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B',
            'TRAC', 'TRBC1', 'TRBC2', 'LCK', 'ZAP70'
        ]
        
        # CXCL13+ T cell signature (TFH-like/exhausted T cells)
        self.gene_signatures['CXCL13_T'] = [
            'CXCL13', 'PDCD1', 'CTLA4', 'TIGIT', 'LAG3',
            'HAVCR2', 'TOX', 'ICOS'
        ]
        
        # Macrophage signature
        self.gene_signatures['Macrophages'] = [
            'CD68', 'CD163', 'CD14', 'MARCO', 'MSR1',
            'MRC1', 'FCGR3A', 'CSF1R'
        ]
        
        # Epithelial/Tumor cell signature
        self.gene_signatures['Epithelial'] = [
            'EPCAM', 'KRT7', 'KRT8', 'KRT18', 'KRT19',
            'CDH1', 'MUC1'
        ]
        
        # Fibroblast general signature
        self.gene_signatures['Fibroblasts'] = [
            'COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA',
            'FAP', 'THY1', 'VIM'
        ]
        
        print(f"  Defined {len(self.gene_signatures)} gene signatures")
        
    def calculate_signature_scores(self):
        """
        Calculate gene signature scores for each sample.
        
        Uses scanpy's score_genes function to calculate signature scores
        for each spot based on the defined gene signatures.
        """
        print("Calculating gene signature scores...")
        
        if not self.gene_signatures:
            self.define_gene_signatures()
        
        for sample_id, adata in self.spatial_data.items():
            print(f"  Processing sample: {sample_id}")
            
            # Calculate signature score for each gene signature
            for sig_name, genes in self.gene_signatures.items():
                # Filter to available genes
                available_genes = [g for g in genes if g in adata.var_names]
                
                if len(available_genes) >= 2:
                    try:
                        sc.tl.score_genes(adata, available_genes, 
                                         score_name=f'{sig_name}_signature_score',
                                         use_raw=False)
                    except Exception as e:
                        print(f"    Warning: Could not calculate {sig_name} signature: {e}")
                        adata.obs[f'{sig_name}_signature_score'] = 0
                elif len(available_genes) == 1:
                    # Single gene - use its expression directly
                    gene_idx = adata.var_names.get_loc(available_genes[0])
                    if hasattr(adata.X, 'toarray'):
                        adata.obs[f'{sig_name}_signature_score'] = adata.X[:, gene_idx].toarray().flatten()
                    else:
                        adata.obs[f'{sig_name}_signature_score'] = adata.X[:, gene_idx].flatten()
                else:
                    adata.obs[f'{sig_name}_signature_score'] = 0
        
        print("  Signature score calculation complete")
    
    def identify_il1r1_caf_and_cxcl13_t(self):
        """
        Identify IL1R1+ CAF and CXCL13+ T cells based on composite scores.
        
        Creates composite scores combining cell type signatures with
        marker gene expression to identify target cell populations.
        """
        print("Identifying IL1R1+ CAF and CXCL13+ T cells...")
        
        for sample_id, adata in self.spatial_data.items():
            print(f"  Processing sample: {sample_id}")
            
            # === IL1R1+ CAF Identification ===
            # Get CAF signature score
            caf_score = adata.obs.get('CAF_signature_score', 
                                     adata.obs.get('Fibroblasts_signature_score', 
                                                  pd.Series(0, index=adata.obs_names)))
            
            # Get IL1R1 expression
            if 'IL1R1' in adata.var_names:
                il1r1_idx = adata.var_names.get_loc('IL1R1')
                if hasattr(adata.X, 'toarray'):
                    il1r1_expr = adata.X[:, il1r1_idx].toarray().flatten()
                else:
                    il1r1_expr = adata.X[:, il1r1_idx].flatten()
            else:
                il1r1_expr = np.zeros(len(adata))
            
            # Normalize IL1R1 expression
            il1r1_normalized = (il1r1_expr - np.min(il1r1_expr)) / (np.max(il1r1_expr) - np.min(il1r1_expr) + 1e-8)
            
            # Calculate composite score
            il1r1_caf_composite = np.array(caf_score) * il1r1_normalized
            
            # Define IL1R1+ CAF threshold
            caf_high = np.array(caf_score) > np.percentile(caf_score, 75)
            il1r1_high = il1r1_expr > np.percentile(il1r1_expr, 75)
            il1r1_caf_mask = caf_high & il1r1_high
            
            # Save to adata
            adata.obs['IL1R1_CAF_composite_score'] = il1r1_caf_composite
            adata.obs['IL1R1_CAF_mask'] = il1r1_caf_mask
            
            # Save to analyzer
            self.il1r1_caf_scores[sample_id] = {
                'composite_score': il1r1_caf_composite,
                'mask': il1r1_caf_mask,
                'n_il1r1_caf': il1r1_caf_mask.sum()
            }
            
            # === CXCL13+ T Cell Identification ===
            # Get T cell signature score
            t_cell_score = adata.obs.get('T_cells_signature_score', 
                                        pd.Series(0, index=adata.obs_names))
            
            # Get CXCL13 expression
            if 'CXCL13' in adata.var_names:
                cxcl13_idx = adata.var_names.get_loc('CXCL13')
                if hasattr(adata.X, 'toarray'):
                    cxcl13_expr = adata.X[:, cxcl13_idx].toarray().flatten()
                else:
                    cxcl13_expr = adata.X[:, cxcl13_idx].flatten()
            else:
                cxcl13_expr = np.zeros(len(adata))
            
            # Normalize CXCL13 expression
            cxcl13_normalized = (cxcl13_expr - np.min(cxcl13_expr)) / (np.max(cxcl13_expr) - np.min(cxcl13_expr) + 1e-8)
            
            # Calculate composite score
            cxcl13_t_composite = np.array(t_cell_score) * cxcl13_normalized
            
            # Define CXCL13+ T cell threshold
            t_high = np.array(t_cell_score) > np.percentile(t_cell_score, 75)
            cxcl13_high = cxcl13_expr > np.percentile(cxcl13_expr, 75)
            cxcl13_t_mask = t_high & cxcl13_high
            
            # Save to adata
            adata.obs['CXCL13_T_composite_score'] = cxcl13_t_composite
            adata.obs['CXCL13_T_mask'] = cxcl13_t_mask
            
            # Save to analyzer
            self.cxcl13_t_scores[sample_id] = {
                'composite_score': cxcl13_t_composite,
                'mask': cxcl13_t_mask,
                'n_cxcl13_t': cxcl13_t_mask.sum()
            }
            
            print(f"    IL1R1+ CAF: {il1r1_caf_mask.sum()} spots ({il1r1_caf_mask.sum()/len(adata)*100:.1f}%)")
            print(f"    CXCL13+ T: {cxcl13_t_mask.sum()} spots ({cxcl13_t_mask.sum()/len(adata)*100:.1f}%)")
    
    def analyze_colocalization(self):
        """
        Analyze spatial colocalization between IL1R1+ CAF and CXCL13+ T cells.
        
        Performs:
        1. Score correlation analysis
        2. Spatial proximity analysis with permutation tests
        3. Sample-wise correlation analysis
        
        Returns:
            dict: Comprehensive colocalization results
        """
        print("Analyzing IL1R1+ CAF and CXCL13+ T cell colocalization...")
        
        all_il1r1_scores = []
        all_cxcl13_scores = []
        sample_labels = []
        sample_correlations = {}
        spatial_proximity_results = {}
        
        for sample_id in self.spatial_data.keys():
            if sample_id in self.il1r1_caf_scores and sample_id in self.cxcl13_t_scores:
                print(f"\n  Analyzing sample: {sample_id}")
                
                adata = self.spatial_data[sample_id]
                
                # Get composite scores
                il1r1_scores = self.il1r1_caf_scores[sample_id]['composite_score']
                cxcl13_scores = self.cxcl13_t_scores[sample_id]['composite_score']
                
                # Calculate score correlation
                if len(il1r1_scores) > 10:
                    pearson_r, pearson_p = pearsonr(il1r1_scores, cxcl13_scores)
                    spearman_r, spearman_p = spearmanr(il1r1_scores, cxcl13_scores)
                    
                    sample_correlations[sample_id] = {
                        'pearson_r': pearson_r,
                        'pearson_p': pearson_p,
                        'spearman_r': spearman_r,
                        'spearman_p': spearman_p,
                        'n_spots': len(il1r1_scores),
                        'n_il1r1_caf': self.il1r1_caf_scores[sample_id]['n_il1r1_caf'],
                        'n_cxcl13_t': self.cxcl13_t_scores[sample_id]['n_cxcl13_t']
                    }
                    
                    print(f"    Score correlation: Pearson r = {pearson_r:.3f} (p = {pearson_p:.3f})")
                
                # Spatial proximity analysis
                proximity_result = self._analyze_spatial_proximity(sample_id)
                if proximity_result:
                    spatial_proximity_results[sample_id] = proximity_result
                
                # Collect for overall analysis
                all_il1r1_scores.extend(il1r1_scores)
                all_cxcl13_scores.extend(cxcl13_scores)
                sample_labels.extend([sample_id] * len(il1r1_scores))
        
        # Calculate overall correlation
        overall_correlation = {}
        if len(all_il1r1_scores) > 10:
            overall_pearson_r, overall_pearson_p = pearsonr(all_il1r1_scores, all_cxcl13_scores)
            overall_spearman_r, overall_spearman_p = spearmanr(all_il1r1_scores, all_cxcl13_scores)
            
            overall_correlation = {
                'pearson_r': overall_pearson_r,
                'pearson_p': overall_pearson_p,
                'spearman_r': overall_spearman_r,
                'spearman_p': overall_spearman_p,
                'n_total_spots': len(all_il1r1_scores),
                'n_samples': len(set(sample_labels))
            }
            
            print(f"\n  Overall correlation analysis:")
            print(f"    Pearson r = {overall_pearson_r:.3f} (p = {overall_pearson_p:.3f})")
            print(f"    Spearman ρ = {overall_spearman_r:.3f} (p = {overall_spearman_p:.3f})")
            print(f"    Total spots: {len(all_il1r1_scores)}, Samples: {len(set(sample_labels))}")
        
        # Store results
        self.colocalization_results = {
            'sample_correlations': sample_correlations,
            'spatial_proximity': spatial_proximity_results,
            'overall_correlation': overall_correlation,
            'sample_details': {
                'all_il1r1_caf_scores': all_il1r1_scores,
                'all_cxcl13_t_scores': all_cxcl13_scores,
                'sample_labels': sample_labels
            }
        }
        
        return self.colocalization_results
    
    def _analyze_spatial_proximity(self, sample_id, radius=50):
        """
        Analyze spatial proximity between IL1R1+ CAF and CXCL13+ T cells.
        
        Args:
            sample_id (str): Sample identifier
            radius (float): Distance threshold for neighbor detection
            
        Returns:
            dict: Spatial proximity analysis results
        """
        print(f"    Analyzing spatial proximity for sample {sample_id}...")
        
        adata = self.spatial_data[sample_id]
        coords = adata.obsm['spatial']
        
        # Get cell masks
        il1r1_caf_mask = adata.obs['IL1R1_CAF_mask'].values
        cxcl13_t_mask = adata.obs['CXCL13_T_mask'].values
        
        if il1r1_caf_mask.sum() == 0 or cxcl13_t_mask.sum() == 0:
            print(f"      Insufficient cells for proximity analysis")
            return None
        
        # Get coordinates
        il1r1_coords = coords[il1r1_caf_mask]
        cxcl13_coords = coords[cxcl13_t_mask]
        
        # Calculate distance matrix
        distances = pairwise_distances(il1r1_coords, cxcl13_coords)
        
        # Calculate proximity metrics
        min_distances = distances.min(axis=1)  # Min distance to nearest CXCL13+ T
        neighbors_within_radius = (distances <= radius).sum(axis=1)  # Neighbors within radius
        
        # Permutation test for random expectation
        n_permutations = 1000
        random_min_distances = []
        random_neighbors = []
        
        for _ in range(n_permutations):
            shuffled_coords = cxcl13_coords[np.random.permutation(len(cxcl13_coords))]
            random_distances = pairwise_distances(il1r1_coords, shuffled_coords)
            random_min_distances.extend(random_distances.min(axis=1))
            random_neighbors.extend((random_distances <= radius).sum(axis=1))
        
        # Statistical tests
        from scipy.stats import mannwhitneyu
        
        stat_distance, p_distance = mannwhitneyu(min_distances, random_min_distances, alternative='less')
        stat_neighbors, p_neighbors = mannwhitneyu(neighbors_within_radius, random_neighbors, alternative='greater')
        
        result = {
            'mean_min_distance': np.mean(min_distances),
            'random_mean_min_distance': np.mean(random_min_distances),
            'mean_neighbors': np.mean(neighbors_within_radius),
            'random_mean_neighbors': np.mean(random_neighbors),
            'p_distance': p_distance,
            'p_neighbors': p_neighbors,
            'n_il1r1_caf': len(il1r1_coords),
            'n_cxcl13_t': len(cxcl13_coords)
        }
        
        print(f"      Mean min distance: {np.mean(min_distances):.2f} vs random {np.mean(random_min_distances):.2f} (p = {p_distance:.3f})")
        print(f"      Mean neighbors: {np.mean(neighbors_within_radius):.2f} vs random {np.mean(random_neighbors):.2f} (p = {p_neighbors:.3f})")
        
        return result
    
    def plot_colocalization_results(self):
        """
        Generate comprehensive visualization of colocalization results.
        
        Creates a 2x3 subplot figure with:
        1. Overall score correlation scatter plot
        2. Per-sample correlation coefficients
        3. Spatial proximity p-values
        4. Cell count distribution
        5. Score distributions
        6. Summary statistics text
        """
        if not self.colocalization_results:
            print("No colocalization results to visualize")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        
        # 1. Overall correlation scatter plot
        ax1 = axes[0, 0]
        all_il1r1 = self.colocalization_results['sample_details']['all_il1r1_caf_scores']
        all_cxcl13 = self.colocalization_results['sample_details']['all_cxcl13_t_scores']
        
        # Subsample if too many points
        if len(all_il1r1) > 3000:
            indices = np.random.choice(len(all_il1r1), 3000, replace=False)
            plot_il1r1 = [all_il1r1[i] for i in indices]
            plot_cxcl13 = [all_cxcl13[i] for i in indices]
        else:
            plot_il1r1 = all_il1r1
            plot_cxcl13 = all_cxcl13
        
        ax1.scatter(plot_il1r1, plot_cxcl13, alpha=0.5, s=8, color='darkred')
        
        # Add regression line
        if len(all_il1r1) > 1:
            z = np.polyfit(all_il1r1, all_cxcl13, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(all_il1r1), max(all_il1r1), 100)
            ax1.plot(x_line, p(x_line), "orange", alpha=0.8, linewidth=3)
        
        # Add statistics
        overall = self.colocalization_results.get('overall_correlation', {})
        if overall:
            stats_text = f'Pearson r = {overall["pearson_r"]:.3f}\n'
            stats_text += f'p = {overall["pearson_p"]:.2e}\n'
            stats_text += f'n = {overall["n_total_spots"]} spots'
            ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                    verticalalignment='top', fontsize=14)
        
        ax1.set_xlabel('IL1R1+ CAF Composite Score', fontsize=20)
        ax1.set_ylabel('CXCL13+ T Cell Composite Score', fontsize=20)
        ax1.set_title('IL1R1+ CAF vs CXCL13+ T Cell\nCorrelation Analysis', fontsize=23)
        ax1.tick_params(labelsize=16)
        ax1.grid(True, alpha=0.3)
        
        # 2. Per-sample correlation coefficients
        ax2 = axes[0, 1]
        sample_corrs = self.colocalization_results.get('sample_correlations', {})
        if sample_corrs:
            samples = list(sample_corrs.keys())
            pearson_rs = [sample_corrs[s]['pearson_r'] for s in samples]
            p_values = [sample_corrs[s]['pearson_p'] for s in samples]
            
            colors = ['darkred' if p < 0.05 else 'gray' for p in p_values]
            bars = ax2.bar(range(len(samples)), pearson_rs, color=colors, alpha=0.7)
            
            ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            ax2.set_xlabel('Sample ID', fontsize=20)
            ax2.set_ylabel('Pearson r', fontsize=20)
            ax2.set_title('Per-Sample Correlations', fontsize=23)
            ax2.set_xticks(range(len(samples)))
            ax2.set_xticklabels(samples, rotation=45, ha='right', fontsize=12)
            ax2.tick_params(labelsize=16)
            
            # Add significance annotations
            for i, (bar, p_val) in enumerate(zip(bars, p_values)):
                height = bar.get_height()
                significance = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                y_pos = height + 0.02 if height > 0 else height - 0.05
                ax2.text(bar.get_x() + bar.get_width()/2., y_pos,
                        significance, ha='center', va='bottom' if height > 0 else 'top', fontsize=12)
        
        # 3. Spatial proximity analysis
        ax3 = axes[0, 2]
        spatial_results = self.colocalization_results.get('spatial_proximity', {})
        if spatial_results:
            samples = list(spatial_results.keys())
            p_distances = [spatial_results[s]['p_distance'] for s in samples]
            p_neighbors = [spatial_results[s]['p_neighbors'] for s in samples]
            
            x = np.arange(len(samples))
            width = 0.35
            
            bars1 = ax3.bar(x - width/2, [-np.log10(p) for p in p_distances], width,
                           label='Distance p-value', alpha=0.7, color='blue')
            bars2 = ax3.bar(x + width/2, [-np.log10(p) for p in p_neighbors], width,
                           label='Neighbors p-value', alpha=0.7, color='red')
            
            ax3.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7, label='p = 0.05')
            ax3.set_xlabel('Sample ID', fontsize=20)
            ax3.set_ylabel('-log10(p-value)', fontsize=20)
            ax3.set_title('Spatial Proximity Analysis', fontsize=23)
            ax3.set_xticks(x)
            ax3.set_xticklabels(samples, rotation=45, ha='right', fontsize=12)
            ax3.tick_params(labelsize=16)
            ax3.legend(fontsize=14)
        
        # 4. Cell count distribution
        ax4 = axes[1, 0]
        if sample_corrs:
            samples = list(sample_corrs.keys())
            il1r1_counts = [sample_corrs[s]['n_il1r1_caf'] for s in samples]
            cxcl13_counts = [sample_corrs[s]['n_cxcl13_t'] for s in samples]
            
            x = np.arange(len(samples))
            width = 0.35
            
            ax4.bar(x - width/2, il1r1_counts, width, label='IL1R1+ CAF', alpha=0.7, color='darkred')
            ax4.bar(x + width/2, cxcl13_counts, width, label='CXCL13+ T cells', alpha=0.7, color='darkblue')
            
            ax4.set_xlabel('Sample ID', fontsize=20)
            ax4.set_ylabel('Cell Count', fontsize=20)
            ax4.set_title('Cell Type Abundance', fontsize=23)
            ax4.set_xticks(x)
            ax4.set_xticklabels(samples, rotation=45, ha='right', fontsize=12)
            ax4.tick_params(labelsize=16)
            ax4.legend(fontsize=16)
        
        # 5. Score distributions
        ax5 = axes[1, 1]
        ax5.hist(all_il1r1, bins=30, alpha=0.7, color='darkred', label='IL1R1+ CAF', density=True)
        ax5.hist(all_cxcl13, bins=30, alpha=0.7, color='darkblue', label='CXCL13+ T cells', density=True)
        ax5.set_xlabel('Composite Score', fontsize=20)
        ax5.set_ylabel('Density', fontsize=20)
        ax5.set_title('Score Distributions', fontsize=23)
        ax5.tick_params(labelsize=16)
        ax5.legend(fontsize=18)
        ax5.grid(True, alpha=0.3)
        
        # 6. Summary text
        ax6 = axes[1, 2]
        ax6.axis('off')
        
        summary_text = f"IL1R1+ CAF vs CXCL13+ T Cell Analysis\n"
        summary_text += f"{'='*40}\n\n"
        
        if overall:
            summary_text += f"Overall Correlation:\n"
            summary_text += f"• Pearson r: {overall['pearson_r']:.3f}\n"
            summary_text += f"• p-value: {overall['pearson_p']:.2e}\n"
            summary_text += f"• Total spots: {overall['n_total_spots']:,}\n"
            summary_text += f"• Samples: {overall['n_samples']}\n\n"
            
            if overall['pearson_p'] < 0.05:
                if overall['pearson_r'] > 0:
                    summary_text += f"✓ Significant POSITIVE correlation\n"
                    summary_text += f"  → Spatial co-localization detected\n\n"
                else:
                    summary_text += f"✓ Significant NEGATIVE correlation\n"
                    summary_text += f"  → Spatial exclusion detected\n\n"
            else:
                summary_text += f"- No significant correlation\n\n"
        
        if sample_corrs:
            sig_samples = sum(1 for corr in sample_corrs.values() if corr['pearson_p'] < 0.05)
            summary_text += f"Sample-wise Analysis:\n"
            summary_text += f"• {sig_samples}/{len(sample_corrs)} samples significant\n"
        
        ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
                fontsize=14, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        plt.savefig('IL1R1_CAF_vs_CXCL13_T_colocalization_analysis.pdf',
                    dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig('IL1R1_CAF_vs_CXCL13_T_colocalization_analysis.png',
                    dpi=300, bbox_inches='tight')
        plt.show()
    
    def save_results(self):
        """Save all analysis results to CSV files."""
        print("Saving analysis results...")
        
        # Save colocalization results
        if self.colocalization_results:
            results_data = []
            
            # Overall results
            overall = self.colocalization_results.get('overall_correlation', {})
            if overall:
                results_data.append({
                    'analysis_type': 'IL1R1_CAF_vs_CXCL13_T_overall',
                    'pearson_r': overall['pearson_r'],
                    'pearson_p': overall['pearson_p'],
                    'spearman_r': overall['spearman_r'],
                    'spearman_p': overall['spearman_p'],
                    'n_spots': overall['n_total_spots'],
                    'n_samples': overall['n_samples']
                })
            
            # Per-sample results
            sample_corrs = self.colocalization_results.get('sample_correlations', {})
            for sample_id, corr_data in sample_corrs.items():
                results_data.append({
                    'analysis_type': f'IL1R1_CAF_vs_CXCL13_T_{sample_id}',
                    'sample_id': sample_id,
                    'pearson_r': corr_data['pearson_r'],
                    'pearson_p': corr_data['pearson_p'],
                    'spearman_r': corr_data['spearman_r'],
                    'spearman_p': corr_data['spearman_p'],
                    'n_spots': corr_data['n_spots'],
                    'n_il1r1_caf': corr_data['n_il1r1_caf'],
                    'n_cxcl13_t': corr_data['n_cxcl13_t']
                })
            
            pd.DataFrame(results_data).to_csv('IL1R1_CAF_CXCL13_T_colocalization_results.csv', index=False)
            print("  Saved: IL1R1_CAF_CXCL13_T_colocalization_results.csv")
        
        # Save spatial proximity results
        spatial_results = self.colocalization_results.get('spatial_proximity', {})
        if spatial_results:
            spatial_data = []
            for sample_id, prox_data in spatial_results.items():
                spatial_data.append({
                    'sample_id': sample_id,
                    **prox_data
                })
            pd.DataFrame(spatial_data).to_csv('IL1R1_CAF_CXCL13_T_spatial_proximity_results.csv', index=False)
            print("  Saved: IL1R1_CAF_CXCL13_T_spatial_proximity_results.csv")
    
    def run_extended_analysis(self):
        """
        Run the complete extended analysis pipeline.
        
        Returns:
            dict: Complete analysis results
        """
        print("=" * 70)
        print("IL1R1+ CAF vs CXCL13+ T Cell Extended Analysis")
        print("=" * 70)
        
        # 1. Load data
        print("\n1. Loading spatial data...")
        self.load_spatial_data()
        if not self.spatial_data:
            print("No data loaded. Exiting.")
            return None
        
        # 2. Define gene signatures
        print("\n2. Defining gene signatures...")
        self.define_gene_signatures()
        
        # 3. Annotate cell types
        print("\n3. Annotating cell types...")
        self.annotate_cell_types()
        
        # 4. Calculate signature scores
        print("\n4. Calculating signature scores...")
        self.calculate_signature_scores()
        
        # 5. Identify target cell populations
        print("\n5. Identifying IL1R1+ CAF and CXCL13+ T cells...")
        self.identify_il1r1_caf_and_cxcl13_t()
        
        # 6. Analyze colocalization
        print("\n6. Analyzing colocalization...")
        results = self.analyze_colocalization()
        
        # 7. Generate visualizations
        print("\n7. Generating visualizations...")
        self.plot_colocalization_results()
        
        # 8. Save results
        print("\n8. Saving results...")
        self.save_results()
        
        print("\nExtended analysis complete!")
        return results


if __name__ == "__main__":
    # Run extended analysis
    analyzer = IL1R1_CXCL13_ColocalizationAnalyzer()
    results = analyzer.run_extended_analysis()

