#!/usr/bin/env python3
"""
====================================================================
Breast Cancer Spatial Transcriptomics - Comprehensive Analysis Module
====================================================================
Comprehensive analysis including:
1. IL1R1+ CAF vs CXCL13+ T cell colocalization
2. CDKN1A+ Macrophage vs CAF colocalization
3. Tumor cell (epithelial) distribution
4. All samples spatial overlay visualization with H&E staining

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
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
import matplotlib.patches as patches
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
from sklearn.preprocessing import StandardScaler
from scipy import stats
from matplotlib.colors import ListedColormap

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


class ComprehensiveSpatialAnalyzer:
    """
    Comprehensive spatial transcriptomics analyzer for breast cancer.
    
    This class provides complete analysis including:
    - Cell type annotation using gene signatures
    - IL1R1+ CAF and CXCL13+ T cell identification
    - CDKN1A+ Macrophage identification
    - Tumor cell (epithelial) identification
    - Spatial colocalization analysis
    - Multi-sample visualization with H&E staining
    
    Attributes:
        data_dir (Path): Data directory path
        spatial_data (dict): Sample AnnData objects
        gene_signatures (dict): Cell type gene signatures
        il1r1_caf_scores (dict): IL1R1+ CAF scores per sample
        cxcl13_t_scores (dict): CXCL13+ T cell scores per sample
        cdkn1a_macrophage_scores (dict): CDKN1A+ Macrophage scores per sample
        tumor_cell_scores (dict): Tumor cell scores per sample
        colocalization_results (dict): IL1R1-CXCL13 colocalization results
        cdkn1a_caf_colocalization_results (dict): CDKN1A-CAF colocalization results
    """
    
    def __init__(self, data_dir="."):
        """Initialize the comprehensive analyzer."""
        self.data_dir = Path(data_dir)
        self.spatial_data = {}
        self.gene_signatures = {}
        self.combined_adata = None
        
        # Cell type scores
        self.il1r1_caf_scores = {}
        self.cxcl13_t_scores = {}
        self.cdkn1a_macrophage_scores = {}
        self.tumor_cell_scores = {}
        
        # Analysis results
        self.colocalization_results = {}
        self.cdkn1a_caf_colocalization_results = {}
        self.abundance_results = {}
    
    # =========================================================================
    # Data Loading and Preprocessing
    # =========================================================================
    
    def load_spatial_data(self):
        """Load all 10X Visium spatial transcriptomics samples."""
        print("Loading spatial transcriptomics data...")
        
        spatial_dirs = glob.glob(str(self.data_dir / "download_page_*_spatial"))
        
        for spatial_dir in spatial_dirs:
            sample_id = Path(spatial_dir).name.replace("download_page_", "").replace("_spatial", "")
            print(f"  Loading sample: {sample_id}")
            
            try:
                h5_files = glob.glob(f"{spatial_dir}/*.h5")
                if h5_files:
                    h5_file = os.path.basename(h5_files[0])
                    adata = sc.read_visium(
                        path=spatial_dir,
                        count_file=h5_file,
                        load_images=True
                    )
                    
                    adata.var_names_make_unique()
                    adata.var['mt'] = adata.var_names.str.startswith('MT-')
                    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
                    adata.obs['sample_id'] = sample_id
                    
                    # Normalize
                    sc.pp.normalize_total(adata, target_sum=1e4)
                    sc.pp.log1p(adata)
                    
                    self.spatial_data[sample_id] = adata
                    print(f"    Loaded: {adata.shape}")
                    
            except Exception as e:
                print(f"    Failed: {e}")
        
        print(f"Total samples loaded: {len(self.spatial_data)}")
        return self.spatial_data
    
    def define_gene_signatures(self):
        """Define gene signatures for all cell types of interest."""
        print("Defining gene signatures...")
        
        self.gene_signatures = {
            # CAF signature
            'CAF': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'POSTN', 
                   'FAP', 'PDGFRA', 'PDGFRB', 'ACTA2', 'FN1', 'VIM', 'THY1'],
            
            # Fibroblast general
            'Fibroblasts': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA', 'FAP', 'THY1', 'VIM'],
            
            # T cell signature
            'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 
                       'TRAC', 'TRBC1', 'TRBC2', 'LCK', 'ZAP70'],
            
            # CXCL13+ T cell (TFH-like/exhausted)
            'CXCL13_T': ['CXCL13', 'PDCD1', 'CTLA4', 'TIGIT', 'LAG3', 
                        'HAVCR2', 'TOX', 'ICOS'],
            
            # Macrophage signature
            'Macrophages': ['CD68', 'CD163', 'CD14', 'MARCO', 'MSR1', 
                           'MRC1', 'FCGR3A', 'CSF1R'],
            
            # Epithelial/Tumor cell signature
            'Epithelial': ['EPCAM', 'KRT7', 'KRT8', 'KRT18', 'KRT19', 
                          'CDH1', 'MUC1'],
            
            # B cell signature
            'B_cells': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5'],
            
            # Endothelial signature
            'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'ERG']
        }
        
        print(f"  Defined {len(self.gene_signatures)} gene signatures")
    
    def calculate_signature_scores(self):
        """Calculate gene signature scores for all samples."""
        print("Calculating gene signature scores...")
        
        if not self.gene_signatures:
            self.define_gene_signatures()
        
        for sample_id, adata in self.spatial_data.items():
            print(f"  Processing sample: {sample_id}")
            
            for sig_name, genes in self.gene_signatures.items():
                available_genes = [g for g in genes if g in adata.var_names]
                
                if len(available_genes) >= 2:
                    try:
                        sc.tl.score_genes(adata, available_genes,
                                         score_name=f'{sig_name}_signature_score',
                                         use_raw=False)
                    except:
                        adata.obs[f'{sig_name}_signature_score'] = 0
                elif len(available_genes) == 1:
                    gene_idx = adata.var_names.get_loc(available_genes[0])
                    if hasattr(adata.X, 'toarray'):
                        adata.obs[f'{sig_name}_signature_score'] = adata.X[:, gene_idx].toarray().flatten()
                    else:
                        adata.obs[f'{sig_name}_signature_score'] = adata.X[:, gene_idx].flatten()
                else:
                    adata.obs[f'{sig_name}_signature_score'] = 0
    
    # =========================================================================
    # Cell Population Identification
    # =========================================================================
    
    def identify_il1r1_caf(self):
        """Identify IL1R1+ CAF cells."""
        print("Identifying IL1R1+ CAF cells...")
        
        for sample_id, adata in self.spatial_data.items():
            # Get CAF signature
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
            
            # Calculate composite score
            il1r1_normalized = (il1r1_expr - np.min(il1r1_expr)) / (np.max(il1r1_expr) - np.min(il1r1_expr) + 1e-8)
            il1r1_caf_composite = np.array(caf_score) * il1r1_normalized
            
            # Define mask (both high)
            caf_high = np.array(caf_score) > np.percentile(caf_score, 75)
            il1r1_high = il1r1_expr > np.percentile(il1r1_expr, 75)
            il1r1_caf_mask = caf_high & il1r1_high
            
            # Save
            adata.obs['IL1R1_CAF_composite_score'] = il1r1_caf_composite
            adata.obs['IL1R1_CAF_mask'] = il1r1_caf_mask
            
            self.il1r1_caf_scores[sample_id] = {
                'composite_score': il1r1_caf_composite,
                'mask': il1r1_caf_mask,
                'n_il1r1_caf': il1r1_caf_mask.sum()
            }
            
            print(f"    Sample {sample_id}: {il1r1_caf_mask.sum()} IL1R1+ CAF spots")
    
    def identify_cxcl13_t_cells(self):
        """Identify CXCL13+ T cells."""
        print("Identifying CXCL13+ T cells...")
        
        for sample_id, adata in self.spatial_data.items():
            # Get T cell signature
            t_score = adata.obs.get('T_cells_signature_score',
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
            
            # Calculate composite score
            cxcl13_normalized = (cxcl13_expr - np.min(cxcl13_expr)) / (np.max(cxcl13_expr) - np.min(cxcl13_expr) + 1e-8)
            cxcl13_t_composite = np.array(t_score) * cxcl13_normalized
            
            # Define mask
            t_high = np.array(t_score) > np.percentile(t_score, 75)
            cxcl13_high = cxcl13_expr > np.percentile(cxcl13_expr, 75)
            cxcl13_t_mask = t_high & cxcl13_high
            
            # Save
            adata.obs['CXCL13_T_composite_score'] = cxcl13_t_composite
            adata.obs['CXCL13_T_mask'] = cxcl13_t_mask
            
            self.cxcl13_t_scores[sample_id] = {
                'composite_score': cxcl13_t_composite,
                'mask': cxcl13_t_mask,
                'n_cxcl13_t': cxcl13_t_mask.sum()
            }
            
            print(f"    Sample {sample_id}: {cxcl13_t_mask.sum()} CXCL13+ T cell spots")
    
    def identify_cdkn1a_macrophages(self):
        """Identify CDKN1A+ Macrophages."""
        print("Identifying CDKN1A+ Macrophages...")
        
        for sample_id, adata in self.spatial_data.items():
            # Get Macrophage signature
            macro_score = adata.obs.get('Macrophages_signature_score',
                                       pd.Series(0, index=adata.obs_names))
            
            # Get CDKN1A expression
            if 'CDKN1A' in adata.var_names:
                cdkn1a_idx = adata.var_names.get_loc('CDKN1A')
                if hasattr(adata.X, 'toarray'):
                    cdkn1a_expr = adata.X[:, cdkn1a_idx].toarray().flatten()
                else:
                    cdkn1a_expr = adata.X[:, cdkn1a_idx].flatten()
            else:
                print(f"    Sample {sample_id}: CDKN1A not found")
                cdkn1a_expr = np.zeros(len(adata))
            
            # Calculate composite score
            cdkn1a_normalized = (cdkn1a_expr - np.min(cdkn1a_expr)) / (np.max(cdkn1a_expr) - np.min(cdkn1a_expr) + 1e-8)
            cdkn1a_macro_composite = np.array(macro_score) * cdkn1a_normalized
            
            # Define mask
            macro_high = np.array(macro_score) > np.percentile(macro_score, 75)
            cdkn1a_high = cdkn1a_expr > np.percentile(cdkn1a_expr, 75)
            cdkn1a_macro_mask = macro_high & cdkn1a_high
            
            # Save
            adata.obs['CDKN1A_Macrophage_composite_score'] = cdkn1a_macro_composite
            adata.obs['CDKN1A_Macrophage_mask'] = cdkn1a_macro_mask
            adata.obs['CDKN1A_expression'] = cdkn1a_expr
            
            self.cdkn1a_macrophage_scores[sample_id] = {
                'composite_score': cdkn1a_macro_composite,
                'mask': cdkn1a_macro_mask,
                'macro_scores': np.array(macro_score),
                'cdkn1a_expr': cdkn1a_expr,
                'n_cdkn1a_macro': cdkn1a_macro_mask.sum()
            }
            
            print(f"    Sample {sample_id}: {cdkn1a_macro_mask.sum()} CDKN1A+ Macrophage spots")
    
    def identify_tumor_cells(self):
        """Identify tumor cells (epithelial cells)."""
        print("Identifying tumor cells (epithelial)...")
        
        for sample_id, adata in self.spatial_data.items():
            # Get Epithelial signature
            epithelial_score = adata.obs.get('Epithelial_signature_score',
                                            pd.Series(0, index=adata.obs_names))
            
            # Define tumor cells as high epithelial signature
            tumor_threshold = np.percentile(epithelial_score, 75)
            tumor_mask = np.array(epithelial_score) > tumor_threshold
            
            # Save
            adata.obs['Tumor_cell_score'] = epithelial_score
            adata.obs['Tumor_cell_mask'] = tumor_mask
            
            self.tumor_cell_scores[sample_id] = {
                'scores': np.array(epithelial_score),
                'mask': tumor_mask,
                'n_tumor': tumor_mask.sum()
            }
            
            print(f"    Sample {sample_id}: {tumor_mask.sum()} tumor cell spots")
    
    # =========================================================================
    # Colocalization Analysis
    # =========================================================================
    
    def analyze_il1r1_cxcl13_colocalization(self):
        """Analyze IL1R1+ CAF vs CXCL13+ T cell colocalization."""
        print("Analyzing IL1R1+ CAF vs CXCL13+ T cell colocalization...")
        
        all_il1r1_scores = []
        all_cxcl13_scores = []
        sample_labels = []
        sample_correlations = {}
        
        for sample_id in self.spatial_data.keys():
            if sample_id in self.il1r1_caf_scores and sample_id in self.cxcl13_t_scores:
                adata = self.spatial_data[sample_id]
                
                il1r1_scores = self.il1r1_caf_scores[sample_id]['composite_score']
                cxcl13_scores = self.cxcl13_t_scores[sample_id]['composite_score']
                
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
                
                all_il1r1_scores.extend(il1r1_scores)
                all_cxcl13_scores.extend(cxcl13_scores)
                sample_labels.extend([sample_id] * len(il1r1_scores))
        
        # Overall correlation
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
            
            print(f"  Overall: r = {overall_pearson_r:.3f}, p = {overall_pearson_p:.2e}")
        
        self.colocalization_results = {
            'sample_correlations': sample_correlations,
            'overall_correlation': overall_correlation,
            'sample_details': {
                'all_il1r1_caf_scores': all_il1r1_scores,
                'all_cxcl13_t_scores': all_cxcl13_scores,
                'sample_labels': sample_labels
            }
        }
        
        return self.colocalization_results
    
    def analyze_cdkn1a_caf_colocalization(self):
        """Analyze CDKN1A+ Macrophage vs CAF colocalization."""
        print("Analyzing CDKN1A+ Macrophage vs CAF colocalization...")
        
        all_cdkn1a_scores = []
        all_caf_scores = []
        sample_labels = []
        sample_correlations = {}
        
        for sample_id in self.spatial_data.keys():
            if sample_id in self.cdkn1a_macrophage_scores:
                adata = self.spatial_data[sample_id]
                
                cdkn1a_scores = self.cdkn1a_macrophage_scores[sample_id]['composite_score']
                
                if 'Fibroblasts_signature_score' in adata.obs.columns:
                    caf_scores = adata.obs['Fibroblasts_signature_score'].values
                else:
                    continue
                
                if len(cdkn1a_scores) > 10:
                    pearson_r, pearson_p = pearsonr(cdkn1a_scores, caf_scores)
                    spearman_r, spearman_p = spearmanr(cdkn1a_scores, caf_scores)
                    
                    caf_mask = caf_scores > np.percentile(caf_scores, 75)
                    
                    sample_correlations[sample_id] = {
                        'pearson_r': pearson_r,
                        'pearson_p': pearson_p,
                        'spearman_r': spearman_r,
                        'spearman_p': spearman_p,
                        'n_spots': len(cdkn1a_scores),
                        'n_cdkn1a_macro': self.cdkn1a_macrophage_scores[sample_id]['n_cdkn1a_macro'],
                        'n_caf': caf_mask.sum()
                    }
                
                all_cdkn1a_scores.extend(cdkn1a_scores)
                all_caf_scores.extend(caf_scores)
                sample_labels.extend([sample_id] * len(cdkn1a_scores))
        
        # Overall correlation
        overall_correlation = {}
        if len(all_cdkn1a_scores) > 10:
            overall_pearson_r, overall_pearson_p = pearsonr(all_cdkn1a_scores, all_caf_scores)
            overall_spearman_r, overall_spearman_p = spearmanr(all_cdkn1a_scores, all_caf_scores)
            
            overall_correlation = {
                'pearson_r': overall_pearson_r,
                'pearson_p': overall_pearson_p,
                'spearman_r': overall_spearman_r,
                'spearman_p': overall_spearman_p,
                'n_total_spots': len(all_cdkn1a_scores),
                'n_samples': len(set(sample_labels))
            }
            
            print(f"  Overall: r = {overall_pearson_r:.3f}, p = {overall_pearson_p:.2e}")
        
        self.cdkn1a_caf_colocalization_results = {
            'sample_correlations': sample_correlations,
            'overall_correlation': overall_correlation,
            'sample_details': {
                'all_cdkn1a_macro_scores': all_cdkn1a_scores,
                'all_caf_scores': all_caf_scores,
                'sample_labels': sample_labels
            }
        }
        
        return self.cdkn1a_caf_colocalization_results
    
    # =========================================================================
    # Visualization
    # =========================================================================
    
    def plot_all_samples_il1r1_cxcl13_with_he(self):
        """
        Plot IL1R1+ CAF vs CXCL13+ T cells for all samples with H&E staining.
        
        Creates a grid visualization with H&E image and spatial overlay for each sample.
        """
        print("Plotting IL1R1+ CAF vs CXCL13+ T cells with H&E staining...")
        
        valid_samples = [s for s in self.spatial_data.keys()
                        if s in self.il1r1_caf_scores and s in self.cxcl13_t_scores]
        
        n_samples = len(valid_samples)
        if n_samples == 0:
            print("  No valid samples for plotting")
            return
        
        n_cols = 6  # 3 samples per row, 2 columns each (H&E + overlay)
        n_rows = (n_samples + 2) // 3
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
        
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        
        for idx, sample_id in enumerate(valid_samples):
            row = idx // 3
            col_start = (idx % 3) * 2
            
            # H&E staining
            ax_he = axes[row, col_start]
            self._plot_he_staining(sample_id, ax_he)
            
            # Spatial overlay
            ax_spatial = axes[row, col_start + 1]
            self._plot_il1r1_cxcl13_overlay(sample_id, ax_spatial)
        
        # Hide unused axes
        for idx in range(n_samples * 2, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].axis('off')
        
        # Unified legend
        self._add_il1r1_cxcl13_legend(fig)
        
        plt.tight_layout()
        plt.savefig('enhanced_all_samples_with_tumor_and_HE.pdf',
                    dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig('enhanced_all_samples_with_tumor_and_HE.png',
                    dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"  Saved visualization for {n_samples} samples")
    
    def plot_all_samples_cdkn1a_caf_with_he(self):
        """
        Plot CDKN1A+ Macrophages vs CAF for all samples with H&E staining.
        """
        print("Plotting CDKN1A+ Macrophages vs CAF with H&E staining...")
        
        valid_samples = [s for s in self.spatial_data.keys()
                        if s in self.cdkn1a_macrophage_scores
                        and 'Fibroblasts_signature_score' in self.spatial_data[s].obs.columns]
        
        n_samples = len(valid_samples)
        if n_samples == 0:
            print("  No valid samples for plotting")
            return
        
        n_cols = 6
        n_rows = (n_samples + 2) // 3
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
        
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        
        for idx, sample_id in enumerate(valid_samples):
            row = idx // 3
            col_start = (idx % 3) * 2
            
            ax_he = axes[row, col_start]
            self._plot_he_staining(sample_id, ax_he)
            
            ax_spatial = axes[row, col_start + 1]
            self._plot_cdkn1a_caf_overlay(sample_id, ax_spatial)
        
        # Hide unused
        for idx in range(n_samples * 2, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].axis('off')
        
        self._add_cdkn1a_caf_legend(fig)
        
        plt.tight_layout()
        plt.savefig('enhanced_all_samples_CDKN1A_Macrophage_CAF_with_HE.pdf',
                    dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig('enhanced_all_samples_CDKN1A_Macrophage_CAF_with_HE.png',
                    dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"  Saved visualization for {n_samples} samples")
    
    def _plot_he_staining(self, sample_id, ax):
        """Plot H&E staining image for a sample."""
        if sample_id not in self.spatial_data:
            ax.text(0.5, 0.5, f'Sample {sample_id}\nH&E Image\nNot Available',
                   ha='center', va='center', transform=ax.transAxes, fontsize=14)
            ax.axis('off')
            return
        
        adata = self.spatial_data[sample_id]
        
        try:
            if hasattr(adata, 'uns') and 'spatial' in adata.uns:
                spatial_key = list(adata.uns['spatial'].keys())[0]
                if 'images' in adata.uns['spatial'][spatial_key]:
                    if 'fullres' in adata.uns['spatial'][spatial_key]['images']:
                        img = adata.uns['spatial'][spatial_key]['images']['fullres']
                    elif 'hires' in adata.uns['spatial'][spatial_key]['images']:
                        img = adata.uns['spatial'][spatial_key]['images']['hires']
                    else:
                        img = None
                    
                    if img is not None:
                        ax.imshow(img)
                        ax.set_title(f'Sample {sample_id}\nH&E Staining', fontsize=16)
                        ax.axis('off')
                        return
        except Exception as e:
            print(f"    Warning loading H&E for {sample_id}: {e}")
        
        ax.text(0.5, 0.5, f'Sample {sample_id}\nH&E Staining\n(Image not available)',
               ha='center', va='center', transform=ax.transAxes, fontsize=14,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    
    def _plot_il1r1_cxcl13_overlay(self, sample_id, ax):
        """Plot IL1R1+ CAF and CXCL13+ T cell spatial overlay."""
        adata = self.spatial_data[sample_id]
        coords = adata.obsm['spatial']
        
        il1r1_mask = adata.obs['IL1R1_CAF_mask'].values
        cxcl13_mask = adata.obs['CXCL13_T_mask'].values
        tumor_mask = adata.obs.get('Tumor_cell_mask', np.zeros(len(adata), dtype=bool))
        
        # Background
        ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=6, alpha=0.3)
        
        # Tumor cells (green)
        if np.sum(tumor_mask) > 0:
            ax.scatter(coords[tumor_mask, 0], coords[tumor_mask, 1],
                      c='green', s=15, alpha=0.6, edgecolors='darkgreen', linewidth=0.2)
        
        # IL1R1+ CAF (red)
        if np.sum(il1r1_mask) > 0:
            ax.scatter(coords[il1r1_mask, 0], coords[il1r1_mask, 1],
                      c='red', s=20, alpha=0.8, edgecolors='darkred', linewidth=0.3)
        
        # CXCL13+ T (blue)
        if np.sum(cxcl13_mask) > 0:
            ax.scatter(coords[cxcl13_mask, 0], coords[cxcl13_mask, 1],
                      c='blue', s=20, alpha=0.8, edgecolors='darkblue', linewidth=0.3)
        
        # Overlap (purple)
        overlap = il1r1_mask & cxcl13_mask
        if np.sum(overlap) > 0:
            ax.scatter(coords[overlap, 0], coords[overlap, 1],
                      c='purple', s=25, alpha=0.9, edgecolors='black', linewidth=0.5)
        
        # Statistics text
        stats_text = f"Tumor: {np.sum(tumor_mask)}\n"
        stats_text += f"IL1R1+CAF: {np.sum(il1r1_mask)}\n"
        stats_text += f"CXCL13+T: {np.sum(cxcl13_mask)}\n"
        stats_text += f"Overlap: {np.sum(overlap)}"
        
        if sample_id in self.colocalization_results.get('sample_correlations', {}):
            corr = self.colocalization_results['sample_correlations'][sample_id]
            stats_text += f"\nr={corr['pearson_r']:.3f}"
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
               fontsize=9, verticalalignment='top', fontfamily='monospace')
        
        ax.set_title(f'Sample {sample_id}\nSpatial Cell Distribution', fontsize=16)
        ax.set_xlabel('Spatial X', fontsize=14)
        ax.set_ylabel('Spatial Y', fontsize=14)
    
    def _plot_cdkn1a_caf_overlay(self, sample_id, ax):
        """Plot CDKN1A+ Macrophage and CAF spatial overlay."""
        adata = self.spatial_data[sample_id]
        coords = adata.obsm['spatial']
        
        cdkn1a_macro_mask = adata.obs['CDKN1A_Macrophage_mask'].values
        caf_scores = adata.obs['Fibroblasts_signature_score'].values
        caf_mask = caf_scores > np.percentile(caf_scores, 75)
        tumor_mask = adata.obs.get('Tumor_cell_mask', np.zeros(len(adata), dtype=bool))
        
        # Background
        ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=6, alpha=0.3)
        
        # Tumor (green)
        if np.sum(tumor_mask) > 0:
            ax.scatter(coords[tumor_mask, 0], coords[tumor_mask, 1],
                      c='green', s=15, alpha=0.6, edgecolors='darkgreen', linewidth=0.2)
        
        # CAF (brown)
        if np.sum(caf_mask) > 0:
            ax.scatter(coords[caf_mask, 0], coords[caf_mask, 1],
                      c='brown', s=18, alpha=0.7, edgecolors='saddlebrown', linewidth=0.3)
        
        # CDKN1A+ Macrophage (orange)
        if np.sum(cdkn1a_macro_mask) > 0:
            ax.scatter(coords[cdkn1a_macro_mask, 0], coords[cdkn1a_macro_mask, 1],
                      c='orange', s=20, alpha=0.8, edgecolors='darkorange', linewidth=0.3)
        
        # Overlap (purple)
        overlap = cdkn1a_macro_mask & caf_mask
        if np.sum(overlap) > 0:
            ax.scatter(coords[overlap, 0], coords[overlap, 1],
                      c='purple', s=25, alpha=0.9, edgecolors='black', linewidth=0.5)
        
        # Statistics
        stats_text = f"Tumor: {np.sum(tumor_mask)}\n"
        stats_text += f"CAF: {np.sum(caf_mask)}\n"
        stats_text += f"CDKN1A+Macro: {np.sum(cdkn1a_macro_mask)}\n"
        stats_text += f"Overlap: {np.sum(overlap)}"
        
        if sample_id in self.cdkn1a_caf_colocalization_results.get('sample_correlations', {}):
            corr = self.cdkn1a_caf_colocalization_results['sample_correlations'][sample_id]
            stats_text += f"\nr={corr['pearson_r']:.3f}"
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
               fontsize=9, verticalalignment='top', fontfamily='monospace')
        
        ax.set_title(f'Sample {sample_id}\nCDKN1A+ Macrophage & CAF', fontsize=16)
        ax.set_xlabel('Spatial X', fontsize=14)
        ax.set_ylabel('Spatial Y', fontsize=14)
    
    def _add_il1r1_cxcl13_legend(self, fig):
        """Add unified legend for IL1R1-CXCL13 plot."""
        legend_elements = [
            plt.scatter([], [], c='lightgray', s=50, alpha=0.6, label='Other spots'),
            plt.scatter([], [], c='green', s=50, alpha=0.8, label='Tumor cells'),
            plt.scatter([], [], c='red', s=50, alpha=0.8, label='IL1R1+ CAF'),
            plt.scatter([], [], c='blue', s=50, alpha=0.8, label='CXCL13+ T cells'),
            plt.scatter([], [], c='purple', s=60, alpha=0.9, label='IL1R1+CAF & CXCL13+T overlap')
        ]
        fig.legend(handles=legend_elements, loc='lower center',
                  bbox_to_anchor=(0.5, -0.02), ncol=5, fontsize=18,
                  frameon=True, fancybox=True, shadow=True)
    
    def _add_cdkn1a_caf_legend(self, fig):
        """Add unified legend for CDKN1A-CAF plot."""
        legend_elements = [
            plt.scatter([], [], c='lightgray', s=50, alpha=0.6, label='Other spots'),
            plt.scatter([], [], c='green', s=50, alpha=0.8, label='Tumor cells'),
            plt.scatter([], [], c='brown', s=50, alpha=0.8, label='CAF'),
            plt.scatter([], [], c='orange', s=50, alpha=0.8, label='CDKN1A+ Macrophages'),
            plt.scatter([], [], c='purple', s=60, alpha=0.9, label='CDKN1A+Macro & CAF overlap')
        ]
        fig.legend(handles=legend_elements, loc='lower center',
                  bbox_to_anchor=(0.5, -0.02), ncol=5, fontsize=18,
                  frameon=True, fancybox=True, shadow=True)
    
    def plot_comparative_analysis(self):
        """
        Plot comparative analysis between IL1R1-CXCL13 and CDKN1A-CAF.
        """
        print("Plotting comparative analysis...")
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. IL1R1-CXCL13 overall correlation
        ax1 = axes[0, 0]
        if self.colocalization_results and self.colocalization_results.get('sample_details'):
            il1r1_scores = self.colocalization_results['sample_details']['all_il1r1_caf_scores']
            cxcl13_scores = self.colocalization_results['sample_details']['all_cxcl13_t_scores']
            
            if len(il1r1_scores) > 3000:
                idx = np.random.choice(len(il1r1_scores), 3000, replace=False)
                plot_il1r1 = [il1r1_scores[i] for i in idx]
                plot_cxcl13 = [cxcl13_scores[i] for i in idx]
            else:
                plot_il1r1 = il1r1_scores
                plot_cxcl13 = cxcl13_scores
            
            ax1.scatter(plot_il1r1, plot_cxcl13, alpha=0.5, s=6, color='darkred')
            
            if len(il1r1_scores) > 1:
                z = np.polyfit(il1r1_scores, cxcl13_scores, 1)
                p = np.poly1d(z)
                x_line = np.linspace(min(il1r1_scores), max(il1r1_scores), 100)
                ax1.plot(x_line, p(x_line), "orange", alpha=0.8, linewidth=2)
            
            overall = self.colocalization_results.get('overall_correlation', {})
            if overall:
                stats_text = f'r = {overall["pearson_r"]:.3f}\np = {overall["pearson_p"]:.2e}'
                ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                        verticalalignment='top', fontsize=12)
        
        ax1.set_xlabel('IL1R1+ CAF Score', fontsize=18)
        ax1.set_ylabel('CXCL13+ T Cell Score', fontsize=18)
        ax1.set_title('IL1R1+ CAF vs CXCL13+ T Cell', fontsize=20)
        ax1.grid(True, alpha=0.3)
        
        # 2. CDKN1A-CAF overall correlation
        ax2 = axes[0, 1]
        if self.cdkn1a_caf_colocalization_results and self.cdkn1a_caf_colocalization_results.get('sample_details'):
            cdkn1a_scores = self.cdkn1a_caf_colocalization_results['sample_details']['all_cdkn1a_macro_scores']
            caf_scores = self.cdkn1a_caf_colocalization_results['sample_details']['all_caf_scores']
            
            if len(cdkn1a_scores) > 3000:
                idx = np.random.choice(len(cdkn1a_scores), 3000, replace=False)
                plot_cdkn1a = [cdkn1a_scores[i] for i in idx]
                plot_caf = [caf_scores[i] for i in idx]
            else:
                plot_cdkn1a = cdkn1a_scores
                plot_caf = caf_scores
            
            ax2.scatter(plot_cdkn1a, plot_caf, alpha=0.5, s=6, color='darkgreen')
            
            if len(cdkn1a_scores) > 1:
                z = np.polyfit(cdkn1a_scores, caf_scores, 1)
                p = np.poly1d(z)
                x_line = np.linspace(min(cdkn1a_scores), max(cdkn1a_scores), 100)
                ax2.plot(x_line, p(x_line), "orange", alpha=0.8, linewidth=2)
            
            overall = self.cdkn1a_caf_colocalization_results.get('overall_correlation', {})
            if overall:
                stats_text = f'r = {overall["pearson_r"]:.3f}\np = {overall["pearson_p"]:.2e}'
                ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                        verticalalignment='top', fontsize=12)
        
        ax2.set_xlabel('CDKN1A+ Macrophage Score', fontsize=18)
        ax2.set_ylabel('CAF Score', fontsize=18)
        ax2.set_title('CDKN1A+ Macrophage vs CAF', fontsize=20)
        ax2.grid(True, alpha=0.3)
        
        # 3. Per-sample correlation comparison
        ax3 = axes[1, 0]
        il1r1_corrs = self.colocalization_results.get('sample_correlations', {})
        cdkn1a_corrs = self.cdkn1a_caf_colocalization_results.get('sample_correlations', {})
        
        common_samples = set(il1r1_corrs.keys()) & set(cdkn1a_corrs.keys())
        if common_samples:
            samples = sorted(list(common_samples))
            il1r1_rs = [il1r1_corrs[s]['pearson_r'] for s in samples]
            cdkn1a_rs = [cdkn1a_corrs[s]['pearson_r'] for s in samples]
            
            x = np.arange(len(samples))
            width = 0.35
            
            ax3.bar(x - width/2, il1r1_rs, width, label='IL1R1+CAF vs CXCL13+T',
                   alpha=0.7, color='darkred')
            ax3.bar(x + width/2, cdkn1a_rs, width, label='CDKN1A+Macro vs CAF',
                   alpha=0.7, color='darkgreen')
            
            ax3.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            ax3.set_xlabel('Sample ID', fontsize=18)
            ax3.set_ylabel('Pearson r', fontsize=18)
            ax3.set_title('Per-Sample Correlation Comparison', fontsize=20)
            ax3.set_xticks(x)
            ax3.set_xticklabels(samples, rotation=45, ha='right', fontsize=12)
            ax3.legend(fontsize=14)
        
        # 4. Summary
        ax4 = axes[1, 1]
        ax4.axis('off')
        
        summary_text = "Comparative Spatial Analysis Summary\n"
        summary_text += "=" * 45 + "\n\n"
        
        if self.colocalization_results and self.colocalization_results.get('overall_correlation'):
            overall = self.colocalization_results['overall_correlation']
            summary_text += "IL1R1+CAF vs CXCL13+T Cells:\n"
            summary_text += f"• Pearson r: {overall['pearson_r']:.3f}\n"
            summary_text += f"• p-value: {overall['pearson_p']:.2e}\n"
            summary_text += f"• Significant: {'Yes' if overall['pearson_p'] < 0.05 else 'No'}\n\n"
        
        if self.cdkn1a_caf_colocalization_results and self.cdkn1a_caf_colocalization_results.get('overall_correlation'):
            overall = self.cdkn1a_caf_colocalization_results['overall_correlation']
            summary_text += "CDKN1A+Macrophage vs CAF:\n"
            summary_text += f"• Pearson r: {overall['pearson_r']:.3f}\n"
            summary_text += f"• p-value: {overall['pearson_p']:.2e}\n"
            summary_text += f"• Significant: {'Yes' if overall['pearson_p'] < 0.05 else 'No'}\n"
        
        ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
                fontsize=12, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        plt.savefig('comparative_spatial_analysis_summary.pdf',
                    dpi=300, bbox_inches='tight', format='pdf')
        plt.savefig('comparative_spatial_analysis_summary.png',
                    dpi=300, bbox_inches='tight')
        plt.show()
    
    # =========================================================================
    # Results Saving
    # =========================================================================
    
    def save_all_results(self):
        """Save all analysis results to CSV files."""
        print("Saving all results...")
        
        # IL1R1-CXCL13 colocalization results
        if self.colocalization_results:
            results_data = []
            
            overall = self.colocalization_results.get('overall_correlation', {})
            if overall:
                results_data.append({
                    'analysis_type': 'IL1R1_CAF_vs_CXCL13_T_overall',
                    **overall
                })
            
            for sample_id, corr in self.colocalization_results.get('sample_correlations', {}).items():
                results_data.append({
                    'analysis_type': f'IL1R1_CAF_vs_CXCL13_T_{sample_id}',
                    'sample_id': sample_id,
                    **corr
                })
            
            pd.DataFrame(results_data).to_csv('IL1R1_CAF_CXCL13_T_colocalization_results.csv', index=False)
        
        # CDKN1A-CAF colocalization results
        if self.cdkn1a_caf_colocalization_results:
            results_data = []
            
            overall = self.cdkn1a_caf_colocalization_results.get('overall_correlation', {})
            if overall:
                results_data.append({
                    'analysis_type': 'CDKN1A_Macro_vs_CAF_overall',
                    **overall
                })
            
            for sample_id, corr in self.cdkn1a_caf_colocalization_results.get('sample_correlations', {}).items():
                results_data.append({
                    'analysis_type': f'CDKN1A_Macro_vs_CAF_{sample_id}',
                    'sample_id': sample_id,
                    **corr
                })
            
            pd.DataFrame(results_data).to_csv('CDKN1A_Macrophage_CAF_colocalization_results.csv', index=False)
        
        print("  Results saved successfully")
    
    # =========================================================================
    # Main Analysis Pipeline
    # =========================================================================
    
    def run_comprehensive_analysis(self):
        """
        Run the complete comprehensive analysis pipeline.
        
        Returns:
            dict: All analysis results
        """
        print("=" * 70)
        print("BRCA Spatial Transcriptomics - Comprehensive Analysis")
        print("=" * 70)
        
        # 1. Load data
        print("\n1. Loading spatial data...")
        self.load_spatial_data()
        if not self.spatial_data:
            print("No data loaded. Exiting.")
            return None
        
        # 2. Define signatures
        print("\n2. Defining gene signatures...")
        self.define_gene_signatures()
        
        # 3. Calculate signature scores
        print("\n3. Calculating signature scores...")
        self.calculate_signature_scores()
        
        # 4. Identify cell populations
        print("\n4. Identifying cell populations...")
        self.identify_il1r1_caf()
        self.identify_cxcl13_t_cells()
        self.identify_cdkn1a_macrophages()
        self.identify_tumor_cells()
        
        # 5. Colocalization analysis
        print("\n5. Analyzing colocalization...")
        il1r1_results = self.analyze_il1r1_cxcl13_colocalization()
        cdkn1a_results = self.analyze_cdkn1a_caf_colocalization()
        
        # 6. Generate visualizations
        print("\n6. Generating visualizations...")
        self.plot_all_samples_il1r1_cxcl13_with_he()
        self.plot_all_samples_cdkn1a_caf_with_he()
        self.plot_comparative_analysis()
        
        # 7. Save results
        print("\n7. Saving results...")
        self.save_all_results()
        
        print("\n" + "=" * 70)
        print("Comprehensive Analysis Complete!")
        print("=" * 70)
        
        # Summary
        print("\nGenerated Output Files:")
        print("  - enhanced_all_samples_with_tumor_and_HE.pdf")
        print("  - enhanced_all_samples_CDKN1A_Macrophage_CAF_with_HE.pdf")
        print("  - comparative_spatial_analysis_summary.pdf")
        print("  - IL1R1_CAF_CXCL13_T_colocalization_results.csv")
        print("  - CDKN1A_Macrophage_CAF_colocalization_results.csv")
        
        return {
            'il1r1_cxcl13_results': il1r1_results,
            'cdkn1a_caf_results': cdkn1a_results
        }


if __name__ == "__main__":
    analyzer = ComprehensiveSpatialAnalyzer()
    results = analyzer.run_comprehensive_analysis()

