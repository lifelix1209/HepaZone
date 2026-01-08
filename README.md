# HepaZone

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**HepaZone** is a computational algorithm designed to assess and classify hepatic zonation patterns from single-cell RNA-sequencing data. By leveraging zone-specific gene expression signatures, HepaZone quantitatively scores cells and assigns them to distinct metabolic zones within liver lobules.

## 🎯 Key Features

- **Zonation Scoring**: Quantitative assessment of zonation status for individual cells
- **Layer Stratification**: Automated classification of cells into periportal, mid-zonal, and pericentral zones
- **Pattern Recognition**: Identification and labeling of cells based on their zonation expression patterns
- **Transcriptome-based**: Purely data-driven approach using scRNA-seq profiles
- **Flexible**: Compatible with multiple single-cell analysis frameworks

## 🧬 Background

Liver zonation refers to the spatial organization of hepatocytes along the porto-central axis, where cells exhibit distinct metabolic functions based on their position. HepaZone enables researchers to:

- Map zonation patterns in healthy and diseased liver tissue
- Track zonation changes during development or pathology
- Identify zone-specific cell populations
- Validate spatial transcriptomics findings

## 📦 Installation

```bash
# Via pip (recommended)
pip install hepazone

# From source
git clone https://github.com/yourusername/HepaZone.git
cd HepaZone
pip install -e .
```

## 🚀 Quick Start

---

**Keywords**: liver zonation, single-cell RNA-seq, hepatocyte, spatial transcriptomics, metabolic zonation, scRNA-seq analysis
