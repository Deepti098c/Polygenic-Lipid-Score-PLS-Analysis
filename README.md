# Polygenic-Lipid-Score-PLS-Analysis
# PLS-Based Allele Identification Pipeline

A Python pipeline for identifying elite breeding lines using **Genome-Wide Association Study (GWAS)** results, **HapMap genotype data**, and **phenotypic lipid content**. The pipeline computes **Polygenic Scores (PLS)**, predicts lipid content, counts beneficial alleles, and ranks breeding lines for selection.

------
GWAS
      \
       \
        --> PLS Calculation --> Beneficial Alleles --> Ranking --> Elite Lines
       /
HapMap/
      \
Phenotype

-------
## Overview

This pipeline integrates GWAS effect sizes with genotype data to estimate the genetic potential of breeding lines. It performs the following analyses:

- Filters significant SNPs from GWAS results
- Converts HapMap genotypes into allele dosage values
- Computes Polygenic Scores (PLS)
- Predicts lipid content using linear regression
- Counts beneficial alleles for each breeding line
- Ranks breeding lines according to genetic potential
- Identifies elite candidate lines for breeding programs

---


# Features

- Automatic filtering of significant GWAS SNPs
- Duplicate SNP aggregation
- HapMap genotype processing
- Allele dosage conversion
- Polygenic score calculation
- Linear regression for phenotype prediction
- Beneficial allele counting
- Automatic ranking of breeding lines
- Candidate line annotation

---

# Requirements

- Python 3.8 or newer

## Python packages

```bash
pip install pandas numpy scikit-learn
```

---

# Input Files

## 1. GWAS SNP File

Example:

```
GWAS_SNP_gene_50kb.csv
```

Required columns:

| Column | Description |
|----------|-------------|
| POS | SNP position |
| REF | Reference allele |
| ALT | Alternate allele |
| Effect | SNP effect size |
| SE | Standard error |
| PC128.FarmCPU | GWAS p-value |
| GeneID or GeneName | Associated gene |

---

## 2. HapMap Genotype File

Example:

```
SAP.cleaned.hmp.txt
```

Required fields:

- SNP position (`pos`)
- genotype columns for each breeding line

---

## 3. Phenotype File

Example:

```
Oil.phe
```

Required columns:

| Column | Description |
|----------|-------------|
| Taxa | Line name |
| PC128 | Lipid content phenotype |

---

# Usage

Update the file names if needed:

```python
snp_file = "GWAS_SNP_gene_50kb.csv"
hmp_file = "SAP.cleaned.hmp.txt"
phenotype_file = "Oil.phe"
```

Run the script:

```bash
python PLS_Alleles_identification.py
```

---

# Output

The pipeline generates:

```
lipid_analysis_final.csv
```

Output columns:

| Column | Description |
|----------|-------------|
| Line | Breeding line |
| PLS | Polygenic Score |
| Predicted_Lipid_Content | Predicted phenotype |
| Beneficial_Allele_Count | Number of favorable alleles |
| Rank | Rank based on PLS |
| Notes | Candidate classification |

---

# Beneficial Allele Definition

For each SNP:

- **Positive effect size** → ALT allele is considered beneficial.
- **Negative effect size** → REF allele is considered beneficial.

The pipeline counts the number of beneficial alleles present in each breeding line.

---

# Candidate Selection

A breeding line is labeled **Top candidate** if it satisfies:

- PLS greater than the 66th percentile
- Beneficial allele count within the top 10% of the population

All other lines are labeled **Low priority**.

---


# Dependencies

- pandas
- numpy
- scikit-learn


---

# Citation

If you use this pipeline in your research, please cite the associated publication or acknowledge this repository.

