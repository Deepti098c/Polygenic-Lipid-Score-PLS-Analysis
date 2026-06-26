# ===============================
# Cleaned PLS Calculation Pipeline
# ===============================

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error

# ----------------------------
# Step 1: Load and filter SNPs
# ----------------------------
def load_significant_snps(snp_file, p_threshold=1e-5):
    snp_data = pd.read_csv(snp_file)
    if 'GeneID' not in snp_data.columns:
        snp_data['GeneID'] = snp_data['GeneName']

    # Aggregate duplicate SNP entries
    snp_data = snp_data.groupby(['POS', 'REF', 'ALT'], as_index=False).agg({
        'Effect':'mean',
        'SE':'mean',
        'GeneID': lambda x: ';'.join(x.astype(str)),
        'PC128.FarmCPU':'mean'
    })

    significant_snps = snp_data[snp_data['PC128.FarmCPU'] < p_threshold]
    if significant_snps.empty:
        raise ValueError(f"No significant SNPs found at p < {p_threshold}")
    print(f"Significant SNPs: {len(significant_snps)}")
    return significant_snps

# ----------------------------
# Step 2: Load HMP genotypes
# ----------------------------
def load_genotypes(hmp_file, significant_snps):
    hmp = pd.read_csv(hmp_file, sep='\t')
    # Keep only SNPs in significant list
    hmp = hmp[hmp['pos'].isin(significant_snps['POS'])]
    
    # Extract line columns
    ignore_cols = ['s#','alleles','chrom','pos','strand','assembly#','center',
                   'protLSID','assayLSID','panelLSID','QCcode']
    line_cols = [c for c in hmp.columns if c not in ignore_cols]
    genotypes = hmp[['pos'] + line_cols].set_index('pos').T
    genotypes.index.name = 'Line'
    
    # Map alleles to 0/1/2 dosage
    snp_dict = significant_snps.set_index('POS')[['REF','ALT','Effect']].to_dict('index')
    
    def allele_to_dosage(x, pos):
        if pd.isna(x): 
            return np.nan
        allele = str(x)[0]
        ref = snp_dict[pos]['REF']
        alt = snp_dict[pos]['ALT']
        if allele == ref:
            return 0
        elif allele == alt:
            return 2
        else:
            return 1

    for pos in genotypes.columns:
        genotypes[pos] = genotypes[pos].apply(lambda x: allele_to_dosage(x, pos))
    
    print(f"Genotypes loaded: {genotypes.shape[0]} lines, {genotypes.shape[1]} SNPs")
    return genotypes, snp_dict

# ----------------------------
# Step 3: Compute PLS
# ----------------------------
def compute_pls(genotypes, snp_dict):
    effect_sizes = np.array([snp_dict[pos]['Effect'] for pos in genotypes.columns])
    # Fill missing genotypes with mean of SNP column
    genotypes_filled = genotypes.fillna(genotypes.mean())
    pls_scores = genotypes_filled.dot(effect_sizes)
    pls_df = pd.DataFrame({'Line': genotypes.index, 'PLS': pls_scores})
    print("PLS scores computed.")
    return pls_df

# ----------------------------
# Step 4: Regression & predicted lipid content
# ----------------------------
def regress_pls(pls_df, phenotype_file, pheno_column='PC128'):
    pheno = pd.read_csv(phenotype_file, sep='\t').rename(columns={'Taxa':'Line'})
    data = pls_df.merge(pheno[['Line', pheno_column]], on='Line', how='inner')
    
    X = data[['PLS']].values
    y = data[pheno_column].values
    
    model = LinearRegression().fit(X, y)
    y_pred = model.predict(X)
    
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    
    print(f"Regression completed: R² = {r2:.3f}, RMSE = {rmse:.3f}")
    
    predicted_df = pd.DataFrame({'Line': data['Line'], 'Predicted_Lipid_Content': y_pred})
    return model, predicted_df

# ----------------------------
# Step 5: Count beneficial alleles
# ----------------------------
def count_beneficial_alleles(genotypes, snp_dict):
    # Count alleles in same direction as effect
    count_dict = {}
    for pos in genotypes.columns:
        effect = snp_dict[pos]['Effect']
        if effect > 0:
            count_dict[pos] = (genotypes[pos] == 2).astype(int)
        else:
            count_dict[pos] = (genotypes[pos] == 0).astype(int)
    count_df = pd.DataFrame(count_dict)
    count_df['Line'] = genotypes.index
    count_df['Beneficial_Allele_Count'] = count_df.drop(columns='Line').sum(axis=1)
    return count_df[['Line', 'Beneficial_Allele_Count']]

# ----------------------------
# Step 6: Rank lines & annotate top candidates
# ----------------------------
def rank_and_label(pls_df, predicted_df, beneficial_df):
    final = pls_df.merge(predicted_df, on='Line').merge(beneficial_df, on='Line')
    final['Rank'] = final['PLS'].rank(ascending=False, method='min')
    
    # Top 10% beneficial alleles + high PLS
    threshold = final['Beneficial_Allele_Count'].quantile(0.9)
    final['Notes'] = np.where(
        (final['PLS'] > final['PLS'].quantile(0.66)) & 
        (final['Beneficial_Allele_Count'] >= threshold),
        'Top candidate',
        'Low priority'
    )
    return final.sort_values('Rank')

# ----------------------------
# Main Execution
# ----------------------------
snp_file = "GWAS_SNP_gene_50kb.csv"
hmp_file = "SAP.cleaned.hmp.txt"
phenotype_file = "Oil.phe"

snps = load_significant_snps(snp_file)
genotypes, snp_dict = load_genotypes(hmp_file, snps)
pls_df = compute_pls(genotypes, snp_dict)
model, predicted_df = regress_pls(pls_df, phenotype_file)
beneficial_df = count_beneficial_alleles(genotypes, snp_dict)
final_df = rank_and_label(pls_df, predicted_df, beneficial_df)

final_df.to_csv("lipid_analysis_final.csv", index=False)
print("PLS pipeline complete! Top candidates:")
print(final_df[final_df['Notes']=='Top candidate'].head(10))
