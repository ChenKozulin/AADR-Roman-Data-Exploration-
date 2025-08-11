import os
os.chdir('Path/to/data')

import pandas as pd
Roman_Meta_Data = pd.read_csv('Roman Metadata.csv')

#'Roman_Meta_Data' holds the contents of the CSV file as a pandas DataFrame
print(Roman_Meta_Data)
!/usr/bin/env python3

# Install admixtools if needed
!pip install admixtools

import pandas as pd
import numpy as np

# EIGENSTRAT format typically comes with 3 files:
# .geno (genotype data - binary)
# .snp (SNP information)  
# .ind (individual information)

# Read the accompanying files first
snp_file = "v62.0_1240k_public.snp"
ind_file = "v62.0_1240k_public.ind"

# Read SNP data
snps = pd.read_csv(snp_file, sep='\s+', header=None, 
                   names=['SNP_ID', 'chromosome', 'genetic_distance', 'position'])

# Read individual data  
individuals = pd.read_csv(ind_file, sep='\s+', header=None,
                         names=['Individual_ID', 'Sex', 'Population'])

print(f"Number of SNPs: {len(snps)}")
print(f"Number of individuals: {len(individuals)}")

# For the samples with damage rate data, analyze patterns
if len(samples_with_damage) > 0:
    # Plot damage rates
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(10, 6))
    plt.hist(samples_with_damage['damage_rate_numeric'], bins=15, alpha=0.7, edgecolor='green')
    plt.xlabel('Damage Rate')
    plt.ylabel('Number of Samples')
    plt.title(f'Ancient DNA Damage Rates (n={len(samples_with_damage)})')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # High damage samples (potentially very ancient)
    high_damage = samples_with_damage[samples_with_damage['damage_rate_numeric'] > 0.1]
    print(f"\nHigh damage samples (>10%): {len(high_damage)}")
    if len(high_damage) > 0:
        print("High damage sample details:")
        print(high_damage[['Genetic.ID', 'damage_rate_numeric', 'Date.mean.in.BP.in.years.before.1950.CE']].head())

# Look at damage rates vs age
import matplotlib.pyplot as plt
import numpy as np

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

# Plot age distribution
ax1.hist(metadata_df['Date.mean.in.BP.in.years.before.1950.CE'].dropna(), bins=20, alpha=0.7)
ax1.set_xlabel('Age (BP years)')
ax1.set_ylabel('Number of samples')
ax1.set_title('Age distribution of samples')

# Plot damage rate vs age
valid_data = metadata_df.dropna(subset=['damage_rate_numeric', 'Date.mean.in.BP.in.years.before.1950.CE'])
if len(valid_data) > 0:
    ax2.scatter(valid_data['Date.mean.in.BP.in.years.before.1950.CE'], 
                valid_data['damage_rate_numeric'], alpha=0.6)
    ax2.set_xlabel('Age (BP years)')
    ax2.set_ylabel('Damage rate')
    ax2.set_title('Damage rate vs Age')

plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set up the analysis
plt.style.use('default')
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Ancient DNA Population Analysis - 197 Samples', fontsize=16, fontweight='bold')

# 1. Age Distribution (BP years)
ax1 = axes[0, 0]
age_data = metadata_df['Date.mean.in.BP.in.years.before.1950.CE'].dropna()
if len(age_data) > 0:
    ax1.hist(age_data, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel('Age (BP years)')
    ax1.set_ylabel('Number of samples')
    ax1.set_title(f'Age Distribution (n={len(age_data)})')
    # Add statistics
    mean_age = age_data.mean()
    median_age = age_data.median()
    ax1.axvline(mean_age, color='red', linestyle='--', label=f'Mean: {mean_age:.0f}')
    ax1.axvline(median_age, color='orange', linestyle='--', label=f'Median: {median_age:.0f}')
    ax1.legend()

# 2. Sex Distribution
ax2 = axes[0, 1]
sex_data = metadata_df['Molecular.Sex'].value_counts()
if len(sex_data) > 0:
    colors = ['lightblue', 'lightpink', 'lightgray']
    wedges, texts, autotexts = ax2.pie(sex_data.values, labels=sex_data.index, autopct='%1.1f%%', 
                                       colors=colors[:len(sex_data)])
    ax2.set_title(f'Sex Distribution (n={sex_data.sum()})')

# 3. Geographic Distribution (Political Entity)
ax3 = axes[0, 2]
geo_data = metadata_df['Political.Entity'].value_counts().head(10)
if len(geo_data) > 0:
    geo_data.plot(kind='barh', ax=ax3, color='lightgreen')
    ax3.set_xlabel('Number of samples')
    ax3.set_title('Top 10 Political Entities')
    ax3.tick_params(axis='y', labelsize=8)

# 4. Data Quality - SNP Coverage
ax4 = axes[1, 0]
snp_data = metadata_df['autosomal.SNPs..1240k.snpset'].dropna()
if len(snp_data) > 0:
    ax4.hist(snp_data, bins=20, alpha=0.7, color='lightcoral', edgecolor='black')
    ax4.set_xlabel('Number of autosomal SNPs (1240k)')
    ax4.set_ylabel('Number of samples')
    ax4.set_title(f'SNP Coverage Distribution (n={len(snp_data)})')
    ax4.set_xscale('log')

# 5. Damage Rate Distribution
ax5 = axes[1, 1]
damage_numeric = metadata_df['damage_rate_numeric'].dropna()
if len(damage_numeric) > 0:
    ax5.hist(damage_numeric, bins=15, alpha=0.7, color='gold', edgecolor='black')
    ax5.set_xlabel('Damage Rate')
    ax5.set_ylabel('Number of samples')
    ax5.set_title(f'Damage Rate Distribution (n={len(damage_numeric)})')

# 6. Age vs Damage Rate
ax6 = axes[1, 2]
combined_data = metadata_df.dropna(subset=['Date.mean.in.BP.in.years.before.1950.CE', 'damage_rate_numeric'])
if len(combined_data) > 5:
    scatter = ax6.scatter(combined_data['Date.mean.in.BP.in.years.before.1950.CE'], 
                         combined_data['damage_rate_numeric'], 
                         alpha=0.6, color='purple')
    ax6.set_xlabel('Age (BP years)')
    ax6.set_ylabel('Damage Rate')
    ax6.set_title(f'Age vs Damage Rate (n={len(combined_data)})')
    
    # Add correlation
    correlation = combined_data['Date.mean.in.BP.in.years.before.1950.CE'].corr(combined_data['damage_rate_numeric'])
    ax6.text(0.05, 0.95, f'Correlation: {correlation:.3f}', transform=ax6.transAxes, 
             bbox=dict(boxstyle="round", facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

# Print detailed statistics
print("=== POPULATION COMPOSITION ANALYSIS ===\n")

print("1. SAMPLE OVERVIEW:")
print(f"   Total samples: {len(metadata_df)}")
print(f"   Samples with age data: {metadata_df['Date.mean.in.BP.in.years.before.1950.CE'].notna().sum()}")
print(f"   Samples with damage rates: {metadata_df['damage_rate_numeric'].notna().sum()}")
print(f"   Samples with SNP data: {metadata_df['autosomal.SNPs..1240k.snpset'].notna().sum()}")

print("\n2. TEMPORAL DISTRIBUTION:")
age_data = metadata_df['Date.mean.in.BP.in.years.before.1950.CE'].dropna()
if len(age_data) > 0:
    print(f"   Age range: {age_data.min():.0f} - {age_data.max():.0f} BP years")
    print(f"   Mean age: {age_data.mean():.0f} BP years")
    print(f"   Median age: {age_data.median():.0f} BP years")
    
    # Time periods
    def categorize_age(bp_years):
   
        if 1500 < bp_years < 1600:
            return "5c.AD (1500 BP)"
        elif 1600 < bp_years < 1700:
            return "4c.AD (1600 BP)"
        elif 1700 < bp_years < 1800:
            return "3c.AD (1700 BP)"
        elif 1800 < bp_years < 1900:
            return "2c.AD (1800 BP)"
        elif 1900 < bp_years < 2000:
            return "1c.AD (1800 BP)"
        else:
            return "BC (>2000 BP)"
    
    age_categories = age_data.apply(categorize_age).value_counts()
    print("\n   Time period distribution:")
    for period, count in age_categories.items():
        print(f"   {period}: {count} ({count/len(age_data)*100:.1f}%)")

print("\n3. SEX DISTRIBUTION:")
sex_counts = metadata_df['Molecular.Sex'].value_counts()
for sex, count in sex_counts.items():
    print(f"   {sex}: {count} ({count/len(metadata_df)*100:.1f}%)")

print("\n4. GEOGRAPHIC DISTRIBUTION:")
print(f"   Unique political entities: {metadata_df['Political.Entity'].nunique()}")
geo_top = metadata_df['Political.Entity'].value_counts().head(5)
print("   Top 5 regions:")
for region, count in geo_top.items():
    print(f"   {region}: {count} ({count/len(metadata_df)*100:.1f}%)")

print("\n5. DATA QUALITY METRICS:")
snp_data = metadata_df['autosomal.SNPs..1240k.snpset'].dropna()
if len(snp_data) > 0:
    print(f"   SNP coverage range: {snp_data.min():.0f} - {snp_data.max():.0f}")
    print(f"   Mean SNP coverage: {snp_data.mean():.0f}")
    
    # Quality categories
    def categorize_snps(snp_count):
        if snp_count > 100000:
            return "High (>100k)"
        elif snp_count > 50000:
            return "Medium-High (50k-100k)"
        elif snp_count > 10000:
            return "Medium (10k-50k)"
        else:
            return "Low (<10k)"
    
    snp_categories = snp_data.apply(categorize_snps).value_counts()
    print("\n   SNP coverage categories:")
    for category, count in snp_categories.items():
        print(f"   {category}: {count} ({count/len(snp_data)*100:.1f}%)")

damage_data = metadata_df['damage_rate_numeric'].dropna()
if len(damage_data) > 0:
    print(f"\n   Damage rate range: {damage_data.min():.4f} - {damage_data.max():.4f}")
    print(f"   Mean damage rate: {damage_data.mean():.4f}")

print("\n6. CORRELATIONS:")
numeric_cols = ['Date.mean.in.BP.in.years.before.1950.CE', 'damage_rate_numeric', 
                'autosomal.SNPs..1240k.snpset', 'Lat.', 'Long.']
print("   Correlation matrix (key variables):")
corr_data = metadata_df[numeric_cols].corr()
print(corr_data.round(3))
