import pandas as pd
import numpy as np
import pickle as pkl
from tqdm import tqdm
from scipy.stats import norm
from statsmodels.stats.proportion import binom_test

# Path to your MAF file
maf_file = "data/mc3.maf"

# Read the MAF file into a DataFrame.
# The 'comment' parameter skips any header lines that start with '#' (if present)
df = pd.read_csv(maf_file, sep="\t", comment="#", low_memory=False)

# Display the first few rows of the DataFrame
print(df.head())

# Save the DataFrame to a pickle file
pkl.dump(df, open("data/mc3.pkl", "wb"))

df = pkl.load(open("data/mc3.pkl", "rb"))
df.head()

# n_ref_count: The number of reads supporting the reference allele
# n_alt_count: The number of reads supporting the alternate allele
# t_alt_count: The number of reads supporting the alternate allele in the tumor sample

columns = ["Hugo_Symbol", "t_alt_count", "t_ref_count", "Segment_Mean", "Purity"]
data = df[columns]
del df

# Remove rows with missing values in any of the columns
# Remove rows with Segment_Mean values outside the range [-0.3, 0.3], cuz we only want to work with allele having copy number 2
length0 = len(data)
data = data.dropna()
data = data[(data["Segment_Mean"]>=-0.3) & (data["Segment_Mean"]<=0.3)]
length1 = len(data)
print(f"Removed {length0 - length1} rows with missing values or Segment_Mean outside the range [-0.3, 0.3]")

# Chunk with longer excution time - the processed data will be saved to a pickle file
# Calculate Variant Allele Frequency (VAF)
# data["VAF"] = data["t_alt_count"] / (data["t_alt_count"] + data["t_ref_count"])

# initialize the tqdm progress bar
tqdm.pandas()

# Calculate the coverage of the mutation 
data['Coverage'] = data['t_alt_count'] + data['t_ref_count']

# Classify mutations using classic VAF v.s purity threshold
# data["Clonality0"] = data.apply(
#     lambda row: "Clonal" if row["VAF"] > row["Purity"] * 0.5 else "Sub-Clonal",
#     axis=1
# )

############################################# calculate using left tail Z test #############################################
# Calculate the Z Test p-value
data["Z_stat"] = ((data['t_alt_count'] / data['Coverage']) - data["Purity"] * 0.5) / ((data["Purity"] * 0.5) * (1 - data["Purity"] * 0.5) / data["Coverage"])**0.5
data["Z_pval"] = norm.cdf(data["Z_stat"])
data["Z_pval"] = data["Z_pval"].fillna(1)  # Fill NaN values with 1

# Classify mutations using a binomial test
data["Clonality"] = data.progress_apply(
    lambda row: "Sub-Clonal" if row["Z_pval"] < 0.05 else "Clonal",
    axis=1
)
# Save the DataFrame to a pickle file
pkl.dump(data, open("data/mc3_Ztest_left.pkl", "wb"))
#####################################################################################################################################

############################################# calculate using left tail binomial test #############################################
# Calculate the Binomial Test p-value
data["Binom_pval"] = data.progress_apply(lambda row: binom_test(row["t_alt_count"], row["Coverage"], row["Purity"]*0.5, alternative="smaller"), axis=1)

# Classify mutations using a binomial test
data["Clonality"] = data.progress_apply(
    lambda row: "Sub-Clonal" if row["Binom_pval"] < 0.05 else "Clonal",
    axis=1
)
# Save the DataFrame to a pickle file
pkl.dump(data, open("data/mc3_binom_left.pkl", "wb"))
#####################################################################################################################################

############################################# calculate using right tail binomial test #############################################
# Calculate the Binomial Test p-value
data["Binom_pval"] = data.progress_apply(lambda row: binom_test(row["t_alt_count"], row["Coverage"], row["Purity"]*0.5, alternative="larger"), axis=1)

# Classify mutations using a binomial test
data["Clonality"] = data.progress_apply(
    lambda row: "Clonal" if row["Binom_pval"] < 0.05 else "Sub-Clonal",
    axis=1
)
# Save the DataFrame to a pickle file
pkl.dump(data, open("data/mc3_binom_right.pkl", "wb"))
#####################################################################################################################################