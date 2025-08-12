import streamlit as st
import pandas as pd
import pickle as pkl
import numpy as np
from scipy.stats import nbinom

@st.cache_data
def load_and_process_data():
    """Loads and processes both datasets."""
    try:
        # Load
        data_binom = pkl.load(open("data/mc3_binom_left.pkl", "rb")).reset_index(drop=True)
        # data_binom = data_binom.reset_index(drop=True)
        # data_binom["Clonality"] = data_binom["Clonality"].astype(str)

        data_ztest = pkl.load(open("data/mc3_Ztest_left.pkl", "rb")).reset_index(drop=True)
        # data_ztest = data_ztest.reset_index(drop=True)
        # data_ztest["Clonality"] = data_ztest["Clonality"].astype(str)

        # Process (Group by Gene)
        summary_binom = data_binom.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
        summary_binom.columns.name = None
        summary_binom = summary_binom.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"}).reset_index()

        summary_ztest = data_ztest.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
        summary_ztest.columns.name = None
        summary_ztest = summary_ztest.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"}).reset_index()

        # Ensure both count columns exist
        if "Clonal Count" not in summary_binom.columns: summary_binom["Clonal Count"] = 0
        if "Sub-Clonal Count" not in summary_binom.columns: summary_binom["Sub-Clonal Count"] = 0
        if "Clonal Count" not in summary_ztest.columns: summary_ztest["Clonal Count"] = 0
        if "Sub-Clonal Count" not in summary_ztest.columns: summary_ztest["Sub-Clonal Count"] = 0

        return data_binom, data_ztest, summary_binom, summary_ztest
    except FileNotFoundError:
        st.error("Error: Data files not found in 'data/' directory (e.g., 'data/mc3_binom_left.pkl'). Please ensure the files exist.")
        return None, None, None, None
    except Exception as e:
        st.error(f"An error occurred during data loading or processing: {e}")
        return None, None, None, None

def extend_clonal_summary(clonal_summary):
    """
    This function extends the clonal summary DataFrame by adding rows for genes with zero mutations.
    It uses a list of HUGO symbols from a CSV file to ensure that all genes are represented.

    The function assume that the data has already gone through the following steps:
    # Load the data from a pickle file
    data = pd.read_pickle("path/to/your/data.pkl")
    data["Clonality"] = data["Clonality"].astype(str)

    # Group by Hugo_Symbol and count occurrences of Clonality
    clonal_summary = data.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)

    # Rename columns for clarity
    clonal_summary.columns.name = None  # Remove column index name
    clonal_summary = clonal_summary.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"})
    """
    hugo_symbols = pd.read_csv("data/zero_mutation_genes.csv")
    hugo_list = hugo_symbols['hugo_symbol'].str.strip().tolist()

    df = pd.DataFrame({
        'Hugo_Symbol': hugo_list,
        'Clonal Count': [0] * len(hugo_list),
        'Sub-Clonal Count': [0] * len(hugo_list)
    })

    # Concatenate and reset index
    clonal_summary = pd.concat([clonal_summary, df], ignore_index=True)

    return clonal_summary

@st.cache_data
def calculate_summary_stats(clonal_summary):
    """Calculates summary statistics from the gene summary dataframe."""
    if clonal_summary is None or clonal_summary.empty:
        return pd.DataFrame()

    from src.config import CONFIG
    ndecimal = CONFIG['ndecimal']
    
    summary_stats = pd.DataFrame({
        "Statistic": ["Mean", 'Median', "Variance"],
        "Clonal": [
            clonal_summary["Clonal Count"].mean().round(ndecimal),
            clonal_summary["Clonal Count"].median().round(ndecimal),
            clonal_summary["Clonal Count"].var().round(ndecimal)
        ],
        "Sub-Clonal": [
            clonal_summary["Sub-Clonal Count"].mean().round(ndecimal),
            clonal_summary["Sub-Clonal Count"].median().round(ndecimal),
            clonal_summary["Sub-Clonal Count"].var().round(ndecimal)
        ]
    })
    return summary_stats

@st.cache_data
def add_overlap_counts(summary):
    """Adds overlap count and log(overlap count) columns to the input DataFrame."""
    summary = summary.copy()
    overlap_counts = summary.groupby(['Clonal Count', 'Sub-Clonal Count']).size().reset_index(name='Overlap Count')
    overlap_counts['Log Overlap Count'] = np.log(overlap_counts['Overlap Count'])
    return overlap_counts 

@st.cache_data
def negbin_quantiles(Y_pred: pd.Series, scale: float ,q_range=0.95):
    """
    Compute (1 - alpha)% quantile interval of a negative binomial 
    given predicted means (mu) and dispersion alpha.

    Parameters:
    - Y_pred: pd.Series of predicted means (μ values)
    - scale: dispersion parameter (1/theta), where Var(Y) = μ + μ²/theta
    - alpha: float, defines the range for quantiles (default is 0.05 for 95% CI)

    Returns:
    - pd.DataFrame with lower and upper quantiles
    """
    theta = 1 / scale  # Dispersion parameter for Negative Binomial
    alpha = 1 - q_range  # Convert to alpha for quantiles
    p = theta / (theta + Y_pred)
    lower = nbinom.ppf(alpha / 2, theta, p)
    upper = nbinom.ppf(1 - alpha / 2, theta, p)
    return lower, upper