import pandas as pd

"""
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

    df = pd.DataFrame(
        0,
        index=hugo_list,
        columns=clonal_summary.columns
    )

    clonal_summary = pd.concat([clonal_summary, df])

    clonal_summary = clonal_summary.sort_index()

    return clonal_summary