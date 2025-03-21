import panel as pn
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

pn.extension('tabulator')

# === Load data ===
data = pkl.load(open("data/mc3_binom.pkl", "rb"))
data = data.reset_index(drop=True)
data["Clonality"] = data["Clonality"].astype(str)

# === Clonality counts table ===
clonal_summary = data.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
clonal_summary.columns.name = None
clonal_summary = clonal_summary.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"})

# === DataFrame display ===
data_table = pn.widgets.Tabulator(data.head(100), pagination='remote', page_size=10, width=1000)

# === Clonality counts plot ===
clonality_counts = data["Clonality"].value_counts()
fig1, ax1 = plt.subplots()
clonality_counts.plot(kind='bar', color=['#4CAF50', '#FF5722'], ax=ax1)
ax1.set_title("Clonality Counts")
ax1.set_xlabel("Clonality")
ax1.set_ylabel("Count")
bar_fig = pn.pane.Matplotlib(fig1, tight=True)
# plt.close(fig1)

# === Hexplot visualization ===
x = clonal_summary["Clonal Count"]
y = clonal_summary["Sub-Clonal Count"]
x_with_const = sm.add_constant(x)
poisson_model = sm.GLM(y, x_with_const, family=sm.families.Poisson()).fit()

x_plot = np.linspace(0, 400, 400)
x_step5 = np.arange(0, 400, 5)
y_pred = poisson_model.predict(sm.add_constant(x_plot))

group_means = [clonal_summary[(x >= i) & (x < i + 5)]["Sub-Clonal Count"].mean() for i in x_step5]

fig2, ax2 = plt.subplots(figsize=(10, 6))
hb = ax2.hexbin(x, y, gridsize=300, cmap="Purples", mincnt=1)
ax2.plot(x_plot, y_pred, color='red', label='Poisson regression')
ax2.plot(x_step5, group_means, color='blue', label='Step-5 Mean')
fig2.colorbar(hb, ax=ax2, label='Count in Bin')
ax2.set_xlim(0, 400)
ax2.set_ylim(0, 2000)
ax2.set_xlabel("Clonal Count")
ax2.set_ylabel("Sub-Clonal Count")
ax2.set_title("Hexbin of Clonal vs Sub-Clonal Counts")
ax2.legend()
hex_fig = pn.pane.Matplotlib(fig2, tight=True)
# plt.close(fig2)

# === Layout ===
dashboard = pn.Column(
    "# 🧬 Tumor Clonality Dashboard",
    "## Overview of the dataset",
    data_table,
    "## Clonality value counts",
    bar_fig,
    "## Clonal vs Sub-Clonal Count Hexplot",
    hex_fig
)

dashboard.servable()