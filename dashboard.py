import panel as pn
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

pn.extension('tabulator', sizing_mode='stretch_width')
# Set default font sizes globally for all plots
plt.rcParams.update({
    'axes.titlesize': 10,    # Titles
    'axes.labelsize': 10,    # X and Y labels
    'xtick.labelsize': 9,    # X-axis tick labels
    'ytick.labelsize': 9,    # Y-axis tick labels
    'legend.fontsize': 9,    # Legend
    'font.size': 10          # General font size (fallback/default)
})
figsize = (5, 4)

# === Load data ===
data_left = pkl.load(open("data/mc3_binom_left.pkl", "rb")).reset_index(drop=True)
data_left["Clonality"] = data_left["Clonality"].astype(str)

data_right = pkl.load(open("data/mc3_binom_right.pkl", "rb")).reset_index(drop=True)
data_right["Clonality"] = data_right["Clonality"].astype(str)

# Function to process each dataset
def process_data(data):
    summary = data.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
    summary.columns.name = None
    summary = summary.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"})
    return summary

clonal_summary_left = process_data(data_left)
clonal_summary_right = process_data(data_right)

# Function to create clonality counts plot
def plot_clonality_counts(data, title):
    orders = ['Clonal', 'Sub-Clonal']
    counts = data["Clonality"].value_counts().reindex(orders, fill_value=0)
    fig, ax = plt.subplots(figsize=figsize)
    counts.plot(kind='bar', color=['#6FCF97', '#EB5757'], ax=ax)
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel("Clonality")
    ax.set_ylabel("Count")
    plt.xticks(rotation=0)
    return pn.pane.Matplotlib(fig, tight=True)

# Function to create pie chart
def pie_chart(data, title):
    orders = data["Hugo_Symbol"].value_counts().index
    counts = data["Hugo_Symbol"].value_counts().reindex(orders, fill_value=0)
    fig, ax = plt.subplots(figsize=figsize)
    counts.plot(kind='pie', ax=ax, autopct='%1.1f%%', startangle=90, colors=['#6FCF97', '#EB5757'])
    ax.set_title(title, fontweight='bold')
    return pn.pane.Matplotlib(fig, tight=True)

# Function to fit Poisson regression
def fit_poisson_model(summary):
    x = summary["Clonal Count"]
    y = summary["Sub-Clonal Count"]
    x_const = sm.add_constant(x)
    poisson_model = sm.GLM(y, x_const, family=sm.families.Poisson()).fit()
    return poisson_model

# Function to fit Negative Binomial regression
def fit_negbin_model(summary):
    x = summary["Clonal Count"]
    y = summary["Sub-Clonal Count"]
    x_const = sm.add_constant(x)
    negbin_model = sm.NegativeBinomial(y, x_const).fit()
    return negbin_model

# Function for hexplot visualization
def plot_hexbin_poisson(summary, poisson_model, title, lims=(400, 400)):
    x = summary["Clonal Count"]
    y = summary["Sub-Clonal Count"]

    xlim, ylim = lims

    x_plot = np.linspace(0, xlim, 400)
    # step = 1
    # x_step = np.arange(0, xlim, step)
    y_pred = poisson_model.predict(sm.add_constant(x_plot))

    # group_means = [summary[(x >= i) & (x < i + step)]["Sub-Clonal Count"].mean() for i in x_step]

    fig, ax = plt.subplots(figsize=figsize)
    # ax.scatter(x_step, group_means, color='#2F80ED', marker='x', label='Mean')
    hb = ax.hexbin(x, y, gridsize=2000, bins='log', cmap="Blues", mincnt=1)
    ax.plot(x_plot, y_pred, color='#FF6F61', linewidth=2, label='Poisson Prediction')
    fig.colorbar(hb, ax=ax, label='Log Count in Bin')
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.set_xlabel("Clonal Count")
    ax.set_ylabel("Sub-Clonal Count")
    ax.set_title(title, fontweight='bold')
    ax.legend()
    return pn.pane.Matplotlib(fig)

# Fit Poisson models
poisson_model_left = fit_poisson_model(clonal_summary_left)
poisson_model_right = fit_poisson_model(clonal_summary_right)

# Fit Negative Binomial models
negbin_model_left = fit_negbin_model(clonal_summary_left)
negbin_model_right = fit_negbin_model(clonal_summary_right)

# Create visualizations
bar_fig_left = plot_clonality_counts(data_left, "Left Data Clonality Counts")
bar_fig_right = plot_clonality_counts(data_right, "Right Data Clonality Counts")

# pie_fig_left = pie_chart(data_left, "Left Data: Hugo Symbol Distribution")
# pie_fig_right = pie_chart(data_right, "Right Data: Hugo Symbol Distribution")

hex_fig_left = plot_hexbin_poisson(clonal_summary_left, poisson_model_left, "Left Data: Clonal vs Sub-Clonal", lims=(300, 100))
hex_fig_right = plot_hexbin_poisson(clonal_summary_right, poisson_model_right, "Right Data: Clonal vs Sub-Clonal", lims=(100, 500))

hex_fig_left_negbin = plot_hexbin_poisson(clonal_summary_left, negbin_model_left, "Left Data: Clonal vs Sub-Clonal (Negative Binomial)", (300, 100))
hex_fig_right_negbin = plot_hexbin_poisson(clonal_summary_right, negbin_model_right, "Right Data: Clonal vs Sub-Clonal (Negative Binomial)", lims=(100, 500))

# === Poisson Model Summaries ===
summary_left = pn.pane.Markdown(f"### Left Data Poisson Regression\n```\n{poisson_model_left.summary()}\n```")
summary_right = pn.pane.Markdown(f"### Right Data Poisson Regression\n```\n{poisson_model_right.summary()}\n```")

# === Negative Binomial Model Summaries ===
summary_left_negbin = pn.pane.Markdown(f"### Left Data Negative Binomial Regression\n```\n{negbin_model_left.summary()}\n```")
summary_right_negbin = pn.pane.Markdown(f"### Right Data Negative Binomial Regression\n```\n{negbin_model_right.summary()}\n```")

# === Data tables ===
data_table_left = pn.widgets.Tabulator(data_left.head(int(len(data_left)*0.3)), pagination='remote', page_size=10, height=300)
data_table_right = pn.widgets.Tabulator(data_right.head(int(len(data_right)*0.3)), pagination='remote', page_size=10, height=300)

# # === Summary sets mean, median and variance ===
# summary_stats_left = pd.DataFrame({
#     "Statistic": ["Mean", 'Median', "Variance"],
#     "Clonal Summary Left": [
#         clonal_summary_left.mean().round(2), 
#         clonal_summary_left.median().round(2), 
#         clonal_summary_left.var().round(2)
#     ],
#     "Clonal Summary Right": [
#         clonal_summary_right.mean().round(2), 
#         clonal_summary_right.median().round(2), 
#         clonal_summary_right.var().round(2)
#     ]
# })

ndecimal = 2

summary_stats_left = pd.DataFrame({
    "Statistic": ["Mean", 'Median', "Variance"],
    "Clonal": [
        clonal_summary_left["Clonal Count"].mean().round(ndecimal),
        clonal_summary_left["Clonal Count"].median().round(ndecimal),
        clonal_summary_left["Clonal Count"].var().round(ndecimal)
    ],
    "Sub-Clonal": [
        clonal_summary_left["Sub-Clonal Count"].mean().round(ndecimal),
        clonal_summary_left["Sub-Clonal Count"].median().round(ndecimal),
        clonal_summary_left["Sub-Clonal Count"].var().round(ndecimal)
    ]
})

summary_stats_right = pd.DataFrame({
    "Statistic": ["Mean", 'Median', "Variance"],
    "Clonal": [
        clonal_summary_right["Clonal Count"].mean().round(ndecimal),
        clonal_summary_right["Clonal Count"].median().round(ndecimal),
        clonal_summary_right["Clonal Count"].var().round(ndecimal)
    ],
    "Sub-Clonal": [
        clonal_summary_right["Sub-Clonal Count"].mean().round(ndecimal),
        clonal_summary_right["Sub-Clonal Count"].median().round(ndecimal),
        clonal_summary_right["Sub-Clonal Count"].var().round(ndecimal)
    ]
})

stats_table_left = pn.widgets.Tabulator(summary_stats_left, pagination='remote', page_size=10, height=300)
stats_table_right = pn.widgets.Tabulator(summary_stats_right, pagination='remote', page_size=10, height=300)

# === Dashboard Layout ===
title = pn.pane.Markdown("# 🧬 Tumor Clonality Dashboard", margin=(0,0,20,0), styles={'color':'#2F80ED'})

header_style = {'background': '#F0F8FF', 'padding': '10px', 'border-radius': '5px'}
comment_style = {'background': '#e6ffe6', 'padding': '10px', 'border-radius': '5px', 'font-size': '16px'}

dashboard = pn.Column(
    title,
    pn.Row(
        pn.Column(pn.pane.Markdown("### Left Dataset Overview", styles=header_style), data_table_left),
        pn.Column(pn.pane.Markdown("### Right Dataset Overview", styles=header_style), data_table_right)
    ),
    # pn.pane.Markdown("## Basic Data Statistics", styles=header_style),
    # pn.widgets.Tabulator(summary_stats, pagination='remote', page_size=10, height=300),
    pn.Row(
        pn.Column(pn.pane.Markdown("### Left Data Summary Statistics", styles=header_style), stats_table_left),
        pn.Column(pn.pane.Markdown("### Right Data Summary Statistics", styles=header_style), stats_table_right),
        pn.pane.Markdown(
            "Both datasets have clear evidence of **overdispersion**, as the variance is much larger than the mean. "
            "This is expected in count data, and is why we are using **Poisson** and **Negative Binomial** regression models.",
            styles=comment_style,
            width=300 # Adjust width to fit text
            )
    ),
    # pn.pane.Markdown("## Hugo Symbol Distribution", styles=header_style),
    # pn.Row(pie_fig_left, pie_fig_right),
    pn.pane.Markdown("## Clonality Counts", styles=header_style),
    pn.Row(
        bar_fig_left, 
        bar_fig_right,
        pn.pane.Markdown(
            "The bar plots on the left show the counts of clonal and sub-clonal mutations in the datasets. "
            "The left dataset has more clonal mutations, while the right dataset has more sub-clonal mutations due to different settings of $H_0$.",
            styles=comment_style,
            width=300 # Adjust width to fit text
            )
        ),
    pn.pane.Markdown("## Poisson Regression Model Summaries", styles=header_style),
    pn.Row(
        summary_left, 
        summary_right,
        pn.pane.Markdown(
            "The Poisson regression models above show a significant relationship between clonal and sub-clonal counts. "
            "Even though the data is overdispersed, the **Pseudo-R-square** values are misleadingly high, **WHY**?",
            styles=comment_style,
            width=300 # Adjust width to fit text
            ),
        ),
    pn.pane.Markdown("## Poisson Regression Hexplots Comparison", styles=header_style),
    pn.Row(hex_fig_left, hex_fig_right),
    pn.pane.Markdown("## Negative Binomial Regression Model Summaries", styles=header_style),
    pn.Row(
        summary_left_negbin, 
        summary_right_negbin,
        pn.pane.Markdown(
            "The Negative Binomial regression models above show a significant relationship between clonal and sub-clonal counts. "
            "The **likelihood ratio test** shows that the Negative Binomial model is a good fit for the data."
            "Not sure if it's justifiable to compare the likelihood values between Poisson and Negative Binomial models, but the Negative Binomial model has indead a likelihood value closer to 1.",
            styles=comment_style,
            width=300 # Adjust width to fit text
            ),
        ),
    pn.pane.Markdown("## Negative Binomial Regression Hexplots Comparison", styles=header_style),
    pn.Row(hex_fig_left_negbin, hex_fig_right_negbin),
    sizing_mode='stretch_width'
)

dashboard.servable()