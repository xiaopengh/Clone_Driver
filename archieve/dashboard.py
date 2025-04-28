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

# === Configuration Section ===
# Define all styling and layout constants here
CONFIG = {
    'figsize': (5, 4),
    'ndecimal': 2,
    'header_style': {'background': '#F0F8FF', 'padding': '10px', 'border-radius': '5px'},
    'comment_style': {'background': '#e6ffe6', 'padding': '10px', 'border-radius': '5px', 'font-size': '16px'},
    'title_style': {'color':'#2F80ED'},
    'table_height': 300,
    'comment_width': 300
}

# === Static Content ===
CONTENT = {
    'title': "# 🧬 Tumor Clonality Dashboard",
    'headers': {
        'left_overview': "### Left Dataset Overview",
        'right_overview': "### Right Dataset Overview",
        'left_stats': "### Left Data Summary Statistics",
        'right_stats': "### Right Data Summary Statistics",
        'clonality_counts': "## Clonality Counts",
        'poisson_summary': "## Poisson Regression Model Summaries",
        'poisson_hexplots': "## Poisson Regression Hexplots Comparison",
        'negbin_summary': "## Negative Binomial Regression Model Summaries",
        'negbin_hexplots': "## Negative Binomial Regression Hexplots Comparison"
    },
    'comments': {
        'stats': "Both datasets have clear evidence of **overdispersion**, as the variance is much larger than the mean. "
                "This is expected in count data, and is why we are using **Poisson** and **Negative Binomial** regression models.",
        'clonality': "The bar plots on the left show the counts of clonal and sub-clonal mutations in the datasets. "
                    "The left dataset has more clonal mutations, while the right dataset has more sub-clonal mutations due to different settings of $H_0$.",
        'poisson': "The Poisson regression models above show a significant relationship between clonal and sub-clonal counts. "
                  "Even though the data is overdispersed, the **Pseudo-R-square** values are misleadingly high, **WHY**?",
        'negbin': "The Negative Binomial regression models above show a significant relationship between clonal and sub-clonal counts. "
                 "The **likelihood ratio test** shows that the Negative Binomial model is a good fit for the data. "
                 "Not sure if it's justifiable to compare the likelihood values between Poisson and Negative Binomial models, but the Negative Binomial model has indead a likelihood value closer to 1."
    }
}

# Helper function to create a header
def create_header(key):
    return pn.pane.Markdown(CONTENT['headers'][key], styles=CONFIG['header_style'])

# Helper function to create a comment
def create_comment(key):
    return pn.pane.Markdown(CONTENT['comments'][key], styles=CONFIG['comment_style'], width=CONFIG['comment_width'])

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
    fig, ax = plt.subplots(figsize=CONFIG['figsize'])
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
    fig, ax = plt.subplots(figsize=CONFIG['figsize'])
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
    y_pred = poisson_model.predict(sm.add_constant(x_plot))

    fig, ax = plt.subplots(figsize=CONFIG['figsize'])
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
data_table_left = pn.widgets.Tabulator(data_left.head(int(len(data_left)*0.3)), pagination='remote', page_size=10, height=CONFIG['table_height'])
data_table_right = pn.widgets.Tabulator(data_right.head(int(len(data_right)*0.3)), pagination='remote', page_size=10, height=CONFIG['table_height'])

summary_stats_left = pd.DataFrame({
    "Statistic": ["Mean", 'Median', "Variance"],
    "Clonal": [
        clonal_summary_left["Clonal Count"].mean().round(CONFIG['ndecimal']),
        clonal_summary_left["Clonal Count"].median().round(CONFIG['ndecimal']),
        clonal_summary_left["Clonal Count"].var().round(CONFIG['ndecimal'])
    ],
    "Sub-Clonal": [
        clonal_summary_left["Sub-Clonal Count"].mean().round(CONFIG['ndecimal']),
        clonal_summary_left["Sub-Clonal Count"].median().round(CONFIG['ndecimal']),
        clonal_summary_left["Sub-Clonal Count"].var().round(CONFIG['ndecimal'])
    ]
})

summary_stats_right = pd.DataFrame({
    "Statistic": ["Mean", 'Median', "Variance"],
    "Clonal": [
        clonal_summary_right["Clonal Count"].mean().round(CONFIG['ndecimal']),
        clonal_summary_right["Clonal Count"].median().round(CONFIG['ndecimal']),
        clonal_summary_right["Clonal Count"].var().round(CONFIG['ndecimal'])
    ],
    "Sub-Clonal": [
        clonal_summary_right["Sub-Clonal Count"].mean().round(CONFIG['ndecimal']),
        clonal_summary_right["Sub-Clonal Count"].median().round(CONFIG['ndecimal']),
        clonal_summary_right["Sub-Clonal Count"].var().round(CONFIG['ndecimal'])
    ]
})

stats_table_left = pn.widgets.Tabulator(summary_stats_left, pagination='remote', page_size=10, height=CONFIG['table_height'])
stats_table_right = pn.widgets.Tabulator(summary_stats_right, pagination='remote', page_size=10, height=CONFIG['table_height'])

# === Dashboard Layout using modular components ===
title = pn.pane.Markdown(CONTENT['title'], margin=(0,0,20,0), styles=CONFIG['title_style'])

# Create sections using helper functions
data_overview_section = pn.Row(
    pn.Column(create_header('left_overview'), data_table_left),
    pn.Column(create_header('right_overview'), data_table_right)
)

stats_section = pn.Row(
    pn.Column(create_header('left_stats'), stats_table_left),
    pn.Column(create_header('right_stats'), stats_table_right),
    create_comment('stats')
)

clonality_section = pn.Column(
    create_header('clonality_counts'),
    pn.Row(bar_fig_left, bar_fig_right, create_comment('clonality'))
)

poisson_summary_section = pn.Column(
    create_header('poisson_summary'),
    pn.Row(summary_left, summary_right, create_comment('poisson'))
)

poisson_hexplots_section = pn.Column(
    create_header('poisson_hexplots'),
    pn.Row(hex_fig_left, hex_fig_right)
)

negbin_summary_section = pn.Column(
    create_header('negbin_summary'),
    pn.Row(summary_left_negbin, summary_right_negbin, create_comment('negbin'))
)

negbin_hexplots_section = pn.Column(
    create_header('negbin_hexplots'),
    pn.Row(hex_fig_left_negbin, hex_fig_right_negbin)
)

# Assemble the dashboard
dashboard = pn.Column(
    title,
    data_overview_section,
    stats_section,
    clonality_section,
    poisson_summary_section,
    poisson_hexplots_section,
    negbin_summary_section,
    negbin_hexplots_section,
    sizing_mode='stretch_width'
)

dashboard.servable()