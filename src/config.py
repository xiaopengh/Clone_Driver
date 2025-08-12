# Configuration settings
CONFIG = {
    'ndecimal': 2,
}

# Static content for the dashboard
CONTENT = {
    'title': "🧬 Tumor Clonality Dashboard",
    'headers': {
        'left_overview': "Binomial Test Generated Dataset Overview (Grouped by Gene)", 
        'right_overview': "Z Test Generated Dataset Overview (Grouped by Gene)", 
        'left_stats': "Binomial Test Generated Data Summary Statistics",
        'right_stats': "Z Test Generated Data Summary Statistics",
        'clonality_counts': "Clonality Counts (Across All Genes)", 
        'poisson_summary': "Poisson Regression Model Summaries ", 
        'poisson_plots': "Poisson Regression Plots Comparison ",
        'left_poisson_plot': "Binomial Data: Sub-Clonal vs Clonal Counts ",
        'right_poisson_plot': "Z-Test Data: Sub-Clonal vs Clonal Counts ",
        'negbin_summary': "Negative Binomial Regression Model Summaries ", 
        'negbin_plots': "Negative Binomial Regression Plots Comparison ",
        'left_negbin_plot': "Binomial Data: Sub-Clonal vs Clonal Counts ",
        'right_negbin_plot': "Z-Test Data: Sub-Clonal vs Clonal Counts "
    },
    'comments': {
        'stats': "Both datasets show potential **overdispersion** (variance > mean) in clonal/sub-clonal counts per gene, suggesting Negative Binomial regression might be more appropriate than Poisson.",
        'clonality': "Comparison of total clonal/sub-clonal counts across all genes based on Binomial vs. Z-test methods.",
        'poisson': "Poisson regression results modeling sub-clonal counts based on clonal counts per gene. Note potential issues if data is overdispersed.",
        'negbin': "Negative Binomial regression results, often better suited for overdispersed count data like clonal counts per gene."
    }
} 