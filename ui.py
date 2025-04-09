import panel as pn
from server import get_dashboard_data

pn.extension('tabulator', sizing_mode='stretch_width')

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
        'left_overview': "### Binomial Test Generated Dataset Overview",
        'right_overview': "### Z Test Generated Dataset Overview",
        'left_stats': "### Binomial Test Generated Data Summary Statistics",
        'right_stats': "### Z Test Generated Data Summary Statistics",
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
                    "The Binomial test generated dataset has more clonal mutations, while the Z Test generated dataset has more sub-clonal mutations due to different settings of $H_0$.",
        'poisson': "The Poisson regression models above show a significant relationship between clonal and sub-clonal counts. "
                  "Even though the data is overdispersed, the **Pseudo-R-square** values are misleadingly high, **WHY**?",
        'negbin': "The Negative Binomial regression models above show a significant relationship between clonal and sub-clonal counts. "
                 "The **likelihood ratio test** shows that the Negative Binomial model is a good fit for the data. "
                 "Not sure if it's justifiable to compare the likelihood values between Poisson and Negative Binomial models, but the Negative Binomial model has indead a likelihood value closer to 1."
    }
}

# Helper functions for UI components
def create_header(key):
    return pn.pane.Markdown(CONTENT['headers'][key], styles=CONFIG['header_style'])

def create_comment(key):
    return pn.pane.Markdown(CONTENT['comments'][key], styles=CONFIG['comment_style'], width=CONFIG['comment_width'])

def create_dashboard(force_recalculate=False):
    # Get all data from server, potentially using the cache
    data = get_dashboard_data(force_recalculate)
    
    # Create tables using pre-calculated data
    data_table_left = pn.widgets.Tabulator(
        data['table_data_left'], 
        pagination='remote', 
        page_size=10, 
        height=CONFIG['table_height']
    )
    
    data_table_right = pn.widgets.Tabulator(
        data['table_data_right'], 
        pagination='remote', 
        page_size=10, 
        height=CONFIG['table_height']
    )
    
    stats_table_left = pn.widgets.Tabulator(
        data['summary_stats_left'], 
        pagination='remote', 
        page_size=10, 
        height=CONFIG['table_height']
    )
    
    stats_table_right = pn.widgets.Tabulator(
        data['summary_stats_right'], 
        pagination='remote', 
        page_size=10, 
        height=CONFIG['table_height']
    )
    
    # Create visualizations using pre-calculated figures
    bar_fig_left = pn.pane.Matplotlib(data['visualizations']['bar_left'], tight=True)
    bar_fig_right = pn.pane.Matplotlib(data['visualizations']['bar_right'], tight=True)
    hex_fig_left = pn.pane.Matplotlib(data['visualizations']['hex_left'])
    hex_fig_right = pn.pane.Matplotlib(data['visualizations']['hex_right'])
    hex_fig_left_negbin = pn.pane.Matplotlib(data['visualizations']['hex_left_negbin'])
    hex_fig_right_negbin = pn.pane.Matplotlib(data['visualizations']['hex_right_negbin'])
    
    # Model summaries using pre-calculated text
    summary_left = pn.pane.Markdown(data['model_summaries']['poisson_left'])
    summary_right = pn.pane.Markdown(data['model_summaries']['poisson_right'])
    summary_left_negbin = pn.pane.Markdown(data['model_summaries']['negbin_left'])
    summary_right_negbin = pn.pane.Markdown(data['model_summaries']['negbin_right'])
    
    # === Create Dashboard Layout ===
    title = pn.pane.Markdown(CONTENT['title'], margin=(0,0,20,0), styles=CONFIG['title_style'])
    
    # Create sections
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
    
    return dashboard

# Create buttons
refresh_button = pn.widgets.Button(name="Refresh Data", button_type="primary")
clear_cache_button = pn.widgets.Button(name="Clear Cache", button_type="warning")

# Function to refresh dashboard with recalculated data
def refresh_dashboard(event):
    # Replace the current dashboard with a new one using recalculated data
    new_dashboard = create_dashboard(force_recalculate=True)
    dashboard_container.objects = [pn.Row(refresh_button, clear_cache_button), new_dashboard]

# Function to clear the cache without refreshing dashboard
def clear_cache(event):
    # This assumes there's a clear_cache function in server module
    # If not, you'll need to implement it based on your caching mechanism
    from server import clear_cache
    clear_cache()
    # Provide feedback that cache was cleared
    clear_cache_button.name = "Cache Cleared!"
    # Reset button text after 2 seconds
    pn.state.add_timeout(2000, lambda: setattr(clear_cache_button, 'name', 'Clear Cache'))

# Attach the functions to the buttons
refresh_button.on_click(refresh_dashboard)
clear_cache_button.on_click(clear_cache)

# Create dashboard with cached data (if available)
dashboard = create_dashboard()

# Container to hold the buttons and dashboard
dashboard_container = pn.Column(
    pn.Row(refresh_button, clear_cache_button),
    dashboard,
    sizing_mode='stretch_width'
)

# Make the container servable
dashboard_container.servable()