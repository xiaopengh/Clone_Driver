import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import statsmodels.api as sm
# from scipy import stats # No longer needed for placeholder binomial test
import pickle as pkl # Added for loading data

# --- Configuration (Adapted from ui_copy.py) ---
# Streamlit handles layout and styling differently, so direct mapping isn't always needed.
# We'll use Streamlit's built-in styling and layout features.
CONFIG = {
    'ndecimal': 2,
}

# --- Static Content (Adapted from ui_copy.py) ---
CONTENT = {
    'title': "🧬 Tumor Clonality Dashboard",
    'headers': {
        'left_overview': "Binomial Test Generated Dataset Overview (Grouped by Gene)", # Updated header
        'right_overview': "Z Test Generated Dataset Overview (Grouped by Gene)", # Updated header
        'left_stats': "Binomial Test Generated Data Summary Statistics",
        'right_stats': "Z Test Generated Data Summary Statistics",
        'clonality_counts': "Clonality Counts (Across All Genes)", # Updated header
        'poisson_summary': "Poisson Regression Model Summaries (Gene Level)", # Updated header
        'poisson_plots': "Poisson Regression Plots Comparison (Gene Level)", # Updated header
        'negbin_summary': "Negative Binomial Regression Model Summaries (Gene Level)", # Updated header
        'negbin_plots': "Negative Binomial Regression Plots Comparison (Gene Level)" # Updated header
    },
    'comments': {
        'stats': "Both datasets show potential **overdispersion** (variance > mean) in clonal/sub-clonal counts per gene, suggesting Negative Binomial regression might be more appropriate than Poisson.", # Updated comment
        'clonality': "Comparison of total clonal/sub-clonal counts across all genes based on Binomial vs. Z-test methods.", # Updated comment
        'poisson': "Poisson regression results modeling sub-clonal counts based on clonal counts per gene. Note potential issues if data is overdispersed.", # Updated comment
        'negbin': "Negative Binomial regression results, often better suited for overdispersed count data like clonal counts per gene." # Updated comment
    }
}

# --- Data Loading and Analysis Functions (Integrated from server_copy logic) ---

@st.cache_data # Cache the results of this function
def load_and_process_data():
    """Loads and processes both datasets."""
    try:
        # Load
        data_binom = pkl.load(open("data/mc3_binom_left.pkl", "rb")).reset_index(drop=True)
        data_binom["Clonality"] = data_binom["Clonality"].astype(str)

        data_ztest = pkl.load(open("data/mc3_Ztest_left.pkl", "rb")).reset_index(drop=True)
        data_ztest["Clonality"] = data_ztest["Clonality"].astype(str)

        # Process (Group by Gene)
        summary_binom = data_binom.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
        summary_binom.columns.name = None
        summary_binom = summary_binom.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"}).reset_index() # Keep Hugo_Symbol

        summary_ztest = data_ztest.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
        summary_ztest.columns.name = None
        summary_ztest = summary_ztest.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"}).reset_index() # Keep Hugo_Symbol

        # Ensure both count columns exist, filling with 0 if one clonality type is missing for a gene
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


@st.cache_data
def calculate_summary_stats(clonal_summary):
    """Calculates summary statistics from the gene summary dataframe."""
    if clonal_summary is None or clonal_summary.empty:
        return pd.DataFrame() # Return empty dataframe if input is invalid

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
def fit_poisson_model(summary):
    """Fits a Poisson regression model using gene summary data and returns summary text and predictions."""
    if summary is None or summary.empty or "Clonal Count" not in summary.columns or "Sub-Clonal Count" not in summary.columns:
        return "Error: Insufficient data for Poisson model.", None
    try:
        # Ensure columns exist and handle potential missing values if necessary
        summary = summary.fillna(0) # Or handle appropriately
        x = summary["Clonal Count"]
        y = summary["Sub-Clonal Count"]
        x_const = sm.add_constant(x)
        model = sm.GLM(y, x_const, family=sm.families.Poisson()).fit()
        predictions = model.predict(x_const)
        # Sort predictions based on x for a smooth line plot
        sorted_indices = x.argsort()
        sorted_x = x.iloc[sorted_indices]
        sorted_predictions = predictions.iloc[sorted_indices]
        
        # Create a dataframe for easier plotting with plotly trace
        pred_df = pd.DataFrame({'Clonal Count': sorted_x, 'Predicted Sub-Clonal Count': sorted_predictions})

        return model.summary().as_text(), pred_df
    except Exception as e:
        return f"Error fitting Poisson model: {e}", None

@st.cache_data
def fit_negbin_model(summary):
    """Fits a Negative Binomial regression model using gene summary data."""
    if summary is None or summary.empty or "Clonal Count" not in summary.columns or "Sub-Clonal Count" not in summary.columns:
        return "Error: Insufficient data for Negative Binomial model."
    try:
        # Ensure columns exist and handle potential missing values if necessary
        summary = summary.fillna(0) # Or handle appropriately
        x = summary["Clonal Count"]
        y = summary["Sub-Clonal Count"]
        x_const = sm.add_constant(x)
        # Use GLM for Negative Binomial for consistency and easier summary access
        # Adjust alpha as needed, starting with 1.0 is common.
        model = sm.GLM(y, x_const, family=sm.families.NegativeBinomial(alpha=1.0)).fit()
        # Alternative: model = sm.NegativeBinomial(y, x_const).fit() # May require different summary handling
        return model.summary().as_text()
    except Exception as e:
        return f"Error fitting Negative Binomial model: {e}"

# --- Main App Logic ---
st.set_page_config(layout="wide")
st.title(CONTENT['title'])

# --- Sidebar for Controls ---
st.sidebar.header("Controls")
# Add button to clear cache
if st.sidebar.button("Clear Cache and Rerun"):
    st.cache_data.clear()
    st.rerun()

# --- Data Loading and Processing ---
# Load and process data using the new functions
data_binom_raw, data_ztest_raw, summary_binom, summary_ztest = load_and_process_data()

# Proceed only if data loading was successful
if data_binom_raw is not None and data_ztest_raw is not None and summary_binom is not None and summary_ztest is not None:

    # Calculate stats using the gene summaries
    stats_binom = calculate_summary_stats(summary_binom)
    stats_ztest = calculate_summary_stats(summary_ztest)

    # Fit models using the gene summaries - unpack predictions
    poisson_summary_binom_text, poisson_pred_binom = fit_poisson_model(summary_binom)
    poisson_summary_ztest_text, poisson_pred_ztest = fit_poisson_model(summary_ztest)
    negbin_summary_binom = fit_negbin_model(summary_binom) # Assuming NegBin doesn't need predictions plotted for now
    negbin_summary_ztest = fit_negbin_model(summary_ztest) # Assuming NegBin doesn't need predictions plotted for now

    # --- Dashboard Layout ---

    # Section 1: Data Overview (Show Gene Summaries)
    st.header("Data Overview (Gene Summaries)")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader(CONTENT['headers']['left_overview'])
        st.dataframe(summary_binom.head(100), height=300) # Display first 100 rows of summary
    with col2:
        st.subheader(CONTENT['headers']['right_overview'])
        st.dataframe(summary_ztest.head(100), height=300) # Display first 100 rows of summary

    # Section 2: Summary Statistics (Based on Gene Summaries)
    st.header("Summary Statistics (Based on Gene Counts)")
    col1, col2, col3 = st.columns([2, 2, 1]) # Adjust column ratios as needed
    with col1:
        st.subheader(CONTENT['headers']['left_stats'])
        st.dataframe(stats_binom)
    with col2:
        st.subheader(CONTENT['headers']['right_stats'])
        st.dataframe(stats_ztest)
    with col3:
        st.info(CONTENT['comments']['stats'])

    # Section 3: Clonality Counts (Total counts from raw data)
    st.header(CONTENT['headers']['clonality_counts'])
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        # Calculate counts from the raw loaded data
        counts_binom = data_binom_raw['Clonality'].value_counts().reset_index()
        counts_binom.columns = ['Clonality', 'count']
        fig_bar_binom = px.bar(counts_binom, x='Clonality', y='count', color='Clonality',
                               title="Binomial Test Total Clonality Counts")
        st.plotly_chart(fig_bar_binom, use_container_width=True)
    with col2:
        # Calculate counts from the raw loaded data
        counts_ztest = data_ztest_raw['Clonality'].value_counts().reset_index()
        counts_ztest.columns = ['Clonality', 'count']
        fig_bar_ztest = px.bar(counts_ztest, x='Clonality', y='count', color='Clonality',
                               title="Z-Test Total Clonality Counts")
        st.plotly_chart(fig_bar_ztest, use_container_width=True)
    with col3:
        st.info(CONTENT['comments']['clonality'])

    # Section 4: Poisson Regression
    st.header(CONTENT['headers']['poisson_summary'])
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        st.subheader("Binomial Data Model")
        st.code(poisson_summary_binom_text) # Use the text summary
    with col2:
        st.subheader("Z-Test Data Model")
        st.code(poisson_summary_ztest_text) # Use the text summary
    with col3:
        st.info(CONTENT['comments']['poisson'])

    # Section 5: Poisson Plots (Using Scatter: Sub-Clonal vs Clonal Counts per Gene)
    st.header(CONTENT['headers']['poisson_plots'])
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Binomial Data: Sub-Clonal vs Clonal Counts (Gene Level)")

        # --- Calculate Overlap Count ---
        summary_binom_counts = summary_binom.groupby(['Clonal Count', 'Sub-Clonal Count']).size().reset_index(name='Overlap Count')
        summary_binom_plot_data = pd.merge(summary_binom, summary_binom_counts, on=['Clonal Count', 'Sub-Clonal Count'])
        summary_binom_plot_data['Log Overlap Count'] = np.log(summary_binom_plot_data['Overlap Count'] + 1) # Log transform for better visualization
        # --- End Calculation ---

        # Base scatter plot - Modified for density and square markers
        fig_poisson_binom = px.scatter(summary_binom_plot_data, # Use the data with overlap counts
                                        x='Clonal Count',
                                        y='Sub-Clonal Count',
                                        title="Poisson Data & Predictions (Binomial - Gene Level)",
                                        hover_data=['Hugo_Symbol', 'Clonal Count', 'Sub-Clonal Count', 'Overlap Count'], # Add Overlap Count to hover
                                        color='Log Overlap Count', # Color by the calculated overlap count
                                        color_continuous_scale=px.colors.sequential.Viridis, # Choose a color scale
                                        symbol_sequence=['square']) # Use square markers
                                        # trendline="ols", trendline_color_override="rgba(255,0,0,0.5)") # Trendline removed from base scatter
        
        # --- Add OLS Trendline Separately ---
        # Fit OLS manually to get coefficients
        x_ols = sm.add_constant(summary_binom['Clonal Count'])
        y_ols = summary_binom['Sub-Clonal Count']
        ols_model = sm.OLS(y_ols, x_ols).fit()
        ols_pred = ols_model.predict(x_ols)
        # Create a dataframe sorted by x for the line
        ols_line_df = pd.DataFrame({'Clonal Count': summary_binom['Clonal Count'], 'OLS Fit': ols_pred}).sort_values('Clonal Count')
        # Add OLS trace
        fig_poisson_binom.add_trace(go.Scatter(x=ols_line_df['Clonal Count'],
                                                y=ols_line_df['OLS Fit'],
                                                mode='lines',
                                                name='OLS Trendline',
                                                line=dict(color='rgba(255,0,0,0.5)', width=2))) # Semi-transparent red
        # --- End OLS Trendline ---


        # Add Poisson prediction line if available
        if poisson_pred_binom is not None:
            fig_poisson_binom.add_trace(go.Scatter(x=poisson_pred_binom['Clonal Count'],
                                y=poisson_pred_binom['Predicted Sub-Clonal Count'],
                                mode='lines',
                                name='Poisson Prediction',
                                line=dict(color='orange', dash='dash', width=5)), # Dashed orange line with increased width
                            secondary_y=False) # Ensure it's on the same y-axis


        st.plotly_chart(fig_poisson_binom, use_container_width=True)
    with col2:
        st.subheader("Z-Test Data: Sub-Clonal vs Clonal Counts (Gene Level)")

        # --- Calculate Overlap Count for Z-Test Data ---
        summary_ztest_counts = summary_ztest.groupby(['Clonal Count', 'Sub-Clonal Count']).size().reset_index(name='Overlap Count')
        summary_ztest_plot_data = pd.merge(summary_ztest, summary_ztest_counts, on=['Clonal Count', 'Sub-Clonal Count'])
        # --- End Calculation ---

        # Base scatter plot - Modified for density and square markers
        fig_poisson_ztest = px.scatter(summary_ztest_plot_data, # Use the data with overlap counts
                                        x='Clonal Count',
                                        y='Sub-Clonal Count',
                                        title="Poisson Data & Predictions (Z-Test - Gene Level)",
                                        hover_data=['Hugo_Symbol', 'Clonal Count', 'Sub-Clonal Count', 'Overlap Count'], # Add Overlap Count to hover
                                        color='Overlap Count', # Color by the calculated overlap count
                                        color_continuous_scale=px.colors.sequential.Viridis, # Choose a color scale
                                        symbol_sequence=['square']) # Use square markers
                                        # trendline="ols", trendline_color_override="rgba(255,0,0,0.5)") # Trendline removed from base scatter
        

        # --- Add OLS Trendline Separately for Z-Test ---
        x_ols_z = sm.add_constant(summary_ztest['Clonal Count'])
        y_ols_z = summary_ztest['Sub-Clonal Count']
        ols_model_z = sm.OLS(y_ols_z, x_ols_z).fit()
        ols_pred_z = ols_model_z.predict(x_ols_z)
        ols_line_df_z = pd.DataFrame({'Clonal Count': summary_ztest['Clonal Count'], 'OLS Fit': ols_pred_z}).sort_values('Clonal Count')
        fig_poisson_ztest.add_trace(go.Scatter(x=ols_line_df_z['Clonal Count'],
                                                y=ols_line_df_z['OLS Fit'],
                                                mode='lines',
                                                name='OLS Trendline',
                                                line=dict(color='rgba(255,0,0,0.5)', width=2)))
        # --- End OLS Trendline ---

        # Add Poisson prediction line if available
        if poisson_pred_ztest is not None:
            fig_poisson_ztest.add_trace(go.Scatter(x=poisson_pred_ztest['Clonal Count'],
                                                    y=poisson_pred_ztest['Predicted Sub-Clonal Count'],
                                                    mode='lines',
                                                    name='Poisson Prediction',
                                                    line=dict(color='orange', dash='dash', width=3))) # Dashed orange line

        st.plotly_chart(fig_poisson_ztest, use_container_width=True)    
    

    # Section 6: Negative Binomial Regression
    st.header(CONTENT['headers']['negbin_summary'])
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        st.subheader("Binomial Data Model")
        st.code(negbin_summary_binom)
    with col2:
        st.subheader("Z-Test Data Model")
        st.code(negbin_summary_ztest)
    with col3:
        st.info(CONTENT['comments']['negbin'])

    # Section 7: Negative Binomial Plots (Using Scatter: Sub-Clonal vs Clonal Counts per Gene)
    st.header(CONTENT['headers']['negbin_plots'])
    col1, col2 = st.columns(2)
    # Note: These plots show the same underlying gene-level count data as the Poisson section.
    # Trendlines here might differ slightly if using a NegBin specific fit, but OLS is used for simplicity.
    with col1:
        st.subheader("Binomial Data: Sub-Clonal vs Clonal Counts (Gene Level)")
        fig_negbin_binom = px.scatter(summary_binom, x='Clonal Count', y='Sub-Clonal Count',
                                      title="NegBin Data (Binomial - Gene Level)",
                                      hover_data=['Hugo_Symbol', 'Clonal Count', 'Sub-Clonal Count'],
                                      trendline="ols", trendline_color_override="blue") # Add trendline
        st.plotly_chart(fig_negbin_binom, use_container_width=True)
    with col2:
        st.subheader("Z-Test Data: Sub-Clonal vs Clonal Counts (Gene Level)")
        fig_negbin_ztest = px.scatter(summary_ztest, x='Clonal Count', y='Sub-Clonal Count',
                                      title="NegBin Data (Z-Test - Gene Level)",
                                      hover_data=['Hugo_Symbol', 'Clonal Count', 'Sub-Clonal Count'],
                                      trendline="ols", trendline_color_override="blue") # Add trendline
        st.plotly_chart(fig_negbin_ztest, use_container_width=True)

    # --- Optional: Display Raw Loaded Data ---
    with st.expander("Show Raw Loaded Data (Binomial)"):
        st.dataframe(data_binom_raw)
    with st.expander("Show Raw Loaded Data (Z-Test)"):
        st.dataframe(data_ztest_raw)

# --- How to Run ---
st.sidebar.markdown("---")
st.sidebar.info("To run this app, save the code as `streamlit_dashboard.py`, ensure your .pkl files are in a 'data' subfolder, and run `streamlit run streamlit_dashboard.py` in your terminal.")

# !!! Reminder: Ensure 'data/mc3_binom_left.pkl' and 'data/mc3_Ztest_left.pkl' exist relative to where you run the script !!!
