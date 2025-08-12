import streamlit as st
from src.config import CONTENT
from src.data_handlers import load_and_process_data, extend_clonal_summary, calculate_summary_stats
from src.models import fit_poisson_model, fit_negbin_model
from src.visualizations import plot_scatter, plot_clonality_counts

def summary_visualization(summary_binom, summary_ztest, suffix=''):
    """
    Display summary visualizations with unique keys for Plotly charts
    Args:
        summary_binom: DataFrame for binomial data
        summary_ztest: DataFrame for z-test data
        suffix: String to make chart keys unique
    """
    # Calculate stats using the gene summaries
    stats_binom = calculate_summary_stats(summary_binom)
    stats_ztest = calculate_summary_stats(summary_ztest)

    # Fit models using the gene summaries
    poisson_summary_binom_text, poisson_pred_binom = fit_poisson_model(summary_binom)
    poisson_summary_ztest_text, poisson_pred_ztest = fit_poisson_model(summary_ztest)
    negbin_summary_binom, negbin_pred_binom, scale_binom = fit_negbin_model(summary_binom)
    negbin_summary_ztest, negbin_pred_ztest, scale_ztest = fit_negbin_model(summary_ztest)

    # --- Dashboard Layout ---

    # Section 1: Data Overview
    st.header("Data Overview (Gene Summaries)")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader(CONTENT['headers']['left_overview'])
        st.dataframe(summary_binom, height=300)
    with col2:
        st.subheader(CONTENT['headers']['right_overview'])
        st.dataframe(summary_ztest, height=300)
    
    # Section 2: Summary Statistics
    st.header("Summary Statistics (Based on Gene Counts)")
    col1, col2, col3 = st.columns([3, 3, 1])
    with col1:
        st.subheader(CONTENT['headers']['left_stats'])
        st.dataframe(stats_binom)
    with col2:
        st.subheader(CONTENT['headers']['right_stats'])
        st.dataframe(stats_ztest)
    with col3:
        st.info(CONTENT['comments']['stats'])

    # Section 3: Clonality Counts
    st.header(CONTENT['headers']['clonality_counts'])
    col1, col2 = st.columns([1, 1])
    with col1:
        fig_bar_binom = plot_clonality_counts(data_binom_raw, "Binomial Test Total Clonality Counts")
        st.plotly_chart(fig_bar_binom, use_container_width=True, key=f"bar_binom{suffix}")
    with col2:
        fig_bar_ztest = plot_clonality_counts(data_ztest_raw, "Z-Test Total Clonality Counts")
        st.plotly_chart(fig_bar_ztest, use_container_width=True, key=f"bar_ztest{suffix}")

    # Section 4: Poisson Regression
    st.header(CONTENT['headers']['poisson_summary'])
    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("Binomial Data Model")
        st.code(poisson_summary_binom_text)
    with col2:
        st.subheader("Z-Test Data Model")
        st.code(poisson_summary_ztest_text)

    # Section 5: Poisson Plots
    st.header(CONTENT['headers']['poisson_plots'])
    col1, col2 = st.columns(2)
    with col1:
        st.subheader(CONTENT['headers']['left_poisson_plot'])
        fig_poisson_binom = plot_scatter(summary_binom, poisson_pred_binom)
        st.plotly_chart(fig_poisson_binom, use_container_width=True, key=f"poisson_binom{suffix}")
    with col2:
        st.subheader(CONTENT['headers']['right_poisson_plot'])
        fig_poisson_ztest = plot_scatter(summary_ztest, poisson_pred_ztest)
        st.plotly_chart(fig_poisson_ztest, use_container_width=True, key=f"poisson_ztest{suffix}")

    # Section 6: Negative Binomial Regression
    st.header(CONTENT['headers']['negbin_summary'])
    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("Binomial Data Model")
        st.code(negbin_summary_binom)
    with col2:
        st.subheader("Z-Test Data Model")
        st.code(negbin_summary_ztest)

    # Section 7: Negative Binomial Plots
    st.header(CONTENT['headers']['negbin_plots'])
    col1, col2 = st.columns(2)
    with col1:
        st.subheader(CONTENT['headers']['left_negbin_plot'])
        fig_negbin_binom = plot_scatter(summary_binom, negbin_pred_binom, scale_binom)
        st.plotly_chart(fig_negbin_binom, use_container_width=True, key=f"negbin_binom{suffix}")
    with col2:
        st.subheader(CONTENT['headers']['right_negbin_plot'])
        fig_negbin_ztest = plot_scatter(summary_ztest, negbin_pred_ztest, scale_ztest)
        st.plotly_chart(fig_negbin_ztest, use_container_width=True, key=f"negbin_ztest{suffix}")

# --- Main App Logic ---
st.set_page_config(layout="wide")
st.title(CONTENT['title'])

# --- Sidebar for Controls ---
st.sidebar.header("Controls")
if st.sidebar.button("Clear Cache and Rerun"):
    st.cache_data.clear()
    st.rerun()

# --- Data Loading and Processing ---
data_binom_raw, data_ztest_raw, summary_binom, summary_ztest = load_and_process_data()
summary_binom_extended = extend_clonal_summary(summary_binom)
summary_ztest_extended = extend_clonal_summary(summary_ztest)

# Proceed only if data loading was successful
if data_binom_raw is not None and data_ztest_raw is not None and summary_binom is not None and summary_ztest is not None:
    
    st.header("This doesn't include genes with zero mutations")
    summary_visualization(summary_binom, summary_ztest)

    st.header("This includes data with 0 clonal count and 0 subclonal count")
    summary_visualization(summary_binom_extended, summary_ztest_extended, "_extended")

    # # --- Optional: Display Raw Loaded Data ---
    # with st.expander("Show Raw Loaded Data (Binomial)"):
    #     st.dataframe(data_binom_raw)
    # with st.expander("Show Raw Loaded Data (Z-Test)"):
    #     st.dataframe(data_ztest_raw)

st.sidebar.markdown("---")