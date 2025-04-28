import streamlit as st
import pandas as pd
import statsmodels.api as sm

@st.cache_data
def fit_poisson_model(summary):
    """Fits a Poisson regression model using gene summary data and returns summary text and predictions."""
    if summary is None or summary.empty or "Clonal Count" not in summary.columns or "Sub-Clonal Count" not in summary.columns:
        return "Error: Insufficient data for Poisson model.", None
    try:
        summary = summary.fillna(0)
        x = summary["Clonal Count"]
        y = summary["Sub-Clonal Count"]
        x_const = sm.add_constant(x)
        model = sm.GLM(y, x_const, family=sm.families.Poisson()).fit()
        predictions = model.predict(x_const)
        
        # Sort predictions based on x for a smooth line plot
        sorted_indices = x.argsort()
        sorted_x = x.iloc[sorted_indices]
        sorted_predictions = predictions.iloc[sorted_indices]
        
        pred_df = pd.DataFrame({'Clonal Count': sorted_x, 'Predicted Sub-Clonal Count': sorted_predictions})
        return model.summary().as_text(), pred_df
    except Exception as e:
        return f"Error fitting Poisson model: {e}", None

@st.cache_data
def fit_negbin_model(summary):
    """Fits a Negative Binomial regression model using gene summary data."""
    if summary is None or summary.empty or "Clonal Count" not in summary.columns or "Sub-Clonal Count" not in summary.columns:
        return "Error: Insufficient data for Negative Binomial model.", None
    try:
        summary = summary.fillna(0)
        x = summary["Clonal Count"]
        y = summary["Sub-Clonal Count"]
        x_const = sm.add_constant(x)
        model = sm.NegativeBinomial(y, x_const).fit(disp=False)
        predictions = model.predict(x_const)

        # Sort predictions based on x for a smooth line plot
        sorted_indices = x.argsort()
        sorted_x = x.iloc[sorted_indices]
        sorted_predictions = predictions.iloc[sorted_indices]

        pred_df = pd.DataFrame({'Clonal Count': sorted_x, 'Predicted Sub-Clonal Count': sorted_predictions})

        # Extract estimated alpha
        alpha_est = model.params_alpha if hasattr(model, "params_alpha") else model.scale
        summary_text = model.summary().as_text()
        summary_text += f"\n\nEstimated alpha (dispersion parameter): {alpha_est:.4f}"

        return summary_text, pred_df
    except Exception as e:
        return f"Error fitting Negative Binomial model: {e}", None 