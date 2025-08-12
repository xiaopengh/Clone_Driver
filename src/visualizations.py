import plotly.express as px
import plotly.graph_objects as go
import statsmodels.api as sm
import pandas as pd
from src.data_handlers import add_overlap_counts, negbin_quantiles

def plot_scatter(summary_df, pred_df, scale=None, plot_title=None):
    """
    Plots a scatter plot of Sub-Clonal vs Clonal Counts with overlap coloring,
    OLS trendline, and Poisson prediction line.
    """
    # Add overlap counts to the summary dataframe
    plot_df = add_overlap_counts(summary_df)

    # Create a new column for the log of the overlap count
    fig = go.Figure()

    # Add the scatter points
    scatter = go.Scatter(
        x=plot_df['Clonal Count'],
        y=plot_df['Sub-Clonal Count'],
        mode='markers',
        marker=dict(
            size=3,
            color=plot_df['Log Overlap Count'],
            colorscale='Viridis',
            colorbar=dict(title='Overlap Count'),
            showscale=True,
            symbol='square'
        ),
        name='Gene',
        text=[f"Clonal: {c}, Sub-Clonal: {s}, Overlap: {o}" for c, s, o in zip(
            plot_df['Clonal Count'], plot_df['Sub-Clonal Count'], plot_df['Overlap Count']
        )],
        hoverinfo='text'
    )
    fig.add_trace(scatter)

    # Add OLS Trendline
    x_ols = sm.add_constant(summary_df['Clonal Count'])
    y_ols = summary_df['Sub-Clonal Count']
    ols_model = sm.OLS(y_ols, x_ols).fit()
    ols_pred = ols_model.predict(x_ols)
    ols_line_df = pd.DataFrame({
        'Clonal Count': summary_df['Clonal Count'],
        'OLS Fit': ols_pred
    }).sort_values('Clonal Count')
    fig.add_trace(go.Scatter(
        x=ols_line_df['Clonal Count'],
        y=ols_line_df['OLS Fit'],
        mode='lines',
        name='OLS Trendline',
        line=dict(color='rgba(255,0,0,0.5)', width=2)
    ))

    # Add prediction line if available
    if pred_df is not None:
        fig.add_trace(go.Scatter(
            x=pred_df['Clonal Count'],
            y=pred_df['Predicted Sub-Clonal Count'],
            mode='lines',
            name='GLM Prediction',
            line=dict(color='green', dash='dash', width=3)
        ))

    if scale:
        lower, upper = negbin_quantiles(pred_df['Predicted Sub-Clonal Count'], scale)
        fig.add_trace(go.Scatter(
            x=list(pred_df['Clonal Count']) + list(pred_df['Clonal Count'])[::-1],
            y=list(lower) + list(upper)[::-1],
            fill='toself',
            fillcolor='rgba(0,100,80,0.2)',
            line=dict(color='rgba(255,255,255,0)'),
            hoverinfo='skip',
            showlegend=True,
            name='Quantile Range'
        ))
    
    # Update layout
    fig.update_layout(
        xaxis=dict(range=[0, 1200]),
        yaxis=dict(range=[0, 400]),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        title=plot_title
    )

    return fig

def plot_clonality_counts(data, title):
    """Creates a bar plot of clonality counts."""
    counts = data['Clonality'].value_counts().reset_index()
    counts.columns = ['Clonality', 'count']
    fig = px.bar(counts, x='Clonality', y='count', color='Clonality', title=title)
    return fig 