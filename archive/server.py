import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import os
from datetime import datetime, timedelta

# === Centralized Plot Titles ===
# All plot titles are defined here for easy management
PLOT_TITLES = {
    # Bar chart titles
    'bar_left': 'Clonality Counts - Binomial Test Dataset',
    'bar_right': 'Clonality Counts - Z Test Dataset',
    
    # Hexplot titles
    'hex_left': 'Clonal vs Subclonal Mutations (Poisson) - Binomial Test Dataset',
    'hex_right': 'Clonal vs Subclonal Mutations (Poisson) - Z Test Dataset',
    'hex_left_negbin': 'Clonal vs Subclonal Mutations (NegBin) - Binomial Test Dataset',
    'hex_right_negbin': 'Clonal vs Subclonal Mutations (NegBin) - Z Test Dataset',
    
    # Model summary titles
    'summary_poisson_left': 'Binomial Test Dataset Poisson Regression',
    'summary_poisson_right': 'Z Test Dataset Poisson Regression',
    'summary_negbin_left': 'Binomial Test Dataset Negative Binomial Regression',
    'summary_negbin_right': 'Z Test Dataset Negative Binomial Regression',
}

# Set default font sizes globally for all plots
plt.rcParams.update({
    'axes.titlesize': 10,    # Titles
    'axes.labelsize': 10,    # X and Y labels
    'xtick.labelsize': 9,    # X-axis tick labels
    'ytick.labelsize': 9,    # Y-axis tick labels
    'legend.fontsize': 9,    # Legend
    'font.size': 10          # General font size (fallback/default)
})

# === Load and process data ===
def load_data():
    data_binom = pkl.load(open("data/mc3_binom_left.pkl", "rb")).reset_index(drop=True)
    data_binom["Clonality"] = data_binom["Clonality"].astype(str)

    data_Ztest = pkl.load(open("data/mc3_Ztest_left.pkl", "rb")).reset_index(drop=True)
    data_Ztest["Clonality"] = data_Ztest["Clonality"].astype(str)
    
    return data_binom, data_Ztest

# Function to process each dataset
def process_data(data):
    summary = data.groupby("Hugo_Symbol")["Clonality"].value_counts().unstack(fill_value=0)
    summary.columns.name = None
    summary = summary.rename(columns={"Clonal": "Clonal Count", "Sub-Clonal": "Sub-Clonal Count"})
    return summary

# Function to create summary statistics
def create_summary_stats(clonal_summary, ndecimal=2):
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

# Function to create clonality counts plot
def plot_clonality_counts(data, title, figsize=(5, 4)):
    orders = ['Clonal', 'Sub-Clonal']
    counts = data["Clonality"].value_counts().reindex(orders, fill_value=0)
    fig, ax = plt.subplots(figsize=figsize)
    counts.plot(kind='bar', color=['#6FCF97', '#EB5757'], ax=ax)
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel("Clonality")
    ax.set_ylabel("Count")
    plt.xticks(rotation=0)
    return fig

# Function for hexplot visualization
def plot_hexbin_poisson(summary, model, title, lims=(400, 400), figsize=(5, 4)):
    x = summary["Clonal Count"]
    y = summary["Sub-Clonal Count"]

    xlim, ylim = lims

    x_plot = np.linspace(0, xlim, 400)
    y_pred = model.predict(sm.add_constant(x_plot))

    fig, ax = plt.subplots(figsize=figsize)
    hb = ax.hexbin(x, y, gridsize=2000, bins='log', cmap="Blues", mincnt=1)
    ax.plot(x_plot, y_pred, color='#FF6F61', linewidth=2, label='Model Prediction')
    fig.colorbar(hb, ax=ax, label='Log Count in Bin')
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.set_xlabel("Clonal Count")
    ax.set_ylabel("Sub-Clonal Count")
    ax.set_title(title, fontweight='bold')
    ax.legend()
    return fig

# Function to prepare table data with specific limits
def prepare_table_data(data, percentage=0.3):
    """Prepare data for display in tables, taking only a percentage of the data"""
    return data.head(int(len(data) * percentage))

# Function to generate model summary text
def get_model_summary_text(model, title):
    """Generate formatted model summary text for display"""
    return f"### {title}\n```\n{model.summary()}\n```"

# Function to create all visualizations
def create_all_visualizations(data_dict, figsize=(5, 4)):
    """Create all visualizations needed for the dashboard"""
    visualizations = {}
    
    # Bar charts
    visualizations['bar_left'] = plot_clonality_counts(
        data_dict['data_binom'], PLOT_TITLES['bar_left'], figsize
    )
    
    visualizations['bar_right'] = plot_clonality_counts(
        data_dict['data_Ztest'], PLOT_TITLES['bar_right'], figsize
    )
    
    # Poisson hexplots
    visualizations['hex_left'] = plot_hexbin_poisson(
        data_dict['clonal_summary_left'],
        data_dict['poisson_model_left'],
        PLOT_TITLES['hex_left'],
        lims=(300, 100),
        figsize=figsize
    )
    
    visualizations['hex_right'] = plot_hexbin_poisson(
        data_dict['clonal_summary_right'],
        data_dict['poisson_model_right'],
        PLOT_TITLES['hex_right'],
        lims=(300, 100),
        figsize=figsize
    )
    
    # Negative binomial hexplots
    visualizations['hex_left_negbin'] = plot_hexbin_poisson(
        data_dict['clonal_summary_left'],
        data_dict['negbin_model_left'],
        PLOT_TITLES['hex_left_negbin'],
        lims=(300, 100),
        figsize=figsize
    )
    
    visualizations['hex_right_negbin'] = plot_hexbin_poisson(
        data_dict['clonal_summary_right'],
        data_dict['negbin_model_right'],
        PLOT_TITLES['hex_right_negbin'],
        lims=(300, 100),
        figsize=figsize
    )
    
    return visualizations

# Function to prepare model summaries
def prepare_model_summaries(data_dict):
    """Prepare model summary texts for display"""
    summaries = {
        'poisson_left': get_model_summary_text(
            data_dict['poisson_model_left'], 
            PLOT_TITLES['summary_poisson_left']
        ),
        'poisson_right': get_model_summary_text(
            data_dict['poisson_model_right'], 
            PLOT_TITLES['summary_poisson_right']
        ),
        'negbin_left': get_model_summary_text(
            data_dict['negbin_model_left'], 
            PLOT_TITLES['summary_negbin_left']
        ),
        'negbin_right': get_model_summary_text(
            data_dict['negbin_model_right'], 
            PLOT_TITLES['summary_negbin_right']
        )
    }
    return summaries

# === Caching functions ===
CACHE_FILE = "data/dashboard_cache.pkl"
CACHE_EXPIRY_HOURS = 0.25  # Cache expires after 24 hours

def save_to_cache(data_dict):
    """Save processed data to cache file"""
    try:
        # Add timestamp to cached data
        cache_data = {
            'timestamp': datetime.now(),
            'data': data_dict
        }
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
        
        # Save to cache file
        with open(CACHE_FILE, 'wb') as f:
            pkl.dump(cache_data, f)
        
        print(f"Data cached successfully at {datetime.now()}")
        return True
    except Exception as e:
        print(f"Error saving to cache: {e}")
        return False

def is_cache_valid():
    """Check if cache exists and is still valid"""
    if not os.path.exists(CACHE_FILE):
        return False
    
    try:
        with open(CACHE_FILE, 'rb') as f:
            cache_data = pkl.load(f)
        
        # Check if cache has timestamp
        if 'timestamp' not in cache_data:
            return False
        
        # Check if cache is expired
        cache_time = cache_data['timestamp']
        expiry_time = datetime.now() - timedelta(hours=CACHE_EXPIRY_HOURS)
        
        return cache_time > expiry_time
    except Exception as e:
        print(f"Error checking cache: {e}")
        return False

def load_from_cache():
    """Load data from cache if it exists and is valid"""
    if not is_cache_valid():
        return None
    
    try:
        with open(CACHE_FILE, 'rb') as f:
            cache_data = pkl.load(f)
        
        print(f"Loaded data from cache created at {cache_data['timestamp']}")
        return cache_data['data']
    except Exception as e:
        print(f"Error loading from cache: {e}")
        return None

def clear_cache():
    """Remove the cache file if it exists"""
    if os.path.exists(CACHE_FILE):
        os.remove(CACHE_FILE)
        print(f"Cache cleared at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    else:
        print("No cache file found to clear")

def get_dashboard_data(force_recalculate=False):
    """Get dashboard data, either from cache or by calculating it"""
    # Try to load from cache first (unless force_recalculate is True)
    if not force_recalculate:
        cached_data = load_from_cache()
        if cached_data is not None:
            return cached_data
    
    # If no valid cache exists or force_recalculate is True, calculate the data
    print("Calculating dashboard data...")
    data_dict = initialize_data_and_models()
    
    # Save the calculated data to cache
    save_to_cache(data_dict)
    
    return data_dict

# Initialize data and models
def initialize_data_and_models():
    # Load data
    data_binom, data_Ztest = load_data()
    
    # Process data
    clonal_summary_left = process_data(data_binom)
    clonal_summary_right = process_data(data_Ztest)
    
    # Create summary statistics
    summary_stats_left = create_summary_stats(clonal_summary_left)
    summary_stats_right = create_summary_stats(clonal_summary_right)
    
    # Fit models
    poisson_model_left = fit_poisson_model(clonal_summary_left)
    poisson_model_right = fit_poisson_model(clonal_summary_right)
    
    negbin_model_left = fit_negbin_model(clonal_summary_left)
    negbin_model_right = fit_negbin_model(clonal_summary_right)
    
    # Create data dictionary
    data_dict = {
        'data_binom': data_binom,
        'data_Ztest': data_Ztest,
        'clonal_summary_left': clonal_summary_left,
        'clonal_summary_right': clonal_summary_right,
        'summary_stats_left': summary_stats_left,
        'summary_stats_right': summary_stats_right,
        'poisson_model_left': poisson_model_left,
        'poisson_model_right': poisson_model_right,
        'negbin_model_left': negbin_model_left,
        'negbin_model_right': negbin_model_right
    }
    
    # Prepare table data
    data_dict['table_data_left'] = prepare_table_data(data_binom)
    data_dict['table_data_right'] = prepare_table_data(data_Ztest)
    
    # Create visualizations
    data_dict['visualizations'] = create_all_visualizations(data_dict)
    
    # Prepare model summaries
    data_dict['model_summaries'] = prepare_model_summaries(data_dict)
    
    return data_dict