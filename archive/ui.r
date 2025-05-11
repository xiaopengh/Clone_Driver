# ui.R
library(shiny)
library(shinythemes)
library(DT)
# library(markdown) # No longer strictly needed for ui if using shiny::markdown

# Helper function to read markdown files safely
read_markdown_file <- function(path) {
  if (file.exists(path)) {
    # Read lines and collapse them into a single string
    paste(readLines(path, warn = FALSE), collapse = "\n")
  } else {
    paste("Error: Markdown file not found at", path) # Provide fallback text
  }
}

# Define UI
fluidPage(
  theme = shinytheme("flatly"), # Apply a theme for better aesthetics
  titlePanel("🧬 Tumor Clonality Dashboard"),
  
  # Add Refresh and Clear Cache buttons (Note: Cache clearing logic is simplified in Shiny)
  fluidRow(
    column(2, actionButton("refresh", "Refresh Data", class = "btn-primary")),
    # Note: Shiny's reactivity model often reduces the need for explicit caching and clearing.
    # A placeholder or simplified clear button can be added if specific non-reactive caching is implemented.
    # column(2, actionButton("clear_cache", "Clear Cache (Placeholder)", class = "btn-warning"))
  ),
  
  hr(), # Add a horizontal rule for separation
  
  # Section 1: Data Overview Tables
  h3("Dataset Overview"),
  fluidRow(
    column(6,
           h4("Binomial Test Generated Dataset Overview", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           DTOutput("data_table_left")
    ),
    column(6,
           h4("Z Test Generated Dataset Overview", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           DTOutput("data_table_right")
    )
  ),
  
  hr(),
  
  # Section 2: Summary Statistics Tables
  h3("Summary Statistics"),
  fluidRow(
    column(4,
           h4("Binomial Test Data Summary Statistics", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           DTOutput("stats_table_left")
    ),
    column(4,
           h4("Z Test Data Summary Statistics", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           DTOutput("stats_table_right")
    ),
    column(4,
           div(style = "background-color: #e6ffe6; padding: 10px; border-radius: 5px; font-size: 14px;", # Adjusted font size
               # Use shiny::markdown with readLines
               shiny::markdown(read_markdown_file("comments/stats_comment.md"))
           )
    )
  ),
  
  hr(),
  
  # Section 3: Clonality Counts Plots
  h2("Clonality Counts"),
  fluidRow(
    column(4, plotOutput("bar_plot_left")),
    column(4, plotOutput("bar_plot_right")),
    column(4,
           div(style = "background-color: #e6ffe6; padding: 10px; border-radius: 5px; font-size: 14px;",
               shiny::markdown(read_markdown_file("comments/clonality_comment.md"))
           )
    )
  ),
  
  hr(),
  
  # Section 4: Poisson Regression
  h2("Poisson Regression Analysis"),
  fluidRow(
    column(8,
           h3("Poisson Regression Model Summaries", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           fluidRow(
             column(6, verbatimTextOutput("summary_poisson_left")),
             column(6, verbatimTextOutput("summary_poisson_right"))
           )
    ),
    column(4,
           div(style = "background-color: #e6ffe6; padding: 10px; border-radius: 5px; font-size: 14px;",
               shiny::markdown(read_markdown_file("comments/poisson_comment.md"))
           )
    )
  ),
  fluidRow(
    column(12, h3("Poisson Regression Hexplots Comparison", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;")),
    column(6, plotOutput("hex_plot_poisson_left")),
    column(6, plotOutput("hex_plot_poisson_right"))
  ),
  
  hr(),
  
  # Section 5: Negative Binomial Regression
  h2("Negative Binomial Regression Analysis"),
  fluidRow(
    column(8,
           h3("Negative Binomial Regression Model Summaries", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;"),
           fluidRow(
             column(6, verbatimTextOutput("summary_negbin_left")),
             column(6, verbatimTextOutput("summary_negbin_right"))
           )
    ),
    column(4,
           div(style = "background-color: #e6ffe6; padding: 10px; border-radius: 5px; font-size: 14px;",
               shiny::markdown(read_markdown_file("comments/negbin_comment.md"))
           )
    )
  ),
  fluidRow(
    column(12, h3("Negative Binomial Regression Hexplots Comparison", style = "background-color: #F0F8FF; padding: 10px; border-radius: 5px;")),
    column(6, plotOutput("hex_plot_negbin_left")),
    column(6, plotOutput("hex_plot_negbin_right"))
  ),
  
  # Add padding at the bottom
  tags$div(style = "padding-bottom: 50px;")
  
) # End fluidPage
