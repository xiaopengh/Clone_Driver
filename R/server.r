# server.R
library(shiny)
library(dplyr)
library(readr) # For readRDS, although base R has it
library(ggplot2)
library(MASS) # For glm.nb
library(tidyr) # For unnesting/pivoting if needed, and unstack equivalent
library(broom) # For tidying model summaries
library(DT)
library(hexbin) # For geom_hex
library(viridis) # Explicitly load viridis for the color scale

# === Server Logic ===
server <- function(input, output, session) {
  
  # Reactive value to trigger recalculations on refresh
  data_trigger <- reactiveVal(0)
  
  observeEvent(input$refresh, {
    data_trigger(runif(1)) # Update trigger with a new random value to force reactivity
  })
  
  # --- Load and Process Data ---
  # Reactive expression to load data, triggered by refresh button
  loaded_data <- reactive({
    data_trigger() # Depend on the trigger
    
    # Load the RDS files (Adjust path if necessary)
    # Assuming ui.R/server.R are in one folder and data is in ../data/
    data_binom <- readRDS("../data/mc3_binom_left.rds")
    data_Ztest <- readRDS("../data/mc3_Ztest_left.rds")
    
    # Ensure Clonality is a factor or character
    data_binom$Clonality <- as.character(data_binom$Clonality)
    data_Ztest$Clonality <- as.character(data_Ztest$Clonality)
    
    return(list(binom = data_binom, ztest = data_Ztest))
  })
  
  # --- Process Data for Summaries and Models ---
  processed_data <- reactive({
    data_list <- loaded_data()
    data_binom <- data_list$binom
    data_Ztest <- data_list$ztest
    
    # Equivalent to Python's groupby and value_counts().unstack()
    process_summary <- function(df) {
      df %>%
        group_by(Hugo_Symbol, Clonality) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        pivot_wider(names_from = Clonality, values_from = Count, values_fill = 0) %>%
        # Ensure columns exist even if no counts are present
        mutate(Clonal = if("Clonal" %in% names(.)) Clonal else 0,
               `Sub-Clonal` = if("Sub-Clonal" %in% names(.)) `Sub-Clonal` else 0) %>%
        rename(`Clonal Count` = Clonal, `Sub-Clonal Count` = `Sub-Clonal`)
    }
    
    clonal_summary_left <- process_summary(data_binom)
    clonal_summary_right <- process_summary(data_Ztest)
    
    # Create summary statistics
    create_stats <- function(summary_df, ndecimal = 2) {
      # Ensure the input columns exist before calculating stats
      req_cols <- c("Clonal Count", "Sub-Clonal Count")
      if (!all(req_cols %in% names(summary_df))) {
        stop("Missing required columns in summary dataframe for stats calculation.")
      }
      
      clonal_mean <- round(mean(summary_df$`Clonal Count`, na.rm = TRUE), ndecimal)
      clonal_median <- round(median(summary_df$`Clonal Count`, na.rm = TRUE), ndecimal)
      clonal_var <- round(var(summary_df$`Clonal Count`, na.rm = TRUE), ndecimal)
      
      subclonal_mean <- round(mean(summary_df$`Sub-Clonal Count`, na.rm = TRUE), ndecimal)
      subclonal_median <- round(median(summary_df$`Sub-Clonal Count`, na.rm = TRUE), ndecimal)
      subclonal_var <- round(var(summary_df$`Sub-Clonal Count`, na.rm = TRUE), ndecimal)
      
      data.frame(
        Statistic = c("Mean", "Median", "Variance"),
        Clonal = c(clonal_mean, clonal_median, clonal_var),
        `Sub-Clonal` = c(subclonal_mean, subclonal_median, subclonal_var),
        check.names = FALSE # Prevent R from changing 'Sub-Clonal' to 'Sub.Clonal'
      )
    }
    
    summary_stats_left <- create_stats(clonal_summary_left)
    summary_stats_right <- create_stats(clonal_summary_right)
    
    # --- Fit Models ---
    # Poisson
    poisson_model_left <- glm(`Sub-Clonal Count` ~ `Clonal Count`, data = clonal_summary_left, family = poisson)
    poisson_model_right <- glm(`Sub-Clonal Count` ~ `Clonal Count`, data = clonal_summary_right, family = poisson)
    
    # Negative Binomial
    negbin_model_left <- tryCatch({
      glm.nb(`Sub-Clonal Count` ~ `Clonal Count`, data = clonal_summary_left)
    }, error = function(e) {
      message("Error fitting NegBin Left: ", e$message); return(NULL)
    })
    
    negbin_model_right <- tryCatch({
      glm.nb(`Sub-Clonal Count` ~ `Clonal Count`, data = clonal_summary_right)
    }, error = function(e) {
      message("Error fitting NegBin Right: ", e$message); return(NULL)
    })
    
    
    return(list(
      clonal_summary_left = clonal_summary_left,
      clonal_summary_right = clonal_summary_right,
      summary_stats_left = summary_stats_left,
      summary_stats_right = summary_stats_right,
      poisson_model_left = poisson_model_left,
      poisson_model_right = poisson_model_right,
      negbin_model_left = negbin_model_left,
      negbin_model_right = negbin_model_right,
      data_binom = data_binom, # Pass original data too if needed
      data_Ztest = data_Ztest
    ))
  })
  
  # --- Render Outputs ---
  
  # Data Tables (showing a sample)
  output$data_table_left <- renderDT({
    datatable(head(loaded_data()$binom, 100), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  output$data_table_right <- renderDT({
    datatable(head(loaded_data()$ztest, 100), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  # Stats Tables
  output$stats_table_left <- renderDT({
    # Use check.names=FALSE in datatable as well if needed, although DT might handle it
    datatable(processed_data()$summary_stats_left, options = list(dom = 't', paging = FALSE), rownames = FALSE)
  })
  output$stats_table_right <- renderDT({
    datatable(processed_data()$summary_stats_right, options = list(dom = 't', paging = FALSE), rownames = FALSE)
  })
  
  # Clonality Count Bar Plots
  output$bar_plot_left <- renderPlot({
    data <- loaded_data()$binom
    ggplot(data, aes(x = factor(Clonality, levels = c('Clonal', 'Sub-Clonal')))) +
      geom_bar(aes(fill = Clonality), show.legend = FALSE) +
      scale_fill_manual(values = c("Clonal" = "#6FCF97", "Sub-Clonal" = "#EB5757")) +
      labs(title = "Clonality Counts - Binomial Test Dataset", x = "Clonality", y = "Count") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size=10),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9))
  })
  output$bar_plot_right <- renderPlot({
    data <- loaded_data()$ztest
    ggplot(data, aes(x = factor(Clonality, levels = c('Clonal', 'Sub-Clonal')))) +
      geom_bar(aes(fill = Clonality), show.legend = FALSE) +
      scale_fill_manual(values = c("Clonal" = "#6FCF97", "Sub-Clonal" = "#EB5757")) +
      labs(title = "Clonality Counts - Z Test Dataset", x = "Clonality", y = "Count") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size=10),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9))
  })
  
  # Model Summaries
  output$summary_poisson_left <- renderPrint({ summary(processed_data()$poisson_model_left) })
  output$summary_poisson_right <- renderPrint({ summary(processed_data()$poisson_model_right) })
  output$summary_negbin_left <- renderPrint({
    model <- processed_data()$negbin_model_left
    if (is.null(model)) "Negative Binomial model fitting failed." else summary(model)
  })
  output$summary_negbin_right <- renderPrint({
    model <- processed_data()$negbin_model_right
    if (is.null(model)) "Negative Binomial model fitting failed." else summary(model)
  })
  
  
  # --- Hex Plots ---
  # CORRECTED Helper function for hex plots
  create_hex_plot <- function(summary_data, model, title, xlim_val, ylim_val) {
    # Ensure input data has the correct columns
    req_cols <- c("Clonal Count", "Sub-Clonal Count")
    if (!all(req_cols %in% names(summary_data))) {
      print(head(summary_data)) # Print head for debugging
      stop("Missing required columns (`Clonal Count`, `Sub-Clonal Count`) in summary_data for hex plot.")
    }
    
    # Generate predictions
    pred_data <- NULL # Initialize pred_data
    if (!is.null(model)) {
      # **FIX 1**: Use backticks for column name in newdata for predict
      x_pred_values <- seq(0, xlim_val, length.out = 100)
      newdata_for_pred <- data.frame(x_pred_values)
      # Name the column EXACTLY as in the model formula
      names(newdata_for_pred) <- "`Clonal Count`"
      
      # Use tryCatch for prediction as well, just in case
      y_pred <- tryCatch({
        predict(model, newdata = newdata_for_pred, type = "response")
      }, error = function(e) {
        message("Error during prediction: ", e$message)
        return(NULL)
      })
      
      if (!is.null(y_pred)) {
        # Use original x values, not the named ones from the data frame
        pred_data <- data.frame(x = x_pred_values, y = y_pred)
      }
    }
    
    # **FIX 2**: Use backticks for column names with spaces in aes()
    p <- ggplot(summary_data, aes(x = `Clonal Count`, y = `Sub-Clonal Count`)) +
      geom_hex(bins = 50) + # Adjust bins as needed
      scale_fill_viridis_c(trans = "log10", name = "Log Count", option = "C") + # Explicitly load viridis and choose an option
      coord_cartesian(xlim = c(0, xlim_val), ylim = c(0, ylim_val)) + # Set limits AFTER stats/geoms
      # **FIX 3**: Use correct spacing in labs()
      labs(title = title, x = "Clonal Count", y = "Sub-Clonal Count") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size=10),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9),
            legend.title = element_text(size=9),
            legend.text = element_text(size=9))
    
    # Add prediction line if model and predictions exist
    if (!is.null(pred_data)) {
      # Ensure pred_data has finite values before plotting
      pred_data_finite <- filter(pred_data, is.finite(x) & is.finite(y))
      if(nrow(pred_data_finite) > 0) {
        p <- p + geom_line(data = pred_data_finite, aes(x = x, y = y), color = "#FF6F61", linewidth = 1, inherit.aes = FALSE)
      } else {
        message("Prediction data contained non-finite values, line not plotted.")
      }
    }
    
    return(p)
  }
  
  # Poisson Hex Plots
  output$hex_plot_poisson_left <- renderPlot({
    req(processed_data()) # Ensure processed_data is available
    p_data <- processed_data()
    create_hex_plot(p_data$clonal_summary_left, p_data$poisson_model_left,
                    "Clonal vs Subclonal Mutations (Poisson) - Binomial Test Dataset",
                    xlim_val = 300, ylim_val = 100) # Match Python limits
  })
  output$hex_plot_poisson_right <- renderPlot({
    req(processed_data())
    p_data <- processed_data()
    create_hex_plot(p_data$clonal_summary_right, p_data$poisson_model_right,
                    "Clonal vs Subclonal Mutations (Poisson) - Z Test Dataset",
                    xlim_val = 300, ylim_val = 100)
  })
  
  # Negative Binomial Hex Plots
  output$hex_plot_negbin_left <- renderPlot({
    req(processed_data())
    p_data <- processed_data()
    create_hex_plot(p_data$clonal_summary_left, p_data$negbin_model_left,
                    "Clonal vs Subclonal Mutations (NegBin) - Binomial Test Dataset",
                    xlim_val = 300, ylim_val = 100)
  })
  output$hex_plot_negbin_right <- renderPlot({
    req(processed_data())
    p_data <- processed_data()
    create_hex_plot(p_data$clonal_summary_right, p_data$negbin_model_right,
                    "Clonal vs Subclonal Mutations (NegBin) - Z Test Dataset",
                    xlim_val = 300, ylim_val = 100)
  })
  
} # End server