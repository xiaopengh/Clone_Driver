library(ggplot2)

# mc3_raw[, .(count_of_rows = .N), by = Transcript_ID]

g <- ggplot(data = mc3_raw)

g + geom_histogram(data = mc3_raw[, .(count_of_rows = .N), by = Transcript_ID]) + 
  aes(x = count_of_rows) +
  labs(title = "Histogram of Count of Rows by Transcript ID",
       x = "Count of Rows",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

