library(ggplot2)
library(dplyr)
library(data.table)


# mc3_raw[, .(count_of_rows = .N), by = Transcript_ID]

summary_dt <- rbindlist(summaries)

summary_dt %>% nrow()
summary_dt %>% na.omit() %>% nrow()

g <- ggplot(data = na.omit(summary_dt), aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal))

g + stat_bin_2d(mapping = aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal), bins = 200) +
  geom_point(x = 1, y = 1, color = "red", size = 3, shape = 17) + # Corrected line
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw()

g + stat_bin_2d(mapping = aes(fill = after_stat(count)), bins = 25) +
  geom_text(stat = "bin_2d", mapping = aes(label = after_stat(count)), size = 3) +
  theme_bw()