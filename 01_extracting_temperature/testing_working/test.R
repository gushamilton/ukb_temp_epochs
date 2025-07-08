pacman::p_load(tidyverse, arrow)
d <- read_csv("/Users/fh6520/Downloads/temperature_epochs_first100.csv.gz")
d2 <- read_parquet("/Users/fh6520/Downloads/temperature_epochs_first100.parquet")
d2
count(d2, id)

d2 %>% 
  filter(id == "1000108_90001_1_0") %>% 
  ggplot(aes(x = time, y = temp)) +
  geom_line()




# Plot the first 10 unique IDs' temperature traces, faceted by ID
library(signal) # for spec.pgram, if not already loaded

first10_ids <- d %>% distinct(id) %>% slice_head(n = 10) %>% pull(id)

# Prepare data for time series and FFT
d_first10 <- d %>%
  filter(id %in% first10_ids) %>%
  select(id, time, temp) %>%
  mutate(time = as.POSIXct(time, origin = "1970-01-01"))

# Plot the temperature traces as before
p <- d_first10 %>%
  ggplot(aes(x = time, y = temp)) +
  geom_line() +
  facet_wrap(~ id, scales = "free_x") +
  ggtitle("Temperature traces for first 10 IDs") +
  theme(strip.text = element_blank())
ggsave("temperature_first10_faceted.png", p, width = 14, height = 8)

# Compute and plot simple Fourier transform (magnitude spectrum) for each ID
library(tibble)
library(patchwork) # for combining plots

fft_plots <- lapply(first10_ids, function(cur_id) {
  df_id <- d_first10 %>% filter(id == cur_id) %>% arrange(time)
  temp_vec <- df_id$temp
  n <- length(temp_vec)
  if (n < 8) return(NULL) # skip if too short

  # Remove NAs for FFT, pad with mean if needed
  temp_vec[is.na(temp_vec)] <- mean(temp_vec, na.rm = TRUE)
  temp_vec <- temp_vec - mean(temp_vec) # remove DC

  # Compute FFT
  fft_res <- fft(temp_vec)
  freq <- seq(0, 0.5, length.out = floor(n/2) + 1) # up to Nyquist
  mag <- Mod(fft_res)[1:length(freq)]
  # Normalize
  mag <- mag / max(mag, na.rm = TRUE)

  tibble(
    id = cur_id,
    freq = freq,
    mag = mag
  )
})

fft_df <- bind_rows(fft_plots)

# Plot the FFT magnitude spectra, faceted by ID
p_fft <- fft_df %>%
  ggplot(aes(x = freq, y = mag)) +
  geom_line(color = "purple") +
  facet_wrap(~ id, scales = "free_y") +
  labs(
    x = "Frequency (cycles per sample)",
    y = "Normalized magnitude",
    title = "FFT magnitude spectra for first 10 IDs"
  )
ggsave("temperature_first10_fft_faceted.png", p_fft, width = 14, height = 8)

# Optionally, show both plots side by side in RStudio
p + p_fft

# Helper: safely convert time to numeric, set non-finite to NA
safe_time_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- NA_real_
  out
}

# Group every 2 rows and average temp
d2 <- d %>%
  mutate(time_num = safe_time_numeric(time)) %>%
  group_by(id, group = (row_number() - 1) %/% 2) %>%
  summarise(
    time = mean(time_num, na.rm = TRUE),  # average the time as numeric, skip NA/Inf
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(time = as.POSIXct(time, origin = "1970-01-01"))

# Only plot if all time values are finite
if (all(is.finite(d2$time))) {
  p2 <- d2 %>%
    ggplot(aes(x = time, y = temp)) +
    geom_line(color = "blue") +
    ggtitle("Averaged every 2 rows (1-min epochs)")
  ggsave("temperature_2row.png", p2, width = 10, height = 5)
} else {
  warning("Non-finite values in d2$time; skipping temperature_2row.png plot.")
}

# Group every 5 rows and average temp
d5 <- d %>%
  mutate(time_num = safe_time_numeric(time)) %>%
  group_by(id, group = (row_number() - 1) %/% 5) %>%
  summarise(
    time = mean(time_num, na.rm = TRUE),
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(time = as.POSIXct(time, origin = "1970-01-01"))

if (all(is.finite(d5$time))) {
  p5 <- d5 %>%
    ggplot(aes(x = time, y = temp)) +
    geom_line(color = "red") +
    ggtitle("Averaged every 5 rows (2.5-min epochs)")
  ggsave("temperature_5row.png", p5, width = 10, height = 5)
} else {
  warning("Non-finite values in d5$time; skipping temperature_5row.png plot.")
}

# Combine all three for comparison, only if all plots are available
library(patchwork)
plots <- list()
plots[[1]] <- p + ggtitle("Original (30s epochs)")
if (exists("p2")) plots[[2]] <- p2
if (exists("p5")) plots[[3]] <- p5

if (length(plots) == 3) {
  p_compare <- plots[[1]] + plots[[2]] + plots[[3]]
  ggsave("temperature_compare.png", p_compare, width = 15, height = 10)
} else {
  warning("Not all plots available for comparison; skipping temperature_compare.png.")
}




d %>% view()
