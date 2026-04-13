# ============================================================
# build_banking_indicator_full_seasonal.R
# Full banking sector indicator model with seasonality
# Colombia
# ------------------------------------------------------------
# Publication-ready baseline model:
# - 6 observed banking indicators
# - common local linear trend (level + slope)
# - monthly seasonal component
# - trend level is interpreted as the synthetic banking indicator
# ============================================================

rm(list = ls())

# -----------------------------
# 1) Packages
# -----------------------------
required_packages <- c("KFAS", "readr", "dplyr", "ggplot2")

to_install <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(required_packages, library, character.only = TRUE))

# -----------------------------
# 2) Project paths
# -----------------------------
project_root <- getwd()

raw_data_path      <- file.path(project_root, "data", "raw", "IndicadoresActualizados.txt")
processed_data_dir <- file.path(project_root, "data", "processed")
figures_dir        <- file.path(project_root, "output", "figures")

dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(raw_data_path))

# -----------------------------
# 3) Read source data
# -----------------------------
bank_data <- read.csv(
  file   = raw_data_path,
  header = TRUE,
  dec    = ",",
  sep    = "\t"
)

# -----------------------------
# 4) Select the six indicators
# -----------------------------
y_raw <- cbind(
  CodCxAdM     = bank_data$CodCxAdM,
  Inv_CadCxAdM = 1 / bank_data$CadCxAdM,
  RdS          = bank_data$RdS,
  ROE          = bank_data$ROE,
  Inv_AdMF     = 1 / bank_data$AdMF,
  Inv_EPxAdM   = 1 / bank_data$EPxAdM
)

valid_rows <- complete.cases(y_raw) & apply(y_raw, 1, function(x) all(is.finite(x)))
y_raw <- y_raw[valid_rows, , drop = FALSE]

date_vector <- NULL
if ("Fecha" %in% names(bank_data)) {
  date_vector <- bank_data$Fecha[valid_rows]
}

# Standardize indicators
Y <- scale(y_raw)

# Convert to time-series matrix
Y_ts <- ts(Y, start = c(1996, 12), frequency = 12)

# -----------------------------
# 5) Build a correct univariate structural template
# -----------------------------
# IMPORTANT:
# - H here must be 1x1 because this is a univariate template
# - degree = 2 means local linear trend (level + slope)
# - seasonal = dummy monthly seasonality

base_model <- SSModel(
  Y_ts[, 1] ~
    SSMtrend(degree = 2, Q = list(matrix(NA), matrix(NA))) +
    SSMseasonal(period = 12, sea.type = "dummy", Q = matrix(NA)),
  H = matrix(NA)
)

# Number of observed series
n_series <- ncol(Y_ts)

# Number of latent states taken directly from the template
state_dim <- dim(base_model$T)[1]

# Extract template observation row
Z_template <- base_model$Z[1, , 1]

# Expand observation equation to all 6 observed series
Z_full <- array(0, dim = c(n_series, state_dim, nrow(Y_ts)))
for (t in 1:nrow(Y_ts)) {
  for (i in 1:n_series) {
    Z_full[i, , t] <- Z_template
  }
}

# -----------------------------
# 6) Build multivariate common-factor model
# -----------------------------
model <- SSModel(
  Y_ts ~ -1 +
    SSMcustom(
      Z     = Z_full,
      T     = base_model$T,
      R     = base_model$R,
      Q     = base_model$Q,
      a1    = base_model$a1,
      P1    = base_model$P1,
      P1inf = base_model$P1inf
    ),
  H = diag(NA, n_series)
)

# -----------------------------
# 7) Update function for MLE
# -----------------------------
# Disturbance variances in Q:
# Q[1,1] = level disturbance variance
# Q[2,2] = slope disturbance variance
# Q[3,3] = seasonal disturbance variance
#
# Measurement error variances:
# diagonal elements of H for the 6 observed indicators

update_model <- function(pars, model) {
  
  # State disturbance variances
  model$Q[1, 1, 1] <- exp(pars[1])  # level variance
  model$Q[2, 2, 1] <- exp(pars[2])  # slope variance
  model$Q[3, 3, 1] <- exp(pars[3])  # seasonal variance
  
  # Measurement error variances
  diag(model$H[, , 1]) <- exp(pars[4:(3 + n_series)])
  
  model
}

# -----------------------------
# 8) Estimate model
# -----------------------------
init_pars <- log(c(
  0.01,  # level variance
  0.01,  # slope variance
  0.01,  # seasonal variance
  rep(0.05, n_series)  # measurement variances
))

fit <- fitSSM(
  inits    = init_pars,
  model    = model,
  updatefn = update_model,
  method   = "BFGS"
)

# Smoothed states
kfs <- KFS(fit$model, smoothing = c("state", "signal"))

# -----------------------------
# 9) Extract components
# -----------------------------
# In the local linear trend model:
# state 1 = level
# state 2 = slope
# remaining states include seasonality

trend_level <- kfs$alphahat[, 1]
trend_slope <- kfs$alphahat[, 2]

# Common signal from the full observation equation
common_signal <- numeric(nrow(Y_ts))
for (t in 1:nrow(Y_ts)) {
  common_signal[t] <- as.numeric(Z_full[1, , t] %*% kfs$alphahat[t, ])
}

# Seasonal component = common signal - level
seasonal_component <- common_signal - trend_level

# Synthetic banking indicator = trend level
banking_indicator <- trend_level

# -----------------------------
# 10) Save outputs
# -----------------------------
results_list <- list(
  dates              = date_vector,
  standardized_data  = Y,
  trend_level        = trend_level,
  trend_slope        = trend_slope,
  seasonal_component = seasonal_component,
  common_signal      = common_signal,
  banking_indicator  = banking_indicator,
  fitted_model       = fit,
  smoothed_results   = kfs
)

saveRDS(
  results_list,
  file = file.path(processed_data_dir, "banking_indicator_full_seasonal.rds")
)

write.csv(
  data.frame(
    date               = if (!is.null(date_vector)) date_vector else seq_along(banking_indicator),
    banking_indicator  = banking_indicator,
    common_signal      = common_signal,
    seasonal_component = seasonal_component,
    trend_slope        = trend_slope
  ),
  file = file.path(processed_data_dir, "banking_indicator_full_seasonal.csv"),
  row.names = FALSE
)

# -----------------------------
# 11) Plot common signal, trend and seasonality
# -----------------------------
png(
  filename = file.path(figures_dir, "banking_indicator_full_seasonal.png"),
  width = 1600,
  height = 1000,
  res = 150
)

ts.plot(
  ts(common_signal, start = c(1996, 12), frequency = 12),
  col = "gray40",
  lwd = 2,
  ylab = "Standardized scale",
  main = "Common Signal, Trend, and Seasonal Component"
)

lines(
  ts(banking_indicator, start = c(1996, 12), frequency = 12),
  col = "black",
  lwd = 3
)

lines(
  ts(seasonal_component, start = c(1996, 12), frequency = 12),
  col = "red",
  lwd = 2,
  lty = 2
)

legend(
  "bottomright",
  legend = c(
    "Common signal (trend + seasonality)",
    "Trend = banking indicator",
    "Seasonal component"
  ),
  col = c("gray40", "black", "red"),
  lty = c(1, 1, 2),
  lwd = c(2, 3, 2)
)

dev.off()

# -----------------------------
# 12) Plot indicator vs observed variables
# -----------------------------
png(
  filename = file.path(figures_dir, "banking_indicator_vs_observed_full_model.png"),
  width = 1600,
  height = 1000,
  res = 150
)

ts.plot(
  ts(banking_indicator, start = c(1996, 12), frequency = 12),
  col = "black",
  lwd = 3,
  ylab = "Standardized scale",
  main = "Synthetic Banking Indicator (Trend Only) vs Observed Indicators"
)

for (j in 1:ncol(Y_ts)) {
  lines(Y_ts[, j], lty = 2, lwd = 1)
}

legend(
  "bottomright",
  legend = c("Banking indicator", colnames(Y_ts)),
  col = c("black", rep("gray50", ncol(Y_ts))),
  lty = c(1, rep(2, ncol(Y_ts))),
  lwd = c(3, rep(1, ncol(Y_ts))),
  cex = 0.8
)

dev.off()

# -----------------------------
# 13) Residual diagnostics
# -----------------------------
residuals_model <- residuals(kfs, type = "recursive")

png(
  filename = file.path(figures_dir, "residual_diagnostics.png"),
  width = 1600,
  height = 1200,
  res = 150
)

par(mfrow = c(2, 2))
plot(residuals_model, main = "Recursive residuals")
acf(residuals_model, main = "ACF of residuals")
pacf(residuals_model, main = "PACF of residuals")
hist(residuals_model, main = "Histogram of residuals", xlab = "Residuals")

dev.off()
par(mfrow = c(1, 1))

# Ljung-Box per series
ljung_box_results <- sapply(1:ncol(residuals_model), function(j) {
  Box.test(na.omit(residuals_model[, j]), lag = 12, type = "Ljung-Box")$p.value
})

ljung_box_table <- data.frame(
  series = colnames(Y_ts),
  ljung_box_pvalue = ljung_box_results
)

write.csv(
  ljung_box_table,
  file = file.path(processed_data_dir, "ljung_box_results.csv"),
  row.names = FALSE
)

cat("Full structural model with seasonality estimated successfully.\n")
print(ljung_box_table)

################################################################################

# =========================================================
# MODEL 2: Comparable common-factor model
# - Common local linear trend
# - Trigonometric seasonality
# - Same multivariate common loading structure as Model 1
# =========================================================

# -------------------------------
# 1) Base univariate model
# -------------------------------
base_model2 <- SSModel(
  Y_ts[, 1] ~
    SSMtrend(degree = 2, Q = list(matrix(NA), matrix(NA))) +
    SSMseasonal(period = 12, sea.type = "trigonometric", Q = matrix(NA)),
  H = matrix(NA)
)

# Dimensions
state_dim2  <- dim(base_model2$T)[1]
Z_template2 <- base_model2$Z[1, , 1]

# Optional quick check
print(dim(base_model2$T))
print(dim(base_model2$R))
print(dim(base_model2$Q))

# -------------------------------
# 2) Expand to multivariate common-factor model
# -------------------------------
Z_full2 <- array(0, dim = c(n_series, state_dim2, nrow(Y_ts)))

for (t in 1:nrow(Y_ts)) {
  for (i in 1:n_series) {
    Z_full2[i, , t] <- Z_template2
  }
}

model2 <- SSModel(
  Y_ts ~ -1 +
    SSMcustom(
      Z     = Z_full2,
      T     = base_model2$T,
      R     = base_model2$R,
      Q     = base_model2$Q,
      a1    = base_model2$a1,
      P1    = base_model2$P1,
      P1inf = base_model2$P1inf
    ),
  H = diag(NA, n_series)
)

# -------------------------------
# 3) Safe update function
# -------------------------------
# pars[1] = level variance
# pars[2] = slope variance
# pars[3] = seasonal variance
# pars[4:(3+n_series)] = measurement variances
update_model2 <- function(pars, model) {
  
  eps <- 1e-6
  
  # Level variance
  model$Q[1, 1, 1] <- exp(pars[1]) + eps
  
  # Slope variance
  model$Q[2, 2, 1] <- exp(pars[2]) + eps
  
  # Seasonal variance:
  # apply the same variance to all seasonal states
  if (dim(model$Q)[1] > 2) {
    seasonal_indices <- 3:dim(model$Q)[1]
    for (i in seasonal_indices) {
      model$Q[i, i, 1] <- exp(pars[3]) + eps
    }
  }
  
  # Measurement noise variances
  diag(model$H[, , 1]) <- exp(pars[4:(3 + n_series)]) + eps
  
  model
}

# -------------------------------
# 4) Safer initialization
# -------------------------------
init_pars_model2 <- log(c(
  0.01,             # level variance
  0.01,             # slope variance
  0.01,             # seasonal variance
  rep(0.1, n_series) # measurement variances
))

# -------------------------------
# 5) Fit model 2
# -------------------------------
fit_model2 <- fitSSM(
  inits    = init_pars_model2,
  model    = model2,
  updatefn = update_model2,
  method   = "BFGS",
  control  = list(maxit = 500)
)

kfs_model2 <- KFS(
  fit_model2$model,
  smoothing = c("state", "signal")
)

# -------------------------------
# 6) Extract components
# -------------------------------
# First state = common trend level
trend_level_model2 <- kfs_model2$alphahat[, 1]

# Second state = common trend slope
trend_slope_model2 <- kfs_model2$alphahat[, 2]

# Common signal from observation equation
common_signal_model2 <- numeric(nrow(Y_ts))
for (t in 1:nrow(Y_ts)) {
  common_signal_model2[t] <- as.numeric(Z_full2[1, , t] %*% kfs_model2$alphahat[t, ])
}

# Seasonal component = common signal - trend level
seasonal_component_model2 <- common_signal_model2 - trend_level_model2

# Synthetic banking indicator = common trend level
banking_indicator_model2 <- trend_level_model2

# -------------------------------
# 7) Save outputs
# -------------------------------
results_list_model2 <- list(
  trend_level        = trend_level_model2,
  trend_slope        = trend_slope_model2,
  common_signal      = common_signal_model2,
  seasonal_component = seasonal_component_model2,
  banking_indicator  = banking_indicator_model2,
  fitted_model       = fit_model2,
  smoothed_results   = kfs_model2
)

saveRDS(
  results_list_model2,
  file = file.path(processed_data_dir, "banking_indicator_model2_trigonometric.rds")
)

write.csv(
  data.frame(
    date               = if (!is.null(date_vector)) date_vector else seq_along(banking_indicator_model2),
    banking_indicator  = banking_indicator_model2,
    common_signal      = common_signal_model2,
    seasonal_component = seasonal_component_model2,
    trend_slope        = trend_slope_model2
  ),
  file = file.path(processed_data_dir, "banking_indicator_model2_trigonometric.csv"),
  row.names = FALSE
)

# -------------------------------
# 8) Plot common signal, trend and seasonality
# -------------------------------
png(
  filename = file.path(figures_dir, "model2_common_factor.png"),
  width = 1600,
  height = 1000,
  res = 150
)

ts.plot(
  ts(common_signal_model2, start = c(1996, 12), frequency = 12),
  col = "gray40",
  lwd = 2,
  ylab = "Standardized scale",
  main = "Model 2: Common Signal, Trend, and Seasonal Component"
)

lines(
  ts(banking_indicator_model2, start = c(1996, 12), frequency = 12),
  col = "black",
  lwd = 3
)

lines(
  ts(seasonal_component_model2, start = c(1996, 12), frequency = 12),
  col = "blue",
  lwd = 2,
  lty = 2
)

legend(
  "bottomright",
  legend = c(
    "Common signal (trend + seasonality)",
    "Trend = banking indicator",
    "Seasonal component"
  ),
  col = c("gray40", "black", "blue"),
  lty = c(1, 1, 2),
  lwd = c(2, 3, 2)
)

dev.off()

# -------------------------------
# 9) Plot indicator vs observed variables
# -------------------------------
png(
  filename = file.path(figures_dir, "model2_vs_observed.png"),
  width = 1600,
  height = 1000,
  res = 150
)

ts.plot(
  ts(banking_indicator_model2, start = c(1996, 12), frequency = 12),
  col = "black",
  lwd = 3,
  ylab = "Standardized scale",
  main = "Model 2: Synthetic Banking Indicator vs Observed Series"
)

for (j in 1:ncol(Y_ts)) {
  lines(Y_ts[, j], lty = 2, lwd = 1)
}

legend(
  "bottomright",
  legend = c("Synthetic indicator", colnames(Y_ts)),
  col = c("black", rep("gray50", ncol(Y_ts))),
  lty = c(1, rep(2, ncol(Y_ts))),
  lwd = c(3, rep(1, ncol(Y_ts))),
  cex = 0.8
)

dev.off()

# -------------------------------
# 10) Residual diagnostics
# -------------------------------
residuals_model2 <- residuals(kfs_model2, type = "recursive")

png(
  filename = file.path(figures_dir, "model2_residual_diagnostics.png"),
  width = 1600,
  height = 1200,
  res = 150
)

par(mfrow = c(2, 2))
plot(residuals_model2, main = "Recursive residuals (Model 2)")
acf(residuals_model2, main = "ACF of residuals (Model 2)")
pacf(residuals_model2, main = "PACF of residuals (Model 2)")
hist(residuals_model2, main = "Histogram of residuals (Model 2)", xlab = "Residuals")
dev.off()

par(mfrow = c(1, 1))

# -------------------------------
# 11) Ljung-Box per series
# -------------------------------
ljung_box_model2_table <- data.frame(
  series = colnames(Y_ts),
  p_value = sapply(1:ncol(residuals_model2), function(j) {
    Box.test(na.omit(residuals_model2[, j]), lag = 12, type = "Ljung-Box")$p.value
  })
)

write.csv(
  ljung_box_model2_table,
  file = file.path(processed_data_dir, "ljung_box_model2.csv"),
  row.names = FALSE
)

print(ljung_box_model2_table)

# -------------------------------
# 12) Model comparison
# -------------------------------
# -------------------------------
# 12) Model comparison (manual AIC/BIC)
# -------------------------------

logLik_m1 <- as.numeric(logLik(fit$model))
logLik_m2 <- as.numeric(logLik(fit_model2$model))

# Number of estimated parameters:
# Model 1:
# 3 state variances + 6 measurement variances = 9
k_m1 <- 3 + n_series

# Model 2:
# 3 state variances + 6 measurement variances = 9
k_m2 <- 3 + n_series

# Sample size
n_obs <- nrow(Y_ts)

# Manual AIC
AIC_m1 <- -2 * logLik_m1 + 2 * k_m1
AIC_m2 <- -2 * logLik_m2 + 2 * k_m2

# Manual BIC
BIC_m1 <- -2 * logLik_m1 + log(n_obs) * k_m1
BIC_m2 <- -2 * logLik_m2 + log(n_obs) * k_m2

comparison_table <- data.frame(
  Model  = c("Model 1", "Model 2"),
  LogLik = c(logLik_m1, logLik_m2),
  AIC    = c(AIC_m1, AIC_m2),
  BIC    = c(BIC_m1, BIC_m2)
)

write.csv(
  comparison_table,
  file = file.path(processed_data_dir, "model_comparison.csv"),
  row.names = FALSE
)

print(comparison_table)

cat("Model 2 with trigonometric seasonality estimated successfully.\n")

fit_model2$optim.out$convergence

################################################################################

# =========================================================
# FINAL PUBLICATION FIGURE: Banking Indicator vs Real Economy
# =========================================================

# -------------------------------
# 1) Load IPI from external TXT
# -------------------------------
ipi_data_path <- file.path(project_root, "data", "raw", "IPI.txt")
stopifnot(file.exists(ipi_data_path))

ipi_data <- read.csv(
  file   = ipi_data_path,
  header = TRUE,
  dec    = ",",
  sep    = "\t"
)

# Inspect column names
print(names(ipi_data))
print(head(ipi_data))

# -------------------------------
# 2) Define the correct IPI column
# -------------------------------
# Replace 'IPI' with the real column name shown by names(ipi_data)
# Examples:
# IPI_raw <- ipi_data$IPI
# IPI_raw <- ipi_data$Valor
# IPI_raw <- ipi_data$indice

IPI_raw <- ipi_data$IPI

# Safety check
if (is.null(IPI_raw)) {
  stop("ERROR: The selected IPI column does not exist in ipi_data. Check names(ipi_data).")
}

IPI_raw <- as.numeric(IPI_raw)
IPI_raw <- IPI_raw[is.finite(IPI_raw)]

# -------------------------------
# 3) Convert to time series
# -------------------------------
IPI_ts <- ts(
  IPI_raw,
  start = c(1996, 12),
  frequency = 12
)

print(length(IPI_ts))
print(start(IPI_ts))
print(end(IPI_ts))

# -------------------------------
# 4) Align banking indicator and IPI
# -------------------------------
min_length <- min(length(banking_indicator), length(IPI_ts))

banking_indicator_aligned <- as.numeric(banking_indicator[1:min_length])
IPI_aligned               <- as.numeric(IPI_ts[1:min_length])
time_aligned              <- as.numeric(time(Y_ts)[1:min_length])

# -------------------------------
# 5) Create plotting data frame
# -------------------------------
df_plot <- data.frame(
  time = time_aligned,
  Banking_Indicator = banking_indicator_aligned,
  IPI = IPI_aligned
)

# Standardize for visual comparison
df_plot$Banking_Indicator <- as.numeric(scale(df_plot$Banking_Indicator))
df_plot$IPI               <- as.numeric(scale(df_plot$IPI))

# -------------------------------
# 6) Publication plot
# -------------------------------
library(ggplot2)

p_final <- ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = Banking_Indicator, color = "Banking Indicator"), linewidth = 1.2) +
  geom_line(aes(y = IPI, color = "Industrial Production Index"), linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Banking Indicator" = "#1f77b4",
      "Industrial Production Index" = "#d62728"
    )
  ) +
  labs(
    title = "Synthetic Banking Indicator vs Real Economic Activity",
    x = "Time",
    y = "Standardized values",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_final)

ggsave(
  filename = file.path(figures_dir, "final_indicator_plot.png"),
  plot = p_final,
  width = 10,
  height = 6,
  dpi = 300
)

# -------------------------------
# 7) Cross-correlation
# -------------------------------
png(
  filename = file.path(figures_dir, "ccf_banking_indicator_vs_ipi.png"),
  width = 1400,
  height = 900,
  res = 150
)

ccf(
  df_plot$Banking_Indicator,
  df_plot$IPI,
  lag.max = 12,
  main = "Cross-correlation: Banking Indicator vs IPI"
)

dev.off()

# -------------------------------
# 8) Simple validation regression
# -------------------------------
validation_model <- lm(IPI ~ Banking_Indicator, data = df_plot)
print(summary(validation_model))

################################################################################

# =========================================================
# VALIDATION BLOCK: Banking Indicator vs Real Economy
# =========================================================

# -------------------------------
# 1) Required packages
# -------------------------------
extra_packages <- c("ggplot2", "mFilter")

to_install_extra <- extra_packages[!extra_packages %in% rownames(installed.packages())]
if (length(to_install_extra) > 0) install.packages(to_install_extra)

invisible(lapply(extra_packages, library, character.only = TRUE))

# -------------------------------
# 2) Smooth IPI to extract trend
# -------------------------------
# For monthly data, lambda = 14400 is a standard choice in HP filtering
hp_ipi <- mFilter::hpfilter(df_plot$IPI, freq = 14400, type = "lambda")

df_plot$IPI_trend <- as.numeric(hp_ipi$trend)
df_plot$IPI_cycle <- as.numeric(hp_ipi$cycle)

# -------------------------------
# 3) Publication-quality clean figure
# -------------------------------
p_clean <- ggplot(df_plot, aes(x = time)) +
  geom_line(
    aes(y = Banking_Indicator, color = "Banking Indicator"),
    linewidth = 1.2
  ) +
  geom_line(
    aes(y = IPI_trend, color = "IPI Trend"),
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Banking Indicator" = "#1f77b4",
      "IPI Trend" = "#d62728"
    )
  ) +
  labs(
    title = "Synthetic Banking Indicator vs Real Economic Activity",
    subtitle = "Comparison with the smoothed Industrial Production Index",
    x = "Time",
    y = "Standardized values",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

print(p_clean)

ggsave(
  filename = file.path(figures_dir, "publication_clean_plot.png"),
  plot = p_clean,
  width = 10,
  height = 6,
  dpi = 300
)

# -------------------------------
# 4) Cross-correlation analysis
# -------------------------------
png(
  filename = file.path(figures_dir, "ccf_banking_indicator_vs_ipi.png"),
  width = 1400,
  height = 900,
  res = 150
)

ccf_result <- ccf(
  df_plot$Banking_Indicator,
  df_plot$IPI,
  lag.max = 12,
  plot = TRUE,
  main = "Cross-correlation: Banking Indicator vs IPI"
)

dev.off()

# Save CCF values to CSV
ccf_table <- data.frame(
  lag = as.numeric(ccf_result$lag),
  ccf = as.numeric(ccf_result$acf)
)

write.csv(
  ccf_table,
  file = file.path(processed_data_dir, "ccf_banking_indicator_vs_ipi.csv"),
  row.names = FALSE
)

print(ccf_table)

# -------------------------------
# 5) Simple correlation analysis
# -------------------------------
cor_raw <- cor(
  df_plot$Banking_Indicator,
  df_plot$IPI,
  use = "complete.obs"
)

cor_trend <- cor(
  df_plot$Banking_Indicator,
  df_plot$IPI_trend,
  use = "complete.obs"
)

cor_table <- data.frame(
  measure = c("Correlation with raw IPI", "Correlation with IPI trend"),
  value = c(cor_raw, cor_trend)
)

write.csv(
  cor_table,
  file = file.path(processed_data_dir, "correlation_results.csv"),
  row.names = FALSE
)

print(cor_table)

cat("Correlation (Banking Indicator vs raw IPI):", cor_raw, "\n")
cat("Correlation (Banking Indicator vs IPI trend):", cor_trend, "\n")

# -------------------------------
# 6) Regression model: raw IPI
# -------------------------------
model_validation_raw <- lm(IPI ~ Banking_Indicator, data = df_plot)
summary_raw <- summary(model_validation_raw)

print(summary_raw)

sink(file.path(processed_data_dir, "regression_raw_ipi.txt"))
cat("Regression: Raw IPI on Banking Indicator\n\n")
print(summary_raw)
sink()

# -------------------------------
# 7) Regression model: IPI trend
# -------------------------------
model_validation_trend <- lm(IPI_trend ~ Banking_Indicator, data = df_plot)
summary_trend <- summary(model_validation_trend)

print(summary_trend)

sink(file.path(processed_data_dir, "regression_ipi_trend.txt"))
cat("Regression: IPI Trend on Banking Indicator\n\n")
print(summary_trend)
sink()

# -------------------------------
# 8) Scatter plot with fitted line
# -------------------------------
p_scatter <- ggplot(df_plot, aes(x = Banking_Indicator, y = IPI_trend)) +
  geom_point(alpha = 0.6, color = "#1f77b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#d62728", linewidth = 1) +
  labs(
    title = "Relationship between the Banking Indicator and IPI Trend",
    x = "Synthetic Banking Indicator",
    y = "Smoothed Industrial Production Index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_scatter)

ggsave(
  filename = file.path(figures_dir, "scatter_banking_indicator_vs_ipi_trend.png"),
  plot = p_scatter,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------
# 9) Save a compact summary table
# -------------------------------
validation_summary <- data.frame(
  metric = c(
    "Correlation with raw IPI",
    "Correlation with IPI trend",
    "R-squared: regression on raw IPI",
    "R-squared: regression on IPI trend",
    "Adj. R-squared: regression on raw IPI",
    "Adj. R-squared: regression on IPI trend"
  ),
  value = c(
    cor_raw,
    cor_trend,
    summary_raw$r.squared,
    summary_trend$r.squared,
    summary_raw$adj.r.squared,
    summary_trend$adj.r.squared
  )
)

write.csv(
  validation_summary,
  file = file.path(processed_data_dir, "validation_summary.csv"),
  row.names = FALSE
)

print(validation_summary)

cat("Validation block completed successfully.\n")

################################################################################

# Cross-correlation interpretación formal

which.max(ccf_result$acf)
ccf_result$lag[which.max(ccf_result$acf)]

# Granger causality (clave)

library(lmtest)

grangertest(IPI ~ Banking_Indicator, order = 2, data = df_plot)
grangertest(Banking_Indicator ~ IPI, order = 2, data = df_plot)

# Modelo dinámico (muy potente)

dyn_model <- lm(IPI ~ lag(Banking_Indicator, 1) + lag(Banking_Indicator, 2), data = df_plot)
summary(dyn_model)

################################################################################

# =========================================================
# EXTENDED VALIDATION BLOCK
# - Granger causality from lag 1 to 12
# - Regressions using IPI differences and growth rates
# - Comparison against IPI trend
# =========================================================

# -------------------------------
# 1) Required package
# -------------------------------
if (!requireNamespace("lmtest", quietly = TRUE)) install.packages("lmtest")
library(lmtest)

# -------------------------------
# 2) Build a clean aligned dataset
# -------------------------------
# Assumes these objects already exist:
# banking_indicator_aligned
# IPI_aligned
# time_aligned

df_validation <- data.frame(
  time = as.numeric(time_aligned),
  Banking_Indicator = as.numeric(banking_indicator_aligned),
  IPI = as.numeric(IPI_aligned)
)

# Remove missing / non-finite values
df_validation <- df_validation[
  is.finite(df_validation$Banking_Indicator) &
    is.finite(df_validation$IPI),
]

# -------------------------------
# 3) Create IPI trend
# -------------------------------
# A simple smooth trend for publication use
# You can later replace this with HP filter or another method if desired

ipi_trend_fit <- stats::lowess(
  x = df_validation$time,
  y = df_validation$IPI,
  f = 0.20
)

df_validation$IPI_trend <- ipi_trend_fit$y

# -------------------------------
# 4) Correlations with raw IPI and trend
# -------------------------------
cor_raw   <- cor(df_validation$Banking_Indicator, df_validation$IPI, use = "complete.obs")
cor_trend <- cor(df_validation$Banking_Indicator, df_validation$IPI_trend, use = "complete.obs")

# -------------------------------
# 5) Regressions: raw IPI and IPI trend
# -------------------------------
model_raw   <- lm(IPI ~ Banking_Indicator, data = df_validation)
model_trend <- lm(IPI_trend ~ Banking_Indicator, data = df_validation)

summary_raw   <- summary(model_raw)
summary_trend <- summary(model_trend)

# -------------------------------
# 6) Granger causality: lags 1 to 12
# -------------------------------
granger_table <- lapply(1:12, function(L) {
  
  test_banking_to_ipi <- tryCatch(
    grangertest(IPI ~ Banking_Indicator, order = L, data = df_validation),
    error = function(e) NULL
  )
  
  test_ipi_to_banking <- tryCatch(
    grangertest(Banking_Indicator ~ IPI, order = L, data = df_validation),
    error = function(e) NULL
  )
  
  p_banking_to_ipi <- if (!is.null(test_banking_to_ipi)) {
    as.numeric(test_banking_to_ipi$`Pr(>F)`[2])
  } else {
    NA_real_
  }
  
  p_ipi_to_banking <- if (!is.null(test_ipi_to_banking)) {
    as.numeric(test_ipi_to_banking$`Pr(>F)`[2])
  } else {
    NA_real_
  }
  
  data.frame(
    lag_order = L,
    p_banking_causes_ipi = p_banking_to_ipi,
    p_ipi_causes_banking = p_ipi_to_banking
  )
})

granger_table <- do.call(rbind, granger_table)

write.csv(
  granger_table,
  file = file.path(processed_data_dir, "granger_lags_1_to_12.csv"),
  row.names = FALSE
)

print(granger_table)

# -------------------------------
# 7) Cross-correlation function
# -------------------------------
ccf_result <- ccf(
  df_validation$Banking_Indicator,
  df_validation$IPI,
  lag.max = 12,
  plot = FALSE
)

ccf_table <- data.frame(
  lag = as.numeric(ccf_result$lag),
  ccf = as.numeric(ccf_result$acf)
)

write.csv(
  ccf_table,
  file = file.path(processed_data_dir, "cross_correlation_table.csv"),
  row.names = FALSE
)

print(ccf_table)

# Best lag according to maximum correlation
best_ccf_index <- which.max(ccf_table$ccf)
best_ccf_lag   <- ccf_table$lag[best_ccf_index]
best_ccf_value <- ccf_table$ccf[best_ccf_index]

cat("Best CCF lag:", best_ccf_lag, "\n")
cat("Best CCF value:", best_ccf_value, "\n")

# -------------------------------
# 8) IPI differences and growth rates
# -------------------------------
# First difference
df_validation$IPI_diff <- c(NA, diff(df_validation$IPI))

# Approximate growth rate (%)
df_validation$IPI_growth <- c(NA, diff(log(df_validation$IPI + abs(min(df_validation$IPI, na.rm = TRUE)) + 1)) * 100)

# Banking indicator first difference
df_validation$Banking_diff <- c(NA, diff(df_validation$Banking_Indicator))

# Remove NAs for differenced/growth analysis
df_diff <- df_validation[complete.cases(df_validation[, c("IPI_diff", "IPI_growth", "Banking_Indicator", "Banking_diff")]), ]

# -------------------------------
# 9) Regressions with differences and growth
# -------------------------------
# Contemporaneous regressions
model_ipi_diff   <- lm(IPI_diff ~ Banking_Indicator, data = df_diff)
model_ipi_growth <- lm(IPI_growth ~ Banking_Indicator, data = df_diff)

# Dynamic regressions with lags of banking indicator
df_diff$Banking_lag1 <- dplyr::lag(df_diff$Banking_Indicator, 1)
df_diff$Banking_lag2 <- dplyr::lag(df_diff$Banking_Indicator, 2)

df_diff_dyn <- df_diff[complete.cases(df_diff[, c("IPI_diff", "IPI_growth", "Banking_lag1", "Banking_lag2")]), ]

model_ipi_diff_dyn <- lm(IPI_diff ~ Banking_lag1 + Banking_lag2, data = df_diff_dyn)
model_ipi_growth_dyn <- lm(IPI_growth ~ Banking_lag1 + Banking_lag2, data = df_diff_dyn)

summary_ipi_diff       <- summary(model_ipi_diff)
summary_ipi_growth     <- summary(model_ipi_growth)
summary_ipi_diff_dyn   <- summary(model_ipi_diff_dyn)
summary_ipi_growth_dyn <- summary(model_ipi_growth_dyn)

# -------------------------------
# 10) Save regression summary table
# -------------------------------
validation_summary <- data.frame(
  metric = c(
    "Correlation with raw IPI",
    "Correlation with IPI trend",
    "R-squared: regression on raw IPI",
    "Adj. R-squared: regression on raw IPI",
    "R-squared: regression on IPI trend",
    "Adj. R-squared: regression on IPI trend",
    "R-squared: regression on IPI difference",
    "Adj. R-squared: regression on IPI difference",
    "R-squared: regression on IPI growth",
    "Adj. R-squared: regression on IPI growth",
    "R-squared: dynamic regression on IPI difference",
    "Adj. R-squared: dynamic regression on IPI difference",
    "R-squared: dynamic regression on IPI growth",
    "Adj. R-squared: dynamic regression on IPI growth",
    "Best CCF lag",
    "Best CCF value"
  ),
  value = c(
    cor_raw,
    cor_trend,
    summary_raw$r.squared,
    summary_raw$adj.r.squared,
    summary_trend$r.squared,
    summary_trend$adj.r.squared,
    summary_ipi_diff$r.squared,
    summary_ipi_diff$adj.r.squared,
    summary_ipi_growth$r.squared,
    summary_ipi_growth$adj.r.squared,
    summary_ipi_diff_dyn$r.squared,
    summary_ipi_diff_dyn$adj.r.squared,
    summary_ipi_growth_dyn$r.squared,
    summary_ipi_growth_dyn$adj.r.squared,
    best_ccf_lag,
    best_ccf_value
  )
)

write.csv(
  validation_summary,
  file = file.path(processed_data_dir, "validation_summary_extended.csv"),
  row.names = FALSE
)

print(validation_summary)

# -------------------------------
# 11) Publication figure: banking indicator vs IPI trend
# -------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Standardize both for visual comparison
df_plot_trend <- df_validation
df_plot_trend$Banking_std <- as.numeric(scale(df_plot_trend$Banking_Indicator))
df_plot_trend$IPI_trend_std <- as.numeric(scale(df_plot_trend$IPI_trend))

p_trend <- ggplot(df_plot_trend, aes(x = time)) +
  geom_line(aes(y = Banking_std, color = "Banking Indicator"), linewidth = 1.2) +
  geom_line(aes(y = IPI_trend_std, color = "IPI Trend"), linewidth = 1.2, linetype = 2) +
  scale_color_manual(values = c("Banking Indicator" = "#1f77b4", "IPI Trend" = "#d62728")) +
  labs(
    title = "Synthetic Banking Indicator vs Smoothed Industrial Production Trend",
    x = "Time",
    y = "Standardized values",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_trend)

ggsave(
  filename = file.path(figures_dir, "banking_indicator_vs_ipi_trend_publication.png"),
  plot = p_trend,
  width = 10,
  height = 6,
  dpi = 300
)

# -------------------------------
# 12) Publication figure: scatter against IPI trend
# -------------------------------
p_scatter <- ggplot(df_plot_trend, aes(x = Banking_std, y = IPI_trend_std)) +
  geom_point(alpha = 0.7, size = 2, color = "#1f77b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#d62728", linewidth = 1) +
  labs(
    title = "Relationship between the Banking Indicator and IPI Trend",
    x = "Synthetic Banking Indicator",
    y = "Smoothed Industrial Production Index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_scatter)


ggsave(
  filename = file.path(figures_dir, "scatter_banking_indicator_vs_ipi_trend.png"),
  plot = p_scatter,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------
# 13) Save detailed model outputs
# -------------------------------
capture.output(summary_raw,
               file = file.path(processed_data_dir, "summary_regression_raw_ipi.txt"))

capture.output(summary_trend,
               file = file.path(processed_data_dir, "summary_regression_ipi_trend.txt"))

capture.output(summary_ipi_diff,
               file = file.path(processed_data_dir, "summary_regression_ipi_difference.txt"))

capture.output(summary_ipi_growth,
               file = file.path(processed_data_dir, "summary_regression_ipi_growth.txt"))

capture.output(summary_ipi_diff_dyn,
               file = file.path(processed_data_dir, "summary_dynamic_regression_ipi_difference.txt"))

capture.output(summary_ipi_growth_dyn,
               file = file.path(processed_data_dir, "summary_dynamic_regression_ipi_growth.txt"))

cat("Extended validation block completed successfully.\n")

################################################################################

# HP FILTER → separar ciclo y tendencia

if (!requireNamespace("mFilter", quietly = TRUE)) install.packages("mFilter")
library(mFilter)

hp <- hpfilter(df_validation$IPI, freq = 129600)  # mensual

df_validation$IPI_cycle <- hp$cycle
df_validation$IPI_trend_hp <- hp$trend

# mismo para banking
hp_banking <- hpfilter(df_validation$Banking_Indicator, freq = 129600)

df_validation$Banking_cycle <- hp_banking$cycle
df_validation$Banking_trend_hp <- hp_banking$trend

cor(df_validation$Banking_cycle, df_validation$IPI_cycle)
cor(df_validation$Banking_trend_hp, df_validation$IPI_trend_hp)

# VAR MODEL (nivel Q1 real)

if (!requireNamespace("vars", quietly = TRUE)) install.packages("vars")
library(vars)

var_data <- na.omit(df_validation[, c("IPI", "Banking_Indicator")])

var_model <- VAR(var_data, p = 2, type = "const")

summary(var_model)

# Impulse Response (esto es oro)

irf_result <- irf(var_model,
                  impulse = "Banking_Indicator",
                  response = "IPI",
                  n.ahead = 12,
                  boot = TRUE)

plot(irf_result)

# Forecast evaluation (esto te diferencia)

# naive model
model_naive <- lm(IPI ~ 1, data = df_validation)

# model with banking
model_pred <- lm(IPI ~ lag(Banking_Indicator, 1), data = df_validation)

AIC(model_naive, model_pred)
BIC(model_naive, model_pred)
