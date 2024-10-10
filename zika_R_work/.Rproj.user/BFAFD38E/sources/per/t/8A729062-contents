# set current folder to be the main folder 
#install.packages("readxl")
library(readxl)
library(EpiEstim)
# Define the file path
file_path <- "../data/dengue.xlsx"

# Read the Excel file
dengue_data <- read_excel(file_path)
names(dengue_data)


incidences_5_year_sum <- c()
for (col_name in c(paste(2014:2023))) {
  # Calculate the column sum
  column_sum <- sum(dengue_data[[col_name]], na.rm = TRUE)
  # Store the column sum in the incidences_5_year vector
  incidences_5_year_sum <- c(incidences_5_year_sum, column_sum)
}

# Output the results
incidences_5_year_sum

mean(incidences_5_year_sum)

incidences_5_year_sum * 0.01

incidences_5_year_sum * 0.05


zika_cases_list <- list()

# Define misdiagnosis rates

incidences_10_year_sum <- c()

# Iterate over columns "2014" to "2023" and extract non-missing values
for (col_name in as.character(2014:2023)) {
  # Calculate the column sum
  column_sum <- sum(dengue_data[[col_name]], na.rm = TRUE)
  # Store the column sum in the incidences_10_year_sum vector
  incidences_10_year_sum <- c(incidences_10_year_sum, column_sum)
}

# Initialize a list to store results of 50 simulations
zika_cases_list <- list()

# Run the simulation 50 times
for (i in 1:50) {
  # Draw 10 random misdiagnosis rates from a uniform distribution between 0.01 and 0.05
  misdiagnosis_rates <- runif(10, min = 0.01, max = 0.05)
  # Initialize an empty list to store the Zika cases for each rate
  zika_cases_for_this_simulation <- list()
  # Calculate Zika cases for each misdiagnosis rate
  for (rate in misdiagnosis_rates) {
    zika_cases <- incidences_10_year_sum * rate
    zika_cases_for_this_simulation <- append(zika_cases_for_this_simulation, list(zika_cases))
  }
  # Combine the cases for this simulation into a matrix and store it in the main list
  zika_cases_list[[i]] <- do.call(cbind, zika_cases_for_this_simulation)
}

# Combine results into a single matrix by stacking all simulations
zika_cases_matrix <- do.call(cbind, zika_cases_list)

# Calculate median, 25th percentile, and 75th percentile across the 50 simulations
zika_median <- apply(zika_cases_matrix, 1, median)
zika_25th <- apply(zika_cases_matrix, 1, quantile, probs = 0.25)
zika_75th <- apply(zika_cases_matrix, 1, quantile, probs = 0.75)

# Combine results into a data frame
zika_statistics <- data.frame(
  Year = 2014:2023,
  Median = zika_median,
  `25th Percentile` = zika_25th,
  `75th Percentile` = zika_75th
)

zika_median <- median(zika_cases_matrix)
zika_25th <- quantile(zika_cases_matrix, probs = 0.25)
zika_75th <- quantile(zika_cases_matrix, probs = 0.75)

# Print the results
cat("Median of Zika cases:", zika_median, "\n")
cat("25th percentile of Zika cases:", zika_25th, "\n")
cat("75th percentile of Zika cases:", zika_75th, "\n")

print(zika_statistics)



incidences_5_year <- c()

# Iterate over columns "2019" to "2023" and extract non-missing values
for (col_name in c(paste(2015:2023))) {
  # Extract non-missing values from the column
  column_values <- dengue_data[[col_name]][!is.na(dengue_data[[col_name]])]
  # Append the non-missing values to the incidences_5_year vector
  incidences_5_year <- c(incidences_5_year, column_values)
}


# View the incidences_5_year vector
print(incidences_5_year)

# SI 15 - 17 DAYS 
vector_length <- 18
SI_vector <- rep(0, vector_length)
days_to_change <- c(16, 17, 19)
SI_vector[days_to_change] <- 1/3
print(SI_vector)
summary(res_weekly)


config <- make_config(list(si_distr = SI_vector))

method <- "non_parametric_si"

dates = length(incidences_5_year)


res_weekly <- estimate_R_agg(incid = incidences_5_year, 
                             dt = 7L, # aggregation window of the data
                             dt_out = 20L, # desired sliding window length
                             iter = 10L,
                             config = config,
                             method = method,
                             grid = list(precision = 0.01, min = -1, max = 1))

plot(res_weekly)

# Generate data using the provided NLS fit with y-center restriction
Rt <- res_weekly$R$`Mean(R)`
days <- 1:length(Rt)
fit <- nls(Rt ~ 1. + b * sin(2 * pi * days / 365 + phi), 
           start = list( b = 0.2, phi = 0))

summary(fit)


# Initialize incidence vector
incidence <- numeric(length(days))
incidence[1] <- 1  # Initial incidence
Rt <- predict(fit)
# Generate incidence data using the iterative formula I[t+1] = I[t] * Rt
for (t in 1:(length(days) - 1)) {
  incidence[t + 1] <- incidence[t] * Rt[t]
}
# Plot incidence
plot(days, incidence, type = "l", xlab = "Days", ylab = "Incidence",
     main = "Incidence over Time")


phi = 0.027881
# Find the date corresponding to the peak by solving for the maximum value of the sine function
peak_date <- -phi * 365 / (2 * pi)  + 365/2

peak_date



mean(Rt)
sd(Rt)

# Initialize incidence vector
incidence <- numeric(length(days))
incidence
length(incidence)
incidence[1] <- 10  # Initial incidence

# Generate incidence data using the iterative formula I[t+1] = I[t] * Rt
for (t in 1:(length(days) - 1)) {
  incidence[t + 1] <- incidence[t] * Rt[t]
}

# Plot incidence
plot(1:365, incidence[1:365], type = "l", xlab = "Days", ylab = "Incidence",
     main = "Incidence over Time")
abline(v = 176 , col = "blue", lty = 2) 
abline(v = 290 , col = "blue", lty = 2) 

365/4
176/30
86 /30
266 /30
# estimate Zika Rt
# Example vector with some values
length_of_vector <- length(incidences_5_year)
# Generate a new vector with uniform distribution between 0.01 and 0.05, with the same length
random_percentage <- runif(length_of_vector, min = 0.01, max = 0.01)
sym_cases <- random_percentage * incidences_5_year
total_cases <- sym_cases * 5


# zika SI
si_distr  <-c(0,0,0.1,0,0.1,0.2,0.1,0.2,0,0,0,0,0.1,0.2)
config <- make_config(list(si_distr = si_distr))
method <- "non_parametric_si"


res_weekly <- estimate_R_agg(incid = total_cases, 
                             dt = 7L, # aggregation window of the data
                             dt_out = 7L, # desired sliding window length
                             iter = 10L,
                             config = config,
                             method = method,
                             grid = list(precision = 0.01, min = -1, max = 1))

plot(res_weekly)
daily_cases = total_cases/7
median(total_cases)/7



length(res_weekly$R$`Mean(R)` )
length(res_weekly$I_local)



# add 13 zero to Rt, now Rt and I_local of 3654
# 3654 is for year 2014 to 2023; simply put 365 for a year and ignore the values from 3651 to 3654
# do two plots, one is days vs Rt, one is days and I_local
# however, for the x axis, we wish to only show  year 2014, year 2015... till year 2023
# arrange two plots vertically, and save it to pdf, name it to be "mis_rate_0.01-0.05"
Rt = res_weekly$R$`Mean(R)`
Rt <- c(rep(0, 13),Rt)
I_local = res_weekly$I_local

Rt = Rt[1:3650]
I_local = I_local[1:3650]

years <- seq(2014, 2023, by = 1)

# Creating plots
#pdf("mis_rate_0.01-0.05.pdf")
pdf("mis_rate_0.01.pdf")
par(mfrow = c(2, 1))  # 2 plots arranged vertically
plot(1:length(Rt), Rt, type = "l", xaxt = "n", xlab = "Days", ylab = "Rt", main = "Rt Over Time")
abline(h = 1.0, col = "red", lty = 2)
axis(1, at = seq(1, length(Rt), by = 365), labels = years)
plot(1:length(I_local), I_local, type = "l", xaxt = "n", xlab = "Days", ylab = "I_local", main = "Incidences Over Time")
axis(1, at = seq(1, length(I_local), by = 365), labels = years)
dev.off()

quartiles <- quantile(daily_cases, probs = c(0.25, 0.75))

# Extract Q1 and Q3 from the result
Q1 <- quartiles[1]
Q3 <- quartiles[2]
Q1
Q3

#plot(1:length(total_cases),total_cases)

barplot(total_cases/7, 
        main = "Total Cases",  # Title of the plot
        xlab = "Time",         # Label for x-axis
        ylab = "Total Cases",  # Label for y-axis
        col = "skyblue",       # Bar color
        border = "black",      # Border color
        ylim = c(0, max(total_cases)/7 * 1.1),  # Adjust ylim to give some space above the highest bar
    
)

hist(total_cases/7, 
     main = "Distribution of Total Cases",  # Title of the plot
     xlab = "Total Cases",                 # Label for x-axis
     ylab = "Frequency",                   # Label for y-axis
     col = "skyblue",                      # Color of bars
     border = "black"                      # Color of border
)


Rt = res_weekly$R$`Mean(R)`

Rt




plot(Rt, type = "l")
mean(Rt)
sd(Rt)
1 + 1.96*0.4
1 - 1.96*0.4

quartiles <- quantile(Rt, probs = c(0.1, 0.9))

# Extract Q1 and Q3 from the result
Q1 <- quartiles[1]
Q3 <- quartiles[2]
Q1
Q3

# wrong method 
hist(Rt, main = "Histogram of Data", xlab = "Value", ylab = "Frequency", col = "lightblue")

shapiro_test <- shapiro.test(Rt)
print(shapiro_test)
install.packages("tseries")
library(tseries)


# Assume your data is stored in a variable called 'data'
# Perform the Augmented Dickey-Fuller test



# Print the results
print(adf_test)


mean(Rt) # 1.03 starionary Rt 
install.packages("forecast")
library(forecast)

# Assume 'Rt' is your time series data

# Plot autocorrelation
acf(Rt, main = "Autocorrelation Plot for Rt Series")


ts_data <- ts(Rt, frequency = 1)
arima_model <- auto.arima(ts_data)
summary(arima_model)
# AIC=-18659.51   AICc=-18659.48   BIC=-18616.12
arima_model$s.e.
coefficients <- coef(arima_model)
standard_errors <- arima_model$var.coef

# daily Rt 
# Mean parameter
mean_estimate <- coefficients["mean"]

# Standard deviation (noise)
sigma_squared <- coefficients["sigma^2"]
noise_sd <- sqrt(sigma_squared)


# Fit ARIMA(0,0,0) model
arima_model <- Arima(Rt, order = c(0, 0, 0))

# Summary of the fitted model
summary(arima_model)



# Estimated mean and standard error
mean_estimate <- 1.0356
se_estimate <- 0.26

# Number of time points in the time series
num_points <- length(Rt)

# Generate white noise
white_noise <- rnorm(num_points, mean = mean_estimate, sd = se_estimate)

# Plot the white noise time series
plot(white_noise, type = "l", main = "White Noise Time Series", xlab = "Time", ylab = "Value")



len_Rt <- length(Rt)
print(len_Rt)

# Create a sequence of x values
x <- seq(1, length(Rt), by = 1)



# Plot the original data and the fitted sine function
plot(x, Rt, type = "l", col = "blue", lwd = 2, main = "Fitted Sine Function", xlab = "x", ylab = "Rt")
lines(x, white_noise, col = "red", lwd = 2)
legend("topleft", legend = c("Original Data", "Fitted Sine Function"), col = c("blue", "red"), lwd = 2, cex = 0.8)


MAPE <- mean(abs((1.0356 - Rt) / Rt) * 100)    # naive prediction 21% error 



Rt_ts = Rt

train_size <- round(length(Rt_ts) * 0.8)

# Split the time series into training and test sets
train_data <- window(Rt_ts, end = train_size)
test_data <- window(Rt_ts, start = train_size + 1)

# Calculate mean and standard error from the training set
train_mean <- mean(train_data)
train_sd <- sd(train_data)

# Generate the estimated time series with noise
estimated_ts <- rnorm(length(test_data), mean = train_mean, sd = train_sd)

# Calculate the mean absolute percentage error (MAPE)
MAPE <- mean(abs((test_data - estimated_ts) / test_data) * 100)

x = 1:length(test_data)
plot(x, test_data, type = "l", col = "blue", lwd = 2, main = "Fitted ", xlab = "x", ylab = "Rt")
lines(x, estimated_ts, col = "red", lwd = 2)


# Print the MAPE
print(paste("Mean Absolute Percentage Error (MAPE):", MAPE)). # 28% error


arima_model <- Arima(Rt, order = c(0, 0, 0))
summary(arima_model)
forecast_values <- forecast(arima_model, h = length(test_data))

result <- sqrt(0.07223)
forecast_ts <- forecast_values$mean
white_noise <- rnorm(length(test_data), mean = 0, sd = result)
forecast_ts <- forecast_ts + white_noise

plot(x, test_data, type = "l", col = "blue", lwd = 2, main = "Fitted ", xlab = "x", ylab = "Rt")
lines(x, forecast_ts, col = "red", lwd = 2)

# modelling the seasonality
train_data_periodic <- ts(train_data, frequency = 365)

plot(train_data_periodic, type = "l", col = "blue", lwd = 2, main = "Fitted ", xlab = "x", ylab = "Rt")
harmonics <- fourier(train_data_periodic, K = 10)
# 
fit <- auto.arima(train_data_periodic, xreg = harmonics, seasonal = FALSE)
newharmonics <- fourier(train_data_periodic, K = 10, h = length(test_data))
fc <- forecast(fit, xreg = newharmonics)
autoplot(fc) + gglines

# Plotting the mean forecast
plot(fc$mean, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Values")

# Adding test data to the plot
lines(test_data, col = "red", lwd = 2)

# Adding a legend
legend("topright", legend = c("Mean Forecast", "Test Data"), col = c("blue", "red"), lty = 1, lwd = 2)
# Plotting the mean forecast

x = 1: length(test_data)
plot(x,fc$mean, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Values", ylim = c(0.5,2.0))

# Adding test data to the plot
lines(x,as.numeric(test_data), col = "red", lwd = 2)

# Adding a legend
legend("topright", legend = c("Mean Forecast", "Test Data"), col = c("blue", "red"), lty = 1, lwd = 2)


random_walk_model <- rwf(train_data)

# Summary of the fitted random walk model
summary(random_walk_model)

plot(train_data, type = "l", col = "blue", lwd = 2, main = "Fitted Random Walk Model", xlab = "Time", ylab = "Value")

# Add the fitted random walk model to the plot
lines(random_walk_model$fitted, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Original Data", "Fitted Random Walk Model"), col = c("blue", "red"), lwd = 2, cex = 0.8)



# Make predictions for the test_data using the random walk model
forecast_values <- forecast(random_walk_model, h = length(test_data))

# Extract the point forecasts from the forecast object
forecast_ts <- forecast_values$mean

# Plot the original time series data and the forecasted values
plot(test_data, type = "l", col = "blue", lwd = 2, main = "Random Walk Model Forecast", xlab = "Time", ylab = "Value")
lines(forecast_ts, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Actual Data", "Forecast"), col = c("blue", "red"), lwd = 2, cex = 0.8)


mean_estimate <- random_walk_model$model$mean
noise_sd <- random_walk_model$model$sd



# Set seed for reproducibility
set.seed(123)

# Number of time periods
n <- length(test_data)

# Generate random noise
noise <- rnorm(n, mean = 0, sd = 0.5 )  # Adjust mean and standard deviation as needed

# Initialize vector to store random walk values
random_walk <- numeric(n)

# Set initial value
random_walk[1] <- 1.03  # Initial value can be any starting point

# Generate random walk values
for (i in 2:n) {
  random_walk[i] <- random_walk[i - 1] + noise[i]
}

# Plot the random walk
plot(random_walk, type = "l", main = "Random Walk with Noise", xlab = "Time", ylab = "Value")



# Define the Kalman filter model using train_data
kalman_model <- KalmanForecast(train_data)

# Make predictions for the test_data using the Kalman filter model
forecast_values <- KalmanForecast(train_data, model = kalman_model, newdata = test_data)

# Extract the point forecasts from the forecast object
forecast_ts <- forecast_values$mean

# Plot the actual test data and the forecasted values
plot(test_data, type = "l", col = "blue", lwd = 2, main = "Kalman Filter Forecast", xlab = "Time", ylab = "Value")

lines(forecast_ts, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Actual Data", "Forecast"), col = c("blue", "red"), lwd = 2, cex = 0.8)




plot(test_data, type = "l", col = "blue", lwd = 2, main = "Kalman Filter Forecast", xlab = "Time", ylab = "Value")



b# Apply Fourier transform
train_fft <- fft(train_data)

# Compute the power spectral density
psd <- abs(train_fft)^2

# Frequency range
freq <- (0:(length(train_data) - 1)) / length(train_data)


# Plot the power spectral density
plot(freq, psd, type = "l", xlab = "Frequency", ylab = "Power", main = "Power Spectral Density")

# Find peaks in the power spectral density to identify potential periodic components
peaks <- findpeaks::findpeaks(psd)$peak[,1]  # You may need to install and load the 'findpeaks' package
periods <- 1 / abs(freq[peaks])  # Calculate periods corresponding to the peaks

cat("Periods (in terms of index):", periods, "\n")



autocorr <- Acf(train_data)

# Plot Autocorrelation Function (ACF)
plot(autocorr, main = "Autocorrelation Function (ACF)")


library(stats)


result <- spec.pgram(Rt_ts[1:train_size], log = "no", fast = TRUE, demean = FALSE, plot = FALSE)

# Extract frequency and PSD values
frequencies <- result$freq # 1458, half the sampling rate
PSD <- result$spec

# Plot the Power Spectral Density
plot(frequencies, PSD, type = "l", xlab = "Frequency", ylab = "PSD", main = "Power Spectral Density",xlim = c(0,0.05))

# kernel smoothing
library(astsa)
k = kernel("daniell", 4)
x = ts(Rt_ts[1:train_size],frequency = 365)
class(
mvspec(x, log="no")
mvspec(x, k, log="no")
soi
astsa::mvspec(soi)

max_indices <- which.max(PSD)



# Retrieve corresponding frequencies
highest_freq <- frequencies[max_indices]

period <- 1 / highest_freq # 0.5 year peoridicity
#


plot(train_data, type = "l", xlab = "Time", ylab = "Value", main = "Time Series with Periodic Vertical Lines")

# Add vertical lines at intervals of the identified period
period <- 182
for (i in 1:10) {
  abline(v = period * i, col = "red", lty = 2)
}

train_data_periodic <- ts(train_data, frequency = 182)

# Print the first few observations of the new periodic time series
head(train_data_periodic)
sarima_model <- auto.arima(train_data_periodic, seasonal = TRUE)

# Print the summary of the SARIMA model
summary(sarima_model)


########
# fit sin
time <- 1: length(train_data)

# Define the sine function with parameters: amplitude (A), frequency (f), phase shift (phi), and vertical shift (c)
sine_function <- function(t, A, f, phi, c) {
  A * cos(2 * pi * f * t + phi) + c
}

# Fit the sine function to 'train_data' using nonlinear least squares (nls)
fit <- nls(train_data ~ sine_function(time, A, f, phi, c), 
           start = list(A = 0.1, f = 1/40, phi = 0, c = mean(train_data)))

# Print the summary of the fit
summary(fit)

plot(time, train_data, type = "l", col = "blue", xlab = "Time", ylab = "Value", main = "Fitting with Sine Function")
white_noise = rnorm(length(time), mean = 0, sd = 0.2561)
lines(time, predict(fit) , col = "red")
legend("topright", legend = c("Original Data", "Fitted Sine Function"), col = c("blue", "red"), lty = 1)


plot(time, train_data, type = "l", col = "blue", xlab = "Time", ylab = "Value", main = "Fitting with Sine Function")
abline(h = 0.8, col = "red")
abline(h = 1.2, col = "red")


install.packages("pracma")
library(pracma)

# Find peaks in the time series data
train_data <-as.numeric(train_data)
peaks <- findpeaks(train_data,minpeakheight = 1.2)


peak_period <- colSums(peaks)[4] - colSums(peaks)[3] + nrow(peaks)

peaks[,1]
hist(peaks[, 1], breaks = 20, col = "skyblue", main = "Histogram of Peak Indices", xlab = "Index")



peak_vec <- vector(mode = "numeric", length = peak_period)

w_start <- 1

for(i in seq_len(nrow(peaks))){
  
  len_i <- length(peaks[i, 3] : peaks[i, 4])
  w_end <- w_start + len_i-1
  peak_vec[w_start:w_end] <- c(peaks[i, 3]:peaks[i, 4])
  w_start <- len_i + w_start
  
}


plot(1:peak_period, train_data[peak_vec], type = "l", lty = 1, ylab = "Force", xlab = "Time")

bottoms <- findpeaks(-train_data,minpeakheight = - 0.8)

# Plot the time series data with bottoms
plot(-train_data, type = "l", xlab = "Time", ylab = "Value", main = "Time Series with Bottoms")
points(bottoms[, 1], bottoms[, 2], col = "blue", pch = 19)
legend("topright", legend = "Bottoms", col = "blue", pch = 19)


hist(bottoms[, 1], breaks = 20, col = "skyblue", main = "Histogram of Peak Indices", xlab = "Index")



# Define the parameters
top <- 1.3
bottom <- 0.7
period <- 40
first_trough <- 15

# Calculate amplitude
A <- (top - bottom) / 2

# Calculate vertical shift
c <- (top + bottom) / 2

# Calculate phase shift
b <- -2 * pi / period * (first_trough - 1)

# Define the sine function
sin_function <- function(x) {
  A * sin((2 * pi / period) * x + b) + c
}

# Generate x values
x_values <- seq(1, length.out = length(train_data))

# Generate y values using the sine function
y_values <- sin_function(x_values)

# Plot the sine function
plot(x_values,train_data, type = "l", xlab = "Time", ylab = "Value", main = "Time Series with Bottoms")
lines(x_values, y_values, type = "l", xlab = "x", ylab = "y", main = "Sine Function", col = "red")






lines(time, sine_function_manual(time), col = "red")
legend("topright", legend = c("Original Data", "Fitted Sine Function"), col = c("blue", "red"), lty = 1)



##########
# fft

# Assuming 'train_data' contains your original time series data

train_data <- train_data - mean(train_data)
# Perform FFT on the time series data
train_data_fft <- fft(train_data)

# Calculate the magnitude of FFT coefficients
train_data_fft_magnitude <- Mod(train_data_fft)

# Plot the magnitude spectrum of the FFT coefficients
plot(train_data_fft_magnitude, type = "l", xlab = "Frequency", ylab = "Magnitude", main = "Magnitude Spectrum")

# Filter out high-frequency components (e.g., above a certain threshold)
threshold <- 50  # Adjust as needed
train_data_fft_filtered <- train_data_fft
train_data_fft_filtered[which(train_data_fft_magnitude > threshold)] <- 0

# Inverse FFT to obtain the denoised time series data
train_data_denoised <- Re(fft(train_data_fft_filtered, inverse = TRUE) / length(train_data_fft_filtered))

# Plot the original and denoised time series data
plot(train_data, type = "l", col = "blue", xlab = "Time", ylab = "Value", main = "Original vs. Denoised Time Series")
lines(train_data_denoised , col = "red")
legend("topright", legend = c("Original", "Denoised"), col = c("blue", "red"), lty = 1)


stl_result <- stl(train_data, s.window = "periodic")

# Plot the decomposition components
plot(stl_result)












if (!requireNamespace("signal", quietly = TRUE)) {
  install.packages("signal")
}
library(signal)

# Assuming 'train_data' contains your original time series data

# Perform FFT on the time series data
train_data_fft <- fft(train_data)

# Calculate the magnitude of FFT coefficients
train_data_fft_magnitude <- Mod(train_data_fft)

# Plot the Power Spectral Density
plot(train_data_fft_magnitude, type = "l", xlab = "Frequency", ylab = "Magnitude", main = "Power Spectral Density")

# Define a cutoff frequency to remove high-frequency noise
cutoff_frequency <- 0.05  # Adjust as needed

# Create a frequency domain filter
filter <- butter(10, cutoff_frequency, type = "low", plane = "z")




# Apply the filter to the FFT coefficients
filtered_fft <- filter$a * train_data_fft

# Perform inverse FFT to obtain the filtered time series data
filtered_data <- Re(fft(filtered_fft, inverse = TRUE) / length(train_data_fft))

# Plot the original and filtered time series data
plot(train_data, type = "l", col = "blue", xlab = "Time", ylab = "Value", main = "Original vs. Filtered Time Series")
lines(filtered_data, col = "red")
legend("topright", legend = c("Original", "Filtered"), col = c("blue", "red"), lty = 1)

# fourier filter 
# fit fourier 

