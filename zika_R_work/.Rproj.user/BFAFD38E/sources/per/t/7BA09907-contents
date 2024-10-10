#install.packages("readxl")
library(readxl)
library(EpiEstim)

# Read the Excel file
data <- read_excel("../data/zika.xlsx")

# Convert start_date column to Date format
data$start_date <- as.Date(data$start_date, format = "%d/%m/%y")

# Filter data based on start_date column
filtered_data <- data[data$start_date >= as.Date("2016-08-21") & data$start_date <= as.Date("2016-11-26"), c("start_date", "case")]

# Print the filtered data
colnames(filtered_data) <- c("date", "value")
print(filtered_data)

plot(filtered_data)

weekly_incid = filtered_data$value
weekly_incid <- weekly_incid[-c(1)]
weekly_incid = c(215 ,107 , 62 , 11 ,  8  , 4  ,11  ,11 , 12)
# configure the parameters 
si_distr  <-c(0,0,0.1,0,0.1,0.2,0.1,0.2,0,0,0,0,0.1,0.2)
config <- make_config(list(si_distr = si_distr))
method <- "non_parametric_si"

Rt_weekly <- estimate_R_agg(incid = weekly_incid, 
                            dt = 7L, # aggregation window of the data
                            dt_out = 7L, # desired sliding window length
                            iter = 30L,
                            config = config,
                            method = method,
                            grid = list(precision = 0.01, min = -1, max = 1))

plot(Rt_weekly)


mean(Rt_weekly$R$`Mean(R)`) # 0.84
plot(Rt_weekly$R$`Mean(R)`)

Rt_weekly$R$`Mean(R)`
# Wallinga and Teunis

amplitude <- 0.05
center <- 1
period <- 365
angular_frequency <- 2 * pi / period

# Create a function to calculate y based on x
sin_wave <- function(x) {
  return(amplitude * sin(angular_frequency * x) + center)
}

# Generate x values for a year (365 days)
x_values <- seq(0, 3650, length.out = 1000)

# Calculate y values using the function
y_values <- sapply(x_values, sin_wave)

plot(x_values, y_values, type = "l", col = "grey", lwd = 2,
     xlab = "Days", ylab = "Value",
     main = "Sine Wave with Rt ,Center at 1.0,\n Amplitude of 0.2, and Period of 365")


# 
# Number of days
n_days <- 365*15
w_s <- c(0, 0, 0.1, 0, 0.1, 0.2, 0.1, 0.2, 0, 0, 0, 0, 0.1, 0.2)
R_t <- 3.62

# Vector to store incidences for each day
I <- numeric(n_days)

# Initial incidence on day 1
I[1:14] <- 10
fc_mean = readRDS("../data/fc_mean.rds") 
R_t <- fc_mean -0.1
length(R_t)
# Loop through each day to compute the new incidence
for (t in 15:n_days) {
  # Determine the infectious pressure for the current day
  infectious_pressure <- 0
  
  # Loop to calculate the infectious pressure based on past incidences and weights
  for (s in 1:(t-1)) {
    if (s <= length(w_s)) {
      infectious_pressure <- infectious_pressure + I[t - s] * w_s[s]
    }
  }
  
  # Compute the expected number of new incidences using the Poisson distribution
  #R_t = sin_wave(t)
  #expected_incidences <- 0.9 * infectious_pressure
  k_t = 0.8
  if (t  < 180) {
    k_t = 3.68
  }
  #expected_incidences <- R_t[t] * infectious_pressure
  expected_incidences <- k_t * infectious_pressure
  
  #Draw from the Poisson distribution to determine new incidences
  #I[t] <- expected_incidences
  I[t] <- rpois(1, expected_incidences)
  if (t %% 30 == 0) {
   
    I[t] <- I[t] + runif(1, 1, 3)
  }
    
}

# Print the incidence vector for 180 days
days = 1:n_days
plot(days, I, type = "l", col = "lightblue", lwd = 2,
     xlab = "Day", ylab = "Incidences",
     main = "Incidences over Time (15 years)")
sum(I)/15

I[1:100]
I[100:500]
I[2000:2100]
# assume 100 - 200 
# "../data/fc_mean.rds"

# 



