

library(readxl)
library(EpiEstim)
file_path <- "../data/dengue.xlsx"

dengue_data <- read_excel(file_path)


incidences_5_year <- c()

# Iterate over columns "2019" to "2023" and extract non-missing values
for (col_name in c(paste(2014:2023))) {
  # Extract non-missing values from the column
  column_values <- dengue_data[[col_name]][!is.na(dengue_data[[col_name]])]
  # Append the non-missing values to the incidences_5_year vector
  incidences_5_year <- c(incidences_5_year, column_values)
}

# View the incidences_5_year vector
print(incidences_5_year) # 522 weeks 


incidences_5_year

# construct generation interval
vector_length <- 18
SI_vector <- rep(0, vector_length)
days_to_change <- c(16, 17, 18)
SI_vector[days_to_change] <- 1/3
print(SI_vector)


config <- make_config(list(si_distr = SI_vector))

method <- "non_parametric_si"

res_weekly <- estimate_R_agg(incid = incidences_5_year, 
                             dt = 7L, # aggregation window of the data
                             dt_out = 20L, # desired sliding window length
                             iter = 30L,
                             config = config,
                             method = method,
                             grid = list(precision = 0.01, min = -1, max = 1))



plot(res_weekly) # 3634
res_weekly$I # from day 15

# Rt is stationary
Rt = res_weekly$R$`Mean(R)`
plot(Rt, type = "l")

# save Rt 
saveRDS(Rt, "rt_values.rds")
Rt = readRDS("../data/rt_values.rds")


37 # daily 3.7 

#####################################
# clean noise first and then fit
fft_result <- fft(Rt)

# # Calculate the frequency values (normalized frequencies)
# fft_freq <- (0:(n - 1)) / n  # Frequency range from 0 to 0.5
# fft_freq <- fft_freq * n
# 
# power_spectrum <- Mod(fft_result)^2
# 
# plot(fft_freq, power_spectrum, type = "l", col = "blue",
#      xlab = "Frequency", ylab = "Power",
#      main = "Frequency Spectrum (Power Spectrum)", ylim = c(0,20))
# 
# # Highlight the relevant frequencies
# abline(v = 0.1, col = "red", lty = 2)
# 
# 
# 
# fft_freq <- (0:(length(fft_result) - 1)) / length(fft_result) 
# frequency_threshold <- 0.05  # This threshold might require tuning based on the data
# 
# # Create a logical mask to keep frequencies below the threshold
# mask <- fft_freq < frequency_threshold
# 
# # Apply the mask to filter high-frequency noise
# fft_filtered <- fft_result
# fft_filtered[!mask] <- 0 
# 
# # Apply IFFT to get the filtered time series
# filtered_time_series <- Re(fft(fft_filtered, inverse = TRUE)) / length(fft_filtered)  # Normalization
# 
# # Plot the original and filtered time series for comparison
# plot(Rt, type = "l", col = "blue", main = "Original vs Filtered Time Series",
#      xlab = "Time", ylab = "Value")
# lines(filtered_time_series, col = "red", lty = 2)  # Plot the filtered data
# legend("topright", legend = c("Original", "Filtered"), col = c("blue", "red"), lty = c(1, 2))
# 
# 
# clean_Rt <- filtered_time_series
# plot(filtered_time_series, type = "l")
# 
# 

## Create a time series object
set.seed(25)
inds <- seq(as.Date("2014-01-15"), as.Date("2015-10-14"), by = "day")
Rt_ts <- ts(Rt,     # random data
           start = c(2014, as.numeric(format(inds[1], "%j"))),
           frequency = 365)

# Split the time series into training and test sets
train_data <- window(Rt_ts,end = 2021)
#test_data <- window(Rt_ts, frequncy = 365,start = train_size + 1)

test_data = window(Rt_ts, start = 2021)

# 
length(Rt_ts)
length(test_data)
length(train_data)
# frequency(train_data)

# do the forecasting 
library(fpp3)
library(forecast)

acf(train_data) # 5,7
pacf(train_data)
harmonics <- fourier(train_data, K = 5)



# fit our own
#my_model <- Arima(train_data,xreg = harmonics, order = c(3,0,0))
#plot(forecast(my_model,xreg = newharmonics,h=length(test_data)))

# auto fit harmonic regression 
fit <- auto.arima(train_data, xreg = harmonics, seasonal = FALSE)
summary(fit)
newharmonics <- fourier(train_data, K = 5, h = length(test_data))

his_fir <-forecast(fit, xreg = harmonics)
fc <- forecast(fit, xreg = newharmonics)
autoplot(fc) 


newharmonics <- fourier(train_data, K = 5, h = 15 * 365)
fc <- forecast(fit, xreg = newharmonics)
autoplot(fc) 
fc$mean
saveRDS(fc$mean,"../data/fc_mean.rds")
# save fc and fitted 
# saveRDS(Rt, "rt_values.rds")

# # Fit a harmonic regression using order 10 for each type of seasonality
# # monthly and yearly 
# multi_fit <- tslm(train_data ~ fourier(train_data, K = c(10, 10)))
# taylor
# str(taylor)
# #
# msts_Rt <- msts(Rt, seasonal.periods = c(30+ 5/12, 365.25))
# msts_train <- window(msts_Rt, end = 9) 
# msts_test <- window(msts_Rt, start = 10) 
# multi_fit <- tslm(msts_train ~ fourier(msts_train, K = c(10, 10)))
# 
# fc <- forecast(multi_fit, newdata = data.frame(fourier(msts_train, K = c(10, 10), h = length((msts_test)))))
# 
# autoplot(fc)


# # Forecast 20 working days ahead
# fc <- forecast(fit, newdata = data.frame(fourier(taylor, K = c(10, 10), h = 20*48)))
# 
# # Plot the forecasts
# autoplot(fc)

# period Rt 


# 1.03, 1.2, 0.9 0.018

summary(fit)
fc$mean

length(train_data)

#fc = forecast(my_model,xreg = newharmonics,h=length(test_data))

fore_df <- data.frame(Date = (1+ length(train_data)):(length(train_data)+length(fc$mean)) ,
                      actual_Vale = as.numeric(test_data),
                      Predicted_Value = as.numeric(fc$mean),
                      LL = as.numeric(fc$lower),
                      UL = as.numeric(fc$upper)
)

# our period is 52 weeks (52*7 = 364 days)
names(fc)

hist_df <- data.frame(Date = 1:length(train_data) ,
                      actual_Vale = as.numeric(train_data),
                      Predicted_Value = as.numeric(his_fir$mean),
                      LL = as.numeric(his_fir$lower),
                      UL = as.numeric(his_fir$upper)
)

hist_df <- as.data.frame(hist_df)
fore_df <- as.data.frame(fore_df)

# Now try rbind()
his_df <- rbind(hist_df, fore_df)

ggplot(data=his_df, aes(x=Date)) +
  geom_line(aes(y=Predicted_Value, colour = "Predicted Rt\n"), colour = "blue",linewidth=1.) +
  geom_line(aes(y=actual_Vale, colour = "actual Rt"), colour = "black",linewidth=1.) +
  theme_bw() +  
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 12))+
  ylab("Rt") + xlab("days") +
  ylim(c(0,3.))+
  geom_ribbon(aes(ymin=LL, ymax=UL), fill="blue", alpha=0.5) +
  ggtitle(paste0("Estimated Rt with 95% Interval vs Actual Rt"))+
  geom_vline(xintercept = length(train_data), col = "red")



# 8 years - 2 years 
 p <- ggplot(data=his_df, aes(x=Date)) +
  geom_line(aes(y=Predicted_Value, colour = "Predicted Rt"), linewidth=1.) +
  geom_line(aes(y=actual_Vale, colour = "Actual Rt"), linewidth=1.) +
  geom_ribbon(aes(ymin=LL, ymax=UL, fill = "95% interval"), alpha=0.5) +
  geom_vline(xintercept =  length(train_data), col = "red") +
  geom_text(aes(x =  length(train_data), y = 2, label = "Validation Years"), 
            hjust = -0.1, vjust = 0, colour = "black", size = 3.5, angle = 0) +
  geom_text(aes(x = length(train_data), y = 2, label = "Training Years"), 
            hjust = 1.1, vjust = 0, colour = "black", size = 3.5, angle = 0) +
  scale_colour_manual(values = c("Predicted Rt" = "blue", "Actual Rt" = "black")) +
  scale_fill_manual(values = c("95% interval" = "blue")) +
  theme_bw() +  
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 12)) +
  ylab("Rt") + 
  xlab("Days") +
  ylim(c(0, 3.)) +
  ggtitle("Estimated Rt with 95% Interval vs Actual Rt")

 ggsave("estimated_rt_plot.png", plot = p, width = 10, height = 7)
 
# 1.02 
 # ggplot(data=his_df, aes(x=Date)) +
 #   geom_line(aes(y=Predicted_Value, colour = "gengue like"), linewidth=1.)
   #geom_line(aes(y=Predicted_Value - 0.1 , colour = "smaller curve"), linewidth=1.)+
   #geom_line(aes(y=Predicted_Value + 0.1 , colour = "big curve"), linewidth=1.)
   #scale_colour_manual(values = c("Predicted Rt" = "blue"))

 plot(his_df$Predicted_Value, type = "l")
 
 save( his_df$Predicted_Value )
 # assume a baseline of 200
 # 200 * 
 
 his_df$Predicted_Value[70]
 
 

# 
# ggplot(data=his_df, aes(x=Date)) +
#   geom_line(aes(y=Predicted_Value, colour = "Predicted Rt\n"), linewidth=1.2) +
#   geom_line(aes(y=actual_Vale, colour = "actual Rt"), linewidth=1.2) +
#   #geom_line(aes(y=Predicted_Value_2, colour = "Driven by post-covid\n houseprice trend\n"), linewidth=1.2)+
#   theme_bw() +  
#   theme(legend.title = element_blank(),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         title = element_text(size = 12))+
#   ylab("Rt") + xlab("days") +
#   #geom_vline(xintercept=as.Date("2024-01-01"), linetype=2) + 
#   geom_ribbon(aes(ymin=LL.95., ymax=UL.95.), fill="grey", alpha=0.5) +
#   #geom_ribbon(aes(ymin=LL_2, ymax=UL_2), fill="grey", alpha=0.5) +
#   ggtitle(paste0("approximate SPI Trend till 2050,  with 50% interval")) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0,size = 12)) #+
# # geom_text(aes(x = as.Date("2012-01-01"), y = 550, label = "Historical Data"), vjust = -0.5) +
# # geom_text(aes(x = as.Date("2040-01-01"), y = 550, label = "Predicted Data"), vjust = -0.5)
 


