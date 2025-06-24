
### Working with tide heights ... 

# Skip the first 6 header lines
df <- read.csv( paste(data_dir, "cedar_island_tides_2020to2028.csv", sep="/"), skip = 6, header = TRUE, stringsAsFactors = FALSE)

# Rename columns for clarity if needed
colnames(df) <- c("datetime", "water_level")

# Parse datetime
df$datetime <- as.POSIXct(df$datetime, format = "%Y/%m/%d %H:%M", tz = "Canada/Pacific")


# Full time series ... 
plot(df$datetime, df$water_level, type = "l",
     main = "Cedar Island Water Level (full record)",
     xlab = "Date", ylab = "Water Level (m)")

# Apr to Sep 2025 ... 
df_sub <- subset(df, datetime >= as.POSIXct("2025-04-01", tz = "Canada/Pacific") &
                   datetime < as.POSIXct("2025-10-01", tz = "Canada/Pacific"))

# Just June 2025 ... 
df_sub <- subset(df, datetime >= as.POSIXct("2025-06-01", tz = "Canada/Pacific") &
                    datetime < as.POSIXct("2025-07-01", tz = "Canada/Pacific"))


plot(df_sub$datetime, df_sub$water_level, type = "l",
     main = "Cedar Island Water Level (June 2025)",
#     main = "Cedar Island Water Level (Apr–Sep 2025)",
     xlab = "Date", ylab = "Water Level (m)")


### HH LL plot

#uses df_sub, from above

library(zoo)
library(dplyr)

# Add logical columns for local max/min
df_sub$high <- rollapply(df_sub$water_level, width=3, 
                          FUN=function(x) x[2] > x[1] & x[2] > x[3], 
                          fill=NA, align="center")
df_sub$low <- rollapply(df_sub$water_level, width=3, 
                         FUN=function(x) x[2] < x[1] & x[2] < x[3], 
                         fill=NA, align="center")

# Now extract highs and lows as before
highs <- df_sub %>% filter(high)
lows  <- df_sub %>% filter(low)

# For each day, pick highest high (HHW) and lowest low (LLW)
hhw <- highs %>%
  group_by(date = as.Date(datetime)) %>%
  filter(water_level == max(water_level)) %>%
  ungroup()

llw <- lows %>%
  group_by(date = as.Date(datetime)) %>%
  filter(water_level == min(water_level)) %>%
  ungroup()


# Plotting 
plot(df_sub$datetime, df_sub$water_level, type = "l",
     main = "Cedar Island - Higher High Water (HHW) & Lower Low Water (LLW) - June 2025",
     xlab = "Date", ylab = "Water Level (m)", col = "grey80")
points(hhw$datetime, hhw$water_level, pch = 19, col = "blue", cex = 1)
points(llw$datetime, llw$water_level, pch = 19, col = "red", cex = 1)
legend("topright", legend = c("HHW", "LLW"), col = c("blue", "red"), pch = 19)


# Plot the full series and the HHW/LLW points as before
# Combine HHW and LLW and order by datetime
hl <- rbind(hhw, llw)
hl <- hl[order(hl$datetime), ]

# Plot as before
plot(df_sub$datetime, df_sub$water_level, type = "l", xaxt = "n",
     main = "Cedar Island tide heights June 2025",
     xlab = "Date", ylab = "Water Level (m)", col = "darkgrey")
#points(hhw$datetime, hhw$water_level, pch = 19, col = "blue", cex = 1)
#points(llw$datetime, llw$water_level, pch = 19, col = "red", cex = 1)

# Connect all HHW and LLW points in order
lines(hl$datetime, hl$water_level, col = "purple", lwd = 1)
abline(h = mean(df_sub$water_level, na.rm = TRUE), col = "black", lwd = 1)

# Daily ticks (no labels)
days <- seq(as.Date(min(df_sub$datetime)), as.Date(max(df_sub$datetime)), by = "day")
axis(1, at = as.POSIXct(days, tz = "Canada/Pacific"), labels = FALSE, tck = -0.015)

# Find every 7th day starting at June 2
start_day <- as.Date("2025-06-02")
label_days <- seq(from = start_day, to = max(days), by = 7)

# Add labels for those days only
axis(1, at = as.POSIXct(label_days, tz = "Canada/Pacific"),
     labels = format(label_days, "%b %d"), las = 1, cex.axis = 1, tick = FALSE)

# Add custom x-axis ticks for all days
#days <- seq(as.Date(min(df_sub$datetime)), as.Date(max(df_sub$datetime)), by = "day")
#axis(1, at = as.POSIXct(days, tz = "Canada/Pacific"), cex.axis=0.7, las=2)

legend("bottomleft", legend = c("Hourly water level", "HHW/LLW trajectory", "Grand mean"),
       col = c("darkgrey", "purple", "black"), lty = c(1,1,1))







### Working with currents, separating time by current turns 

lines <- readLines(paste(data_dir, "2025_Blackney_tides_AprtoSep.csv", sep="/"))
lines <- trimws(lines)
lines <- lines[lines != ""]  # Remove blank lines

month_map <- setNames(1:12, c("January","February","March","April","May","June",
                              "July","August","September","October","November","December"))
year <- NA
current_month <- NA
day_vecs <- list()
cur_lines <- c()
i <- 1
while (i <= length(lines)) {
  line <- lines[i]
  if (grepl("^\\d{4}$", line)) {
    year <- as.integer(line)
    i <- i + 1
    next
  }
  if (line %in% names(month_map)) {
    if (!is.null(current_month) && length(cur_lines) > 0) {
      # Split all lines for this month on spaces and flatten into a vector
      tokens <- unlist(strsplit(cur_lines, "\\s+"))
      tokens <- tokens[tokens != ""]
      day_vecs[[current_month]] <- tokens
    }
    current_month <- line
    cur_lines <- c()
    i <- i + 1
    next
  }
  cur_lines <- c(cur_lines, line)
  i <- i + 1
}
# Handle last month
if (!is.null(current_month) && length(cur_lines) > 0) {
  tokens <- unlist(strsplit(cur_lines, "\\s+"))
  tokens <- tokens[tokens != ""]
  day_vecs[[current_month]] <- tokens
}


# Example: tokens for a selected month
day_vecs[["July"]]


#---- Examining tokens 

tokens <- day_vecs[[1]]  # Or use a specific month, e.g., day_vecs[["June"]]

print(tokens[1:20])
grep("^\\d+$", tokens[1:20], value=TRUE)
grepl("^\\d+$", tokens[1:20])

x <- day_vecs[["April"]]
valid_tokens <- grep("^\\d+$", x, value=TRUE)


#---- NOW DO THE VECTORS ... 

tide_df <- data.frame(
  year = integer(), month = integer(), day = integer(), time = character(), stringsAsFactors = FALSE
)

for (mon in names(day_vecs)) {
  tokens <- day_vecs[[mon]]
  # Keep only digit tokens (e.g., days and tide times)
  valid_tokens <- grep("^\\d+$", tokens, value=TRUE)
  current_day <- NA
  tide_times <- c()
  tide_days  <- c()
  i <- 1
  while (i <= length(valid_tokens)) {
    token <- valid_tokens[i]
    if (nchar(token) <= 2) {
      # Day number (assumes days 1-31, possibly with leading zero)
      current_day <- as.integer(token)
      i <- i + 1
      next
    }
    if (nchar(token) == 4 && !is.na(current_day)) {
      # Tide time for current day
      tide_times <- c(tide_times, token)
      tide_days  <- c(tide_days,  current_day)
    }
    i <- i + 1
  }
  if (length(tide_times) > 0) {
    sel_idx <- seq(1, length(tide_times), by = 2)  # every second tide time
    tide_df <- rbind(
      tide_df,
      data.frame(
        year = year,
        month = month_map[mon],
        day = tide_days[sel_idx],
        time = tide_times[sel_idx],
        stringsAsFactors = FALSE
      )
    )
  }
}




# Create POSIXct datetimes
tide_df$datetime <- as.POSIXct(
  sprintf("%04d-%02d-%02d %s", tide_df$year, tide_df$month, tide_df$day, 
          paste0(substr(tide_df$time,1,2), ":", substr(tide_df$time,3,4))),
  tz = "America/Vancouver"
)

head(tide_df)


#---- Some playing around. 
# Will need to examine the days that seem too low or too high if we go ahead with this time step.

# Add a column for the duration of each tide 
tide_df$duration_hours <- c(NA, diff(tide_df$datetime))

x <- tide_df[ tide_df$month == 4,]
plot( x$duration_hours ~ x$datetime, type="l")

# And sum for days as simple check of numbers ... 
# NOTE: This includes the interval to the next day’s first event so a minor boundary issue. 
x <- aggregate(duration_hours ~ year + month + day, data = tide_df, FUN = sum, na.rm = TRUE)
daily_hours <- x[order(x$month, x$day), ]


y <- as.POSIXct( sprintf("%04d-%02d-%02d", x$year, x$month, x$day), tz = "America/Vancouver" )
plot(x$duration_hours ~ y)


head(x)

###################





#---- OLDER Stuff ... 

# Data structure will track a plant over 12 months. 
# All other results will be scaled to/from this. 

# An initial plant - baby sporophyte
iPlant <- dataframe( "frond_l" = 1.0,   # cm
                     "frond_w" = 0.05,  # g
                     "frond_n" = 50     # number of fronds
                     "stype_d" = 0.1)   # stype diameter for mass scaling


# Known parameters
B0 <- 0.01    # Initial biomass (kg)
K <- 10       # Carrying capacity (kg)
B_final <- 10 # Final biomass after t months
t <- 6        # Time (months)

# Logistic growth function
logistic_growth <- function(B0, K, t, r) {
  B_predicted <- K / (1 + ((K - B0) / B0) * exp(-r * t))
  return(B_predicted)
}

# Monthly growth rate ESTIMATE for 6 months - optimal
x <- NULL
for (t in 1:6) {
y <- logistic_growth( 0.01, 10, t, 1.565)
x <- c( x, y)
}
plot(x)




# Parameters
B0 <- 0.01    # Initial biomass (kg)
K <- 10       # Carrying capacity (kg)
t <- 1:6      # Time steps (months)
T <- c(10, 12, 15, 18, 20, 22)  # Example temperatures (°C)

# Environmental influence on growth rate
r_max <- 0.8  # Maximum growth rate
T_opt <- 20   # Optimal temperature (°C)
T_range <- 10 # Temperature range for growth (°C)


# Temperature-dependent growth rate function
r_temp <- function(T) {
  r_max * (1 - ((T - T_opt)^2) / (T_range^2))
}


# Logistic growth model with temperature influence
logistic_growth_temp <- function(t, T) {
  r_T <- r_temp(T)  # Calculate growth rate for each temperature
  B_t <- K / (1 + ((K - B0) / B0) * exp(-r_T * t))
  return(B_t)
}

# Apply model over time and temperatures
biomass <- sapply(T, function(temp) logistic_growth_temp(t, temp))


# Plot results
matplot(t, biomass, type = "l", lty = 1, col = rainbow(length(T)),
        xlab = "Time (months)", ylab = "Biomass (kg)",
        main = "Temperature-Dependent Logistic Growth")
legend("bottomright", legend = paste("T =", T, "°C"), col = rainbow(length(T)), lty = 1)


# And then we  scale back up to aPatch using our own allometric scaling ... whew. 



From Pontier et al. 2024:
Measured growth at 3 sites across 4 years
  
Mean growth rates per sampling events across all sites and years:
Blades : 0.83 ± 0.20 cm/day to 8.25 ± 0.78 cm/day
stipes : 0.23 ± 0.15 cm/day to 9.31 ± 0.75 cm/day


aPatch



