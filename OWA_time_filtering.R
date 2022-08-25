# Clean up
rm(list = ls(all.names = TRUE))

setwd("~/RWorkspace/")

#load libraries 
library(tidyverse)
library(lubridate)
library(tdr)
library(ggplot2)

# Mast location
#https://onlinelibrary.wiley.com/doi/pdf/10.1002/we.2757
mlat <- 40. + 55./60 + 28.88/3600.
mlon <- 141. + 24./60. + 25.48/3600.

#---------------------------------------------------
# Read MAST data
tbl_mast <- as_tibble(read_csv("~/OWA/Dataset_All_ver1.0.csv",
                               col_names = TRUE, show_col_types = FALSE)) %>%
  rename(orbit="Asc/Dsc") %>%
  rename(u15=SAT_15m_WSave, u25=CUP_25m_WSave, d15=SAT_15m_WDave, 
         u50=CUP_50m_WSave, d50=VANE_50m_WDave, d61=VANE_61m_WDave,
         u63=CUP_63m_WSave)

# Round time to match for satellite observations  
T_org_10 <- round_date(tbl_mast[['Time_org']], unit = "10 minute")
T_org_00 <- round_date(tbl_mast[['Time_org']], unit = "hour")

# Filtering: obs at the top of the hour  
tbl_mast10 <- tbl_mast %>%
  mutate(T_org_10 = T_org_10, time = T_org_00) %>%
  filter(T_org_10 == Time_UTC) %>%
  select(-T_org_10)

# Norrowing down
tbl_mast_hr <- tbl_mast %>%
  mutate(T_org_10 = T_org_10, time = T_org_00) %>%
  filter(abs(difftime(T_org_10, tbl_mast[['Time_UTC']], units = "mins")) < 60) %>%
  select(-T_org_10)

#---------------------------------------------------
# Read SAR Wind Speed based on ERA5 and extrapolate

# Simple log extrapolation
# Uz = Uh * (np.log((z - d) / z0) / np.log((h - d) / z0))

h   <- 10.   # base height (m)
z15 <- 15.   # target height  = 15 (m)
z25 <- 25.   # target height  = 25 (m)
z50 <- 50.   # target height  = 50 (m)
z63 <- 63.   # target height  = 63 (m)
z0 <- 0.0002 # surface roughness (m)
d <- 0.      # zero plane displacement (m)


tbl_era5 <- as_tibble(read_csv("~/OWA/era5_40.92468889_141.40707778_new.csv",
                               col_names = TRUE, show_col_types = FALSE)) %>%
  mutate(kind = "ERA5", dist=sqrt((lat-mlat)^2 + (lon-mlon)^2)) %>%
  select(lat, lon, time, u10=speed10m_inverted, dist, kind) %>%
  filter(!is.na(u10), .preserve = TRUE) %>%
  mutate(u15=u10*(log((z15-d)/z0)/log((h-d)/z0)),
         u25=u10*(log((z25-d)/z0)/log((h-d)/z0)),
         u50=u10*(log((z50-d)/z0)/log((h-d)/z0)),
         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))

#---------------------------------------------------
# Find the nearest point to the MAST and extract the data
tvalid <- tbl_era5 %>%
  group_by(time) %>%
  summarise(mindist = min(dist))

target_era5 <- tbl_era5 %>%
  filter(dist %in% tvalid$mindist, .preserve = TRUE) %>%
  distinct(time, .keep_all = TRUE)

#---------------------------------------------------
# Combine the mast and SAR wind data
tbl_tot <- inner_join(tbl_mast10, target_era5, by=c("time")) %>%
  rename(u10.ERA5=u10, u15.MAST=u15.x, u15.ERA5=u15.y,
         u25.MAST=u25.x, u25.ERA5=u25.y,
         u50.MAST=u50.x, u50.ERA5=u50.y,
         u63.MAST=u63.x, u63.ERA5=u63.y) %>%
  select(-dist, -kind) 

#---------------------------------------------------
# Basic statistics retrieval

# at 15 m 
wind15 <- tbl_tot %>%
  filter(u15.MAST != 999.99)    # Exclude invalid data

n15 <- nrow(wind15)
mbe0   <- round(tdStats(wind15$u10.ERA5, wind15$u15.MAST, "mbe"), 2)
mbe15  <- round(tdStats(wind15$u15.ERA5, wind15$u15.MAST, "mbe"), 2)
rmse0  <- round(tdStats(wind15$u10.ERA5, wind15$u15.MAST, "rmse"), 2)
rmse15 <- round(tdStats(wind15$u15.ERA5, wind15$u15.MAST, "rmse"), 2)
corcoef0  <- round(cor(wind15$u10.ERA5, wind15$u15.MAST), 2)
corcoef15 <- round(cor(wind15$u15.ERA5, wind15$u15.MAST), 2)
r2_0  <- round(tdStats(wind15$u10.ERA5, wind15$u15.MAST, "r2"), 2)
r2_15 <- round(tdStats(wind15$u15.ERA5, wind15$u15.MAST, "r2"), 2)
vt0 <- var.test(wind15$u10.ERA5, wind15$u15.MAST)
vt15 <- var.test(wind15$u15.ERA5, wind15$u15.MAST)

ln_lwst <- data.frame(as.integer(n15), mbe0, rmse0, corcoef0, r2_0, 
                      sprintf("%.4f", vt0$p.value))
ln_15m <- data.frame(as.integer(n15), mbe15, rmse15, corcoef15, r2_15, 
                     sprintf("%.4f", vt15$p.value))
colnames(ln_lwst) <- c("num","MBE","RMSE","Cor.coeff","R^2","p-value")
colnames(ln_15m) <- c("num","MBE","RMSE","Cor.coeff","R^2","p-value")

timestamps15 <- wind15 %>%
  select(time_sat_obs = Time_org, time_reanalysis=time, time_mast_obs=Time_UTC,
         latitude=lat, longitude=lon)
write.csv(timestamps15, "./csv_files/timestamps_15m.csv", quote=F, row.names=F)

#---
# at 25 m 
wind25 <- tbl_tot %>%
  filter(u25.MAST != 999.99)

n25 <- nrow(wind25)
mbe25  <- round(tdStats(wind25$u25.ERA5, wind25$u25.MAST, "mbe"), 2)
rmse25 <- round(tdStats(wind25$u25.ERA5, wind25$u25.MAST, "rmse"), 2)
corcoef25 <- round(cor(wind25$u25.MAST, wind25$u25.ERA5), 2)
r2_25 <- round(tdStats(wind25$u25.ERA5, wind25$u25.MAST, "r2"), 2)
vt25 <- var.test(wind25$u25.ERA5, wind25$u25.MAST)
ln_25m <- data.frame(as.integer(n25), mbe25, rmse25, corcoef25, r2_25, 
                     sprintf("%.4f",vt25$p.value))
colnames(ln_25m) <- c("num","MBE","RMSE","Cor.coeff","R^2","p-value")

timestamps25 <- wind25 %>%
  select(time_sat_obs = Time_org, time_reanalysis=time, time_mast_obs=Time_UTC,
         latitude=lat, longitude=lon)
write.csv(timestamps25, "./csv_files/timestamps_25m.csv", quote=F, row.names=F)

#---
# at 50 m 
wind50 <- tbl_tot %>%
  filter(u50.MAST != 999.99)

n50 <- nrow(wind50)
mbe50  <- round(tdStats(wind50$u50.ERA5, wind50$u50.MAST, "mbe"), 2)
rmse50 <- round(tdStats(wind50$u50.ERA5, wind50$u50.MAST, "rmse"), 2)
corcoef50 <- round(cor(wind50$u50.MAST, wind50$u50.ERA5), 2)
r2_50 <- round(tdStats(wind50$u50.ERA5, wind50$u50.MAST, "r2"), 2)
vt50 <- var.test(wind50$u50.ERA5, wind50$u50.MAST)

ln_50m <- data.frame(as.integer(n50), mbe50, rmse50, corcoef50, r2_50, 
                     sprintf("%.4f",vt50$p.value))
colnames(ln_50m) <- c("num","MBE","RMSE","Cor.coeff","R^2","p-value")

timestamps50 <- wind50 %>%
  select(time_sat_obs = Time_org, time_reanalysis=time, time_mast_obs=Time_UTC,
         latitude=lat, longitude=lon)
write.csv(timestamps50, "./csv_files/timestamps_50m.csv", quote=F, row.names=F)

#---
# at 63 m 
wind63 <- tbl_tot %>%
  filter(u63.MAST != 999.99)

n63 <- nrow(wind63)
mbe63  <- round(tdStats(wind63$u63.ERA5, wind63$u63.MAST, "mbe"), 2)
rmse63 <- round(tdStats(wind63$u63.ERA5, wind63$u63.MAST, "rmse"), 2)
corcoef63 <- round(cor(wind63$u63.MAST, wind63$u63.ERA5), 2)
r2_63 <- round(tdStats(wind63$u63.ERA5, wind63$u63.MAST, "r2"), 2)
vt63 <- var.test(wind63$u63.ERA5, wind63$u63.MAST)
ln_63m <- data.frame(as.integer(n63), mbe63, rmse63, corcoef63, r2_63, 
                     sprintf("%.4f",vt63$p.value))
colnames(ln_63m) <- c("num","MBE","RMSE","Cor.coeff","R^2","p-value")

timestamps63 <- wind63 %>%
  select(time_sat_obs = Time_org, time_reanalysis=time, time_mast_obs=Time_UTC,
         latitude=lat, longitude=lon)
write.csv(timestamps63, "./csv_files/timestamps_61m.csv", quote=F, row.names=F)

#---------------------------------------------------
# Basic Statistics
stats <- bind_rows(ln_lwst, ln_15m, ln_25m, ln_50m, ln_63m) 
rownames(stats) <- c("Lowest", "15m", "25m", "50m", "63m")
write.csv(stats, "./csv_files/basic_stats.csv", quote=F)

#---------------------------------------------------
# Plotting

# at the lowest level
line1 <- as.character(sprintf("N: %3d", n15))
line2 <- as.character(sprintf("RMSE: %.2f", rmse0))
line3 <- as.character(sprintf("Bias: %.2f", mbe0))
line4 <- as.character(sprintf("Cor: %.2f", corcoef0))
line5 <- as.character(sprintf("p-value: %.4f", vt0$p.value))
g0 <- ggplot(wind15, aes(x = u15.MAST, y = u10.ERA5)) +
  labs(
    title    = "Wind speed at the lowest level", 
    x        = "MAST @ 15 m",
    y        = "SAR (ERA5) @ 10 m",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=1, y=20, label=line1, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=19, label=line2, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=18, label=line3, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=17, label=line4, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=16, label=line5, size = 4, hjust=0) 
plot(g0)
ggsave("./figures/wind_lowest.png", width = 5, height = 5)

#-----
# at 15 m
line1 <- as.character(sprintf("N: %3d", n15))
line2 <- as.character(sprintf("RMSE: %.2f", rmse15))
line3 <- as.character(sprintf("Bias: %.2f", mbe15))
line4 <- as.character(sprintf("Cor: %.2f", corcoef15))
line5 <- as.character(sprintf("p-value: %.4f", vt15$p.value))

g1 <- ggplot(wind15, aes(x = u15.MAST, y = u15.ERA5)) +
  labs(
    title    = "Wind speed at 15 m", 
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=1, y=20, label=line1, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=19, label=line2, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=18, label=line3, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=17, label=line4, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=16, label=line5, size = 4, hjust=0)
plot(g1)
ggsave("./figures/wind_15m.png", width = 5, height = 5)
#-----
# at 25 m
line1 <- as.character(sprintf("N: %3d", n25))
line2 <- as.character(sprintf("RMSE: %.2f", rmse25))
line3 <- as.character(sprintf("Bias: %.2f", mbe25))
line4 <- as.character(sprintf("Cor: %.2f", corcoef25))
line5 <- as.character(sprintf("p-value: %.4f", vt25$p.value))

g2 <- ggplot(wind25, aes(x = u25.MAST, y = u25.ERA5)) +
  labs(
    title    = "Wind speed at 25 m", 
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=1, y=20, label=line1, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=19, label=line2, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=18, label=line3, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=17, label=line4, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=16, label=line5, size = 4, hjust=0)
plot(g2)
ggsave("./figures/wind_25m.png", width = 5, height = 5)
#-----
# at 50 m
line1 <- as.character(sprintf("N: %3d", n50))
line2 <- as.character(sprintf("RMSE: %.2f", rmse50))
line3 <- as.character(sprintf("Bias: %.2f", mbe50))
line4 <- as.character(sprintf("Cor: %.2f", corcoef50))
line5 <- as.character(sprintf("p-value: %.4f", vt50$p.value))

g3 <- ggplot(wind50, aes(x = u50.MAST, y = u50.ERA5)) +
  labs(
    title    = "Wind speed at 50 m", 
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=1, y=20, label=line1, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=19, label=line2, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=18, label=line3, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=17, label=line4, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=16, label=line5, size = 4, hjust=0)
plot(g3)
ggsave("./figures/wind_50m.png", width = 5, height = 5)
#-----
# at 63 m
line1 <- as.character(sprintf("N: %3d", n63))
line2 <- as.character(sprintf("RMSE: %.2f", rmse63))
line3 <- as.character(sprintf("Bias: %.2f", mbe63))
line4 <- as.character(sprintf("Cor: %.2f", corcoef63))
line5 <- as.character(sprintf("p-value: %.4f", vt63$p.value))

g4 <- ggplot(wind63, aes(x = u63.MAST, y = u63.ERA5)) +
  labs(
    title    = "Wind speed at 63 m", 
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=1, y=20, label=line1, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=19, label=line2, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=18, label=line3, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=17, label=line4, size = 4, hjust=0) +
  annotate(geom="text", x=1, y=16, label=line5, size = 4, hjust=0)
plot(g4)
ggsave("./figures/wind_63m.png", width = 5, height = 5)

# end