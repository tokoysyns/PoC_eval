rm(list = ls(all.names = TRUE))
setwd("~/RWorkspace/")

#load libraries 
library(tidyverse)
library(lubridate)
#library(geosphere)
library(tdr)
library(lattice)
library(ggplot2)
library(gridExtra)

# Mast location
#https://onlinelibrary.wiley.com/doi/pdf/10.1002/we.2757
mlat <- 40. + 55./60 + 28.88/3600.
mlon <- 141. + 24./60. + 25.48/3600.

#---------------------------------------------------
# Read MAST data
tbl_mast <- as_tibble(read_csv("~/OWA/Dataset_All_ver1.0.csv",
                               col_names = TRUE, show_col_types = FALSE)) %>%
  rename(orbit="Asc/Dsc") %>%
  #filter(flag_SAT==0, .preserve = TRUE) %>%
  #select(Time_org, Time_UTC, SAT_15m_WDave, SAT_15m_WSave, 
  #       CUP_50m_WSave, VANE_50m_WDave, orbit) %>%
  rename(u15=SAT_15m_WSave, u25=CUP_25m_WSave, d15=SAT_15m_WDave, 
         u50=CUP_50m_WSave, d50=VANE_50m_WDave, d61=VANE_61m_WDave,
         u63=CUP_63m_WSave)

T_org_10 <- round_date(tbl_mast[['Time_org']], unit = "10 minute")
T_org_00 <- round_date(tbl_mast[['Time_org']], unit = "hour")

tbl_mast10 <- tbl_mast %>%
  mutate(T_org_10 = T_org_10, time = T_org_00) %>%
  filter(T_org_10 == Time_UTC) %>%
  select(-T_org_10)
#filter(abs(difftime(T_org_10, tbl_mast[['Time_UTC']], units = "mins")) < 10) %>%
#select(time, everything())

tbl_mast_hr <- tbl_mast %>%
  mutate(T_org_10 = T_org_10, time = T_org_00) %>%
  filter(abs(difftime(T_org_10, tbl_mast[['Time_UTC']], units = "mins")) < 60) %>%
  select(-T_org_10)

#---------------------------------------------------
# Read SAR Wind Speed based on ERA5 and extrapolate

#def log_wind(Uh, h=10, z=15, z0=0.0002, d=0):
#  '''
#    Uh = wind speed at h m
#    Uz = target wind speed at z m
#    h = base height (m)
#    z = target height (m)
#    d = zero plane displacement (m)
#    z0 = surface roughness (m)
#    '''
#Uz = Uh * (np.log((z - d) / z0) / np.log((h - d) / z0))

#def log_wind(Uh, h=10, z=15, z0=0.0002, d=0):
#  '''
#    Uh = wind speed at h m
#    Uz = target wind speed at z m
#    h = base height (m)
#    z = target height (m)
#    d = zero plane displacement (m)
#    z0 = surface roughness (m)
#    '''
#Uz = Uh * (np.log((z - d) / z0) / np.log((h - d) / z0))

h <- 10.
z15 <- 15.
z25 <- 25.
z50 <- 50.
z63 <- 63.
z0 <- 0.0002
d <- 0.

tbl_era5 <- as_tibble(read_csv("~/OWA/era5_40.92468889_141.40707778_new.csv",
                               col_names = TRUE, show_col_types = FALSE)) %>%
  mutate(kind = "ERA5", dist=sqrt((lat-mlat)^2 + (lon-mlon)^2)) %>%
  select(lat, lon, time, u10=speed10m_inverted, dist, kind) %>%
  filter(!is.na(u10), .preserve = TRUE) %>%
  mutate(u15=u10*(log((z15-d)/z0)/log((h-d)/z0)),
         u25=u10*(log((z25-d)/z0)/log((h-d)/z0)),
         u50=u10*(log((z50-d)/z0)/log((h-d)/z0)),
         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    # log_wind

#tbl_era5 <- as_tibble(read_csv("~/OWA/era5_40.92468889_141.40707778.csv",
#                               col_names = TRUE, show_col_types = FALSE)) %>%
#  mutate(kind = "ERA5", dist=sqrt((lat-mlat)^2 + (lon-mlon)^2)) %>%
#  select(lat, lon, time, u10=speed10m_inverted, 
#         u15=extrapolated_15, u20=extrapolated_20,
#         u50=extrapolated_50, u58=extrapolated_58, 
#         dist, kind) %>%
#  mutate(u25=u10*(log((z25-d)/z0)/log((h-d)/z0)),
#         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    # log_wind

#---------------------------------------------------
# Find the nearest point to the MAST
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
  select(-dist, -kind) %>%
  mutate(u15.MAST = if_else(u15.MAST == 999.99, NA_real_, u15.MAST)) %>%
  mutate(u25.MAST = if_else(u25.MAST == 999.99, NA_real_, u25.MAST)) %>%
  mutate(u50.MAST = if_else(u50.MAST == 999.99, NA_real_, u50.MAST)) %>%
  mutate(u63.MAST = if_else(u63.MAST == 999.99, NA_real_, u63.MAST)) 
#---------------------------------------------------
# Time series
setwd("./figures/")

temp <- tbl_tot %>%
  filter(!is.na(u15.MAST), .preserve = TRUE)

tmin <- min(temp$time)
tmax <- max(temp$time)

g1 <- ggplot(temp, aes(x = time), colour="Height") +
  labs(
    #title    = "Wind speed", 
    x        = "Time",
    y        = "m/s",
    color    = "Wind Speed"
  ) + 
  xlim(tmin, tmax) + ylim(0, 20) +
  geom_line(aes(y = u10.ERA5, colour = "SAR @ 10 m")) +
  geom_line(aes(y = u15.MAST, colour = "MAST @ 15 m")) +
  theme(legend.position = "top", legend.justification = c("left"))
plot(g1)
ggsave(file = "wind_time_series_lowest.png", plot = g1, width = 9, height = 4)


g2 <- ggplot(temp, aes(x = time), colour="Height") +
  labs(
    #title    = "Wind speed", 
    x        = "Time",
    y        = "m/s",
    color    = "Wind Speed"
  ) + 
  xlim(tmin, tmax) + ylim(0, 20) +
  geom_line(aes(y = u15.ERA5, colour = "SAR @ 15 m")) +
  geom_line(aes(y = u15.MAST, colour = "MAST @ 15 m")) +
  theme(legend.position = "top", legend.justification = c("left"))
plot(g2)
ggsave(file = "wind_time_series_15m.png", plot = g2, width = 9, height = 4)


temp <- tbl_tot %>%
  filter(!is.na(u25.MAST), .preserve = TRUE)

g3 <- ggplot(temp, aes(x = time), colour="Height") +
  labs(
    #title    = "Wind speed", 
    x        = "Time",
    y        = "m/s",
    color    = "Wind Speed"
  ) + 
  xlim(tmin, tmax) + ylim(0, 20) +
  geom_line(aes(y = u25.ERA5, colour = "SAR @ 25 m")) +
  geom_line(aes(y = u25.MAST, colour = "MAST @ 25 m")) +
  theme(legend.position = "top", legend.justification = c("left"))
plot(g3)
ggsave(file = "wind_time_series_25m.png", plot = g3, width = 9, height = 4)


temp <- tbl_tot %>%
  filter(!is.na(u50.MAST), .preserve = TRUE)

g4 <- ggplot(temp, aes(x = time), colour="Height") +
  labs(
    #    title    = "Wind speed", 
    x        = "Time",
    y        = "m/s",
    color    = "Wind Speed"
  ) + 
  xlim(tmin, tmax) + ylim(0, 20) +
  geom_line(aes(y = u50.ERA5, colour = "SAR @ 50 m")) +
  geom_line(aes(y = u50.MAST, colour = "MAST @ 50 m")) +
  theme(legend.position = "top", legend.justification = c("left"))
plot(g4)
ggsave(file = "wind_time_series_50m.png", plot = g4, width = 9, height = 4)


temp <- tbl_tot %>%
  filter(!is.na(u63.MAST), .preserve = TRUE)

g5 <- ggplot(temp, aes(x = time), colour="Height") +
  labs(
    #    title    = "Wind speed", 
    x        = "Time",
    y        = "m/s",
    color    = "Wind Speed"
  ) + 
  xlim(tmin, tmax) + ylim(0, 20) +
  geom_line(aes(y = u63.ERA5, colour = "SAR @ 63 m")) +
  geom_line(aes(y = u63.MAST, colour = "MAST @ 63 m")) +
  theme(legend.position = "top", legend.justification = c("left"))
plot(g5)
ggsave(file = "wind_time_series_63m.png", plot = g5, width = 9, height = 4)

