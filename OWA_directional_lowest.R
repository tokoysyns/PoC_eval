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

h <- 10.
z15 <- 15.
z25 <- 25.
z50 <- 50.
z63 <- 63.
z0 <- 0.0002
d <- 0.

#tbl_era5 <- as_tibble(read_csv("~/OWA/era5_40.92468889_141.40707778.csv",
#                               col_names = TRUE, show_col_types = FALSE)) %>%
#  mutate(kind = "ERA5", dist=sqrt((lat-mlat)^2 + (lon-mlon)^2)) %>%
#  select(lat, lon, time, u10=speed10m_inverted, 
#         u15=extrapolated_15, u20=extrapolated_20,
#         u50=extrapolated_50, u58=extrapolated_58, 
#         dist, kind) %>%
#  mutate(u25=u10*(log((z25-d)/z0)/log((h-d)/z0)),
#         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    # log_wind
tbl_era5 <- as_tibble(read_csv("~/OWA/era5_40.92468889_141.40707778_new.csv",
                               col_names = TRUE, show_col_types = FALSE)) %>%
  mutate(kind = "ERA5", dist=sqrt((lat-mlat)^2 + (lon-mlon)^2)) %>%
  select(lat, lon, time, u10=speed10m_inverted, dist, kind) %>%
  filter(!is.na(u10), .preserve = TRUE) %>%
  mutate(u15=u10*(log((z15-d)/z0)/log((h-d)/z0)),
         u25=u10*(log((z25-d)/z0)/log((h-d)/z0)),
         u50=u10*(log((z50-d)/z0)/log((h-d)/z0)),
         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    # log_wind
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
  select(-dist, -kind) 
#---------------------------------------------------
# lowest
wind15 <- tbl_tot %>%
  filter(u15.MAST != 999.99)

# North
nwind15 <- wind15 %>%
  filter(d15<45|d15>315)

nn15 <- nrow(nwind15)
nmbe0   <- round(tdStats(nwind15$u10.ERA5, nwind15$u15.MAST, "mbe"), 2)
nmbe15  <- round(tdStats(nwind15$u15.ERA5, nwind15$u15.MAST, "mbe"), 2)
nrmse0  <- round(tdStats(nwind15$u10.ERA5, nwind15$u15.MAST, "rmse"), 2)
nrmse15 <- round(tdStats(nwind15$u15.ERA5, nwind15$u15.MAST, "rmse"), 2)
ncorcoef0  <- round(cor(nwind15$u10.ERA5, nwind15$u15.MAST), 2)
ncorcoef15 <- round(cor(nwind15$u15.ERA5, nwind15$u15.MAST), 2)
nr2_0  <- round(tdStats(nwind15$u10.ERA5, nwind15$u15.MAST, "r2"), 2)
nr2_15 <- round(tdStats(nwind15$u15.ERA5, nwind15$u15.MAST, "r2"), 2)
nvt0 <- var.test(nwind15$u10.ERA5, nwind15$u15.MAST)
nvt15 <- var.test(nwind15$u15.ERA5, nwind15$u15.MAST)

nln_lwst <- data.frame(as.integer(nn15), nmbe0, nrmse0, ncorcoef0, nr2_0, 
                       sprintf("%.4f", nvt0$p.value))
nln_15m <- data.frame(as.integer(nn15), nmbe15, nrmse15, ncorcoef15, nr2_15, 
                      sprintf("%.4f", nvt15$p.value))

setwd("./figures/")

ng0 <- ggplot(nwind15, aes(x = u15.MAST, y = u10.ERA5)) +
  labs(
    title    = "Wind speed at the lowest level", 
    subtitle = "North",
    x        = "MAST @ 15 m",
    y        = "SAR (ERA5) @ 10 m",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", nln_lwst[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", nln_lwst[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", nln_lwst[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", nln_lwst[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", nln_lwst[1,6])), size = 4, hjust=0)
plot(ng0)
ggsave("wind_lowest_north.png", width = 5, height = 5)

ng1 <- ggplot(nwind15, aes(x = u15.MAST, y = u15.ERA5)) +
  labs(
    title    = "Wind speed at 15 m", 
    subtitle = "North",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", nln_15m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", nln_15m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", nln_15m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", nln_15m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", nln_15m[1,6])), size = 4, hjust=0)
plot(ng1)
ggsave("wind_15m_north.png", width = 5, height = 5)

#---------------------------------------------------
# East
ewind15 <- wind15 %>%
  filter(d15>45&d15<135)

en15 <- nrow(ewind15)
embe0   <- round(tdStats(ewind15$u10.ERA5, ewind15$u15.MAST, "mbe"), 2)
embe15  <- round(tdStats(ewind15$u15.ERA5, ewind15$u15.MAST, "mbe"), 2)
ermse0  <- round(tdStats(ewind15$u10.ERA5, ewind15$u15.MAST, "rmse"), 2)
ermse15 <- round(tdStats(ewind15$u15.ERA5, ewind15$u15.MAST, "rmse"), 2)
ecorcoef0  <- round(cor(ewind15$u10.ERA5, ewind15$u15.MAST), 2)
ecorcoef15 <- round(cor(ewind15$u15.ERA5, ewind15$u15.MAST), 2)
er2_0  <- round(tdStats(ewind15$u10.ERA5, ewind15$u15.MAST, "r2"), 2)
er2_15 <- round(tdStats(ewind15$u15.ERA5, ewind15$u15.MAST, "r2"), 2)
evt0 <- var.test(ewind15$u10.ERA5, ewind15$u15.MAST)
evt15 <- var.test(ewind15$u15.ERA5, ewind15$u15.MAST)

eln_lwst <- data.frame(as.integer(en15), embe0, ermse0, ecorcoef0, er2_0, 
                       sprintf("%.4f", evt0$p.value))
eln_15m <- data.frame(as.integer(en15), embe15, ermse15, ecorcoef15, er2_15, 
                      sprintf("%.4f", evt15$p.value))

eg0 <- ggplot(ewind15, aes(x = u15.MAST, y = u10.ERA5)) +
  labs(
    title    = "Wind speed at the lowest level", 
    subtitle = "East",
    x        = "MAST @ 15 m",
    y        = "SAR (ERA5) @ 10 m",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", eln_lwst[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", eln_lwst[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", eln_lwst[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", eln_lwst[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", eln_lwst[1,6])), size = 4, hjust=0)
plot(eg0)
ggsave("wind_lowest_east.png", width = 5, height = 5)


eg1 <- ggplot(ewind15, aes(x = u15.MAST, y = u15.ERA5)) +
  labs(
    title    = "Wind speed at 15 m", 
    subtitle = "East",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", eln_15m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", eln_15m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", eln_15m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", eln_15m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", eln_15m[1,6])), size = 4, hjust=0)
plot(eg1)
ggsave("wind_15m_east.png", width = 5, height = 5)


#---------------------------------------------------
# South
swind15 <- wind15 %>%
  filter(d15>135&d15<2255)

sn15 <- nrow(swind15)
smbe0   <- round(tdStats(swind15$u10.ERA5, swind15$u15.MAST, "mbe"), 2)
smbe15  <- round(tdStats(swind15$u15.ERA5, swind15$u15.MAST, "mbe"), 2)
srmse0  <- round(tdStats(swind15$u10.ERA5, swind15$u15.MAST, "rmse"), 2)
srmse15 <- round(tdStats(swind15$u15.ERA5, swind15$u15.MAST, "rmse"), 2)
scorcoef0  <- round(cor(swind15$u10.ERA5, swind15$u15.MAST), 2)
scorcoef15 <- round(cor(swind15$u15.ERA5, swind15$u15.MAST), 2)
sr2_0  <- round(tdStats(swind15$u10.ERA5, swind15$u15.MAST, "r2"), 2)
sr2_15 <- round(tdStats(swind15$u15.ERA5, swind15$u15.MAST, "r2"), 2)
svt0 <- var.test(swind15$u10.ERA5, swind15$u15.MAST)
svt15 <- var.test(swind15$u15.ERA5, swind15$u15.MAST)

sln_lwst <- data.frame(as.integer(sn15), smbe0, srmse0, scorcoef0, sr2_0, 
                       sprintf("%.4f", svt0$p.value))
sln_15m <- data.frame(as.integer(sn15), smbe15, srmse15, scorcoef15, sr2_15, 
                      sprintf("%.4f", svt15$p.value))

sg0 <- ggplot(swind15, aes(x = u15.MAST, y = u10.ERA5)) +
  labs(
    title    = "Wind speed at the lowest level", 
    subtitle = "South",
    x        = "MAST @ 15 m",
    y        = "SAR (ERA5) @ 10 m",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", sln_lwst[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", sln_lwst[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", sln_lwst[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", sln_lwst[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", sln_lwst[1,6])), size = 4, hjust=0)
plot(sg0)
ggsave("wind_lowest_south.png", width = 5, height = 5)

sg1 <- ggplot(swind15, aes(x = u15.MAST, y = u15.ERA5)) +
  labs(
    title    = "Wind speed at 15 m", 
    subtitle = "South",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", sln_15m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", sln_15m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", sln_15m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", sln_15m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", sln_15m[1,6])), size = 4, hjust=0)
plot(sg1)
ggsave("wind_15m_south.png", width = 5, height = 5)

#---------------------------------------------------
# West
wwind15 <- wind15 %>%
  filter(d15>225 & d15<315)

wn15 <- nrow(wwind15)
wmbe0   <- round(tdStats(wwind15$u10.ERA5, wwind15$u15.MAST, "mbe"), 2)
wmbe15  <- round(tdStats(wwind15$u15.ERA5, wwind15$u15.MAST, "mbe"), 2)
wrmse0  <- round(tdStats(wwind15$u10.ERA5, wwind15$u15.MAST, "rmse"), 2)
wrmse15 <- round(tdStats(wwind15$u15.ERA5, wwind15$u15.MAST, "rmse"), 2)
wcorcoef0  <- round(cor(wwind15$u10.ERA5, wwind15$u15.MAST), 2)
wcorcoef15 <- round(cor(wwind15$u15.ERA5, wwind15$u15.MAST), 2)
wr2_0  <- round(tdStats(wwind15$u10.ERA5, wwind15$u15.MAST, "r2"), 2)
wr2_15 <- round(tdStats(wwind15$u15.ERA5, wwind15$u15.MAST, "r2"), 2)
wvt0 <- var.test(wwind15$u10.ERA5, wwind15$u15.MAST)
wvt15 <- var.test(wwind15$u15.ERA5, wwind15$u15.MAST)

wln_lwst <- data.frame(as.integer(wn15), wmbe0, wrmse0, wcorcoef0, wr2_0, 
                       sprintf("%.4f", wvt0$p.value))
wln_15m <- data.frame(as.integer(wn15), wmbe15, wrmse15, wcorcoef15, wr2_15, 
                      sprintf("%.4f", wvt15$p.value))

wg0 <- ggplot(wwind15, aes(x = u15.MAST, y = u10.ERA5)) +
  labs(
    title    = "Wind speed at the lowest level", 
    subtitle = "West",
    x        = "MAST @ 15 m",
    y        = "SAR (ERA5) @ 10 m",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", wln_lwst[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", wln_lwst[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", wln_lwst[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", wln_lwst[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", wln_lwst[1,6])), size = 4, hjust=0)
plot(wg0)
ggsave("wind_lowest_west.png", width = 5, height = 5)

wg1 <- ggplot(wwind15, aes(x = u15.MAST, y = u15.ERA5)) +
  labs(
    title    = "Wind speed at 15 m", 
    subtitle = "West",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", wln_15m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", wln_15m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", wln_15m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", wln_15m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", wln_15m[1,6])), size = 4, hjust=0)
plot(wg1)
ggsave("wind_15m_west.png", width = 5, height = 5)

