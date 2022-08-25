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
# at 50 m
wind50 <- tbl_tot %>%
  filter(u50.MAST != 999.99)

# North
nwind50 <- wind50 %>%
  filter(d50<45|d50>315)

nn50 <- nrow(nwind50)
nmbe50  <- round(tdStats(nwind50$u50.ERA5, nwind50$u50.MAST, "mbe"), 2)
nrmse50 <- round(tdStats(nwind50$u50.ERA5, nwind50$u50.MAST, "rmse"), 2)
ncorcoef50 <- round(cor(nwind50$u50.ERA5, nwind50$u50.MAST), 2)
nr2_50 <- round(tdStats(nwind50$u50.ERA5, nwind50$u50.MAST, "r2"), 2)
nvt50 <- var.test(nwind50$u50.ERA5, nwind50$u50.MAST)

nln_50m <- data.frame(as.integer(nn50), nmbe50, nrmse50, ncorcoef50, nr2_50, 
                      sprintf("%.4f", nvt50$p.value))

setwd("./figures/")

ng1 <- ggplot(nwind50, aes(x = u50.MAST, y = u50.ERA5)) +
  labs(
    title    = "Wind speed at 50 m", 
    subtitle = "North",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", nln_50m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", nln_50m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", nln_50m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", nln_50m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", nln_50m[1,6])), size = 4, hjust=0)
plot(ng1)
ggsave("wind_50m_north.png", width = 5, height = 5)

#---------------------------------------------------
# East
ewind50 <- wind50 %>%
  filter(d50>45&d50<135)

en50 <- nrow(ewind50)
embe50  <- round(tdStats(ewind50$u50.ERA5, ewind50$u50.MAST, "mbe"), 2)
ermse50 <- round(tdStats(ewind50$u50.ERA5, ewind50$u50.MAST, "rmse"), 2)
ecorcoef50 <- round(cor(ewind50$u50.ERA5, ewind50$u50.MAST), 2)
er2_50 <- round(tdStats(ewind50$u50.ERA5, ewind50$u50.MAST, "r2"), 2)
evt50 <- var.test(ewind50$u50.ERA5, ewind50$u50.MAST)

eln_50m <- data.frame(as.integer(en50), embe50, ermse50, ecorcoef50, er2_50, 
                      sprintf("%.4f", evt50$p.value))

eg1 <- ggplot(ewind50, aes(x = u50.MAST, y = u50.ERA5)) +
  labs(
    title    = "Wind speed at 50 m", 
    subtitle = "East",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", eln_50m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", eln_50m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", eln_50m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", eln_50m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", eln_50m[1,6])), size = 4, hjust=0)
plot(eg1)
ggsave("wind_50m_east.png", width = 5, height = 5)

#---------------------------------------------------
# South
swind50 <- wind50 %>%
  filter(d50>135&d50<225)

sn50 <- nrow(swind50)
smbe50  <- round(tdStats(swind50$u50.ERA5, swind50$u50.MAST, "mbe"), 2)
srmse50 <- round(tdStats(swind50$u50.ERA5, swind50$u50.MAST, "rmse"), 2)
scorcoef50 <- round(cor(swind50$u50.ERA5, swind50$u50.MAST), 2)
sr2_50 <- round(tdStats(swind50$u50.ERA5, swind50$u50.MAST, "r2"), 2)
svt50 <- var.test(swind50$u50.ERA5, swind50$u50.MAST)

sln_50m <- data.frame(as.integer(sn50), smbe50, srmse50, scorcoef50, sr2_50, 
                      sprintf("%.4f", svt50$p.value))

sg1 <- ggplot(swind50, aes(x = u50.MAST, y = u50.ERA5)) +
  labs(
    title    = "Wind speed at 50 m", 
    subtitle = "South",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", sln_50m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", sln_50m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", sln_50m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", sln_50m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", sln_50m[1,6])), size = 4, hjust=0)
plot(sg1)
ggsave("wind_50m_south.png", width = 5, height = 5)

#---------------------------------------------------
# West
wwind50 <- wind50 %>%
  filter(d50>225 & d50<315)

wn50 <- nrow(wwind50)
wmbe50  <- round(tdStats(wwind50$u50.ERA5, wwind50$u50.MAST, "mbe"), 2)
wrmse50 <- round(tdStats(wwind50$u50.ERA5, wwind50$u50.MAST, "rmse"), 2)
wcorcoef50 <- round(cor(wwind50$u50.ERA5, wwind50$u50.MAST), 2)
wr2_50 <- round(tdStats(wwind50$u50.ERA5, wwind50$u50.MAST, "r2"), 2)
wvt50 <- var.test(wwind50$u50.ERA5, wwind50$u50.MAST)

wln_50m <- data.frame(as.integer(wn50), wmbe50, wrmse50, wcorcoef50, wr2_50, 
                      sprintf("%.4f", wvt50$p.value))

wg1 <- ggplot(wwind50, aes(x = u50.MAST, y = u50.ERA5)) +
  labs(
    title    = "Wind speed at 50 m", 
    subtitle = "West",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", wln_50m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", wln_50m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", wln_50m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", wln_50m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", wln_50m[1,6])), size = 4, hjust=0)
plot(wg1)
ggsave("wind_50m_west.png", width = 5, height = 5)

