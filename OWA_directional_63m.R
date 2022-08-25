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
  #       CUP_25m_WSave, VANE_25m_WDave, orbit) %>%
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
# at 63 m
wind63 <- tbl_tot %>%
  filter(u63.MAST != 999.99)

# North
nwind63 <- wind63 %>%
  filter(d61<45|d61>315)

nn63 <- nrow(nwind63)
nmbe63  <- round(tdStats(nwind63$u63.ERA5, nwind63$u63.MAST, "mbe"), 2)
nrmse63 <- round(tdStats(nwind63$u63.ERA5, nwind63$u63.MAST, "rmse"), 2)
ncorcoef63 <- round(cor(nwind63$u63.ERA5, nwind63$u63.MAST), 2)
nr2_63 <- round(tdStats(nwind63$u63.ERA5, nwind63$u63.MAST, "r2"), 2)
nvt63 <- var.test(nwind63$u63.ERA5, nwind63$u63.MAST)

nln_63m <- data.frame(as.integer(nn63), nmbe63, nrmse63, ncorcoef63, nr2_63, 
                      sprintf("%.4f", nvt63$p.value))

setwd("./figures/")

ng1 <- ggplot(nwind63, aes(x = u63.MAST, y = u63.ERA5)) +
  labs(
    title    = "Wind speed at 63 m", 
    subtitle = "North",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", nln_63m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", nln_63m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", nln_63m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", nln_63m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", nln_63m[1,6])), size = 4, hjust=0)
plot(ng1)
ggsave("wind_63m_north.png", width = 5, height = 5)

#---------------------------------------------------
# East
ewind63 <- wind63 %>%
  filter(d61>45&d61<135)

en63 <- nrow(ewind63)
embe63  <- round(tdStats(ewind63$u63.ERA5, ewind63$u63.MAST, "mbe"), 2)
ermse63 <- round(tdStats(ewind63$u63.ERA5, ewind63$u63.MAST, "rmse"), 2)
ecorcoef63 <- round(cor(ewind63$u63.ERA5, ewind63$u63.MAST), 2)
er2_63 <- round(tdStats(ewind63$u63.ERA5, ewind63$u63.MAST, "r2"), 2)
evt63 <- var.test(ewind63$u63.ERA5, ewind63$u63.MAST)

eln_63m <- data.frame(as.integer(en63), embe63, ermse63, ecorcoef63, er2_63, 
                      sprintf("%.4f", evt63$p.value))

eg1 <- ggplot(ewind63, aes(x = u63.MAST, y = u63.ERA5)) +
  labs(
    title    = "Wind speed at 63 m", 
    subtitle = "East",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", eln_63m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", eln_63m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", eln_63m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", eln_63m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", eln_63m[1,6])), size = 4, hjust=0)
plot(eg1)
ggsave("wind_63m_east.png", width = 5, height = 5)

#---------------------------------------------------
# South
swind63 <- wind63 %>%
  filter(d61>135&d61<225)

sn63 <- nrow(swind63)
smbe63  <- round(tdStats(swind63$u63.ERA5, swind63$u63.MAST, "mbe"), 2)
srmse63 <- round(tdStats(swind63$u63.ERA5, swind63$u63.MAST, "rmse"), 2)
scorcoef63 <- round(cor(swind63$u63.ERA5, swind63$u63.MAST), 2)
sr2_63 <- round(tdStats(swind63$u63.ERA5, swind63$u63.MAST, "r2"), 2)
svt63 <- var.test(swind63$u63.ERA5, swind63$u63.MAST)

sln_63m <- data.frame(as.integer(sn63), smbe63, srmse63, scorcoef63, sr2_63, 
                      sprintf("%.4f", svt63$p.value))

sg1 <- ggplot(swind63, aes(x = u63.MAST, y = u63.ERA5)) +
  labs(
    title    = "Wind speed at 63 m", 
    subtitle = "South",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", sln_63m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", sln_63m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", sln_63m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", sln_63m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", sln_63m[1,6])), size = 4, hjust=0)
plot(sg1)
ggsave("wind_63m_south.png", width = 5, height = 5)

#---------------------------------------------------
# West
wwind63 <- wind63 %>%
  filter(d61>225 & d61<315)

wn63 <- nrow(wwind63)
wmbe63  <- round(tdStats(wwind63$u63.ERA5, wwind63$u63.MAST, "mbe"), 2)
wrmse63 <- round(tdStats(wwind63$u63.ERA5, wwind63$u63.MAST, "rmse"), 2)
wcorcoef63 <- round(cor(wwind63$u63.ERA5, wwind63$u63.MAST), 2)
wr2_63 <- round(tdStats(wwind63$u63.ERA5, wwind63$u63.MAST, "r2"), 2)
wvt63 <- var.test(wwind63$u63.ERA5, wwind63$u63.MAST)

wln_63m <- data.frame(as.integer(wn63), wmbe63, wrmse63, wcorcoef63, wr2_63, 
                      sprintf("%.4f", wvt63$p.value))

wg1 <- ggplot(wwind63, aes(x = u63.MAST, y = u63.ERA5)) +
  labs(
    title    = "Wind speed at 63 m", 
    subtitle = "West",
    x        = "MAST",
    y        = "SAR (ERA5)",
    caption  = "m/s"
  ) + 
  xlim(0, 20) + ylim(0,20) +
  geom_abline(linetype="dotted") +
  geom_point() +
  annotate(geom="text", x=0.2, y=20, 
           label=as.character(sprintf("N: %3d", wln_63m[1,1])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=19, 
           label=as.character(sprintf("RMSE: %.2f", wln_63m[1,3])), size = 4, hjust=0)  +
  annotate(geom="text", x=0.2, y=18, 
           label=as.character(sprintf("Bias: %.2f", wln_63m[1,2])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=17, 
           label=as.character(sprintf("Cor.: %.2f", wln_63m[1,4])), size = 4, hjust=0) +
  annotate(geom="text", x=0.2, y=16, 
           label=as.character(sprintf("p-value: %s", wln_63m[1,6])), size = 4, hjust=0)
plot(wg1)
ggsave("wind_63m_west.png", width = 5, height = 5)

