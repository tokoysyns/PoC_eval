rm(list = ls(all.names = TRUE))
setwd("~/RWorkspace/")

#load libraries 
library(tidyverse)
library(lubridate)
library(tdr)
library(ggplot2)
library(knitr)

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

T_org_10 <- round_date(tbl_mast[['Time_org']], unit = "10 minute")
T_org_00 <- round_date(tbl_mast[['Time_org']], unit = "hour")

tbl_mast10 <- tbl_mast %>%
  mutate(T_org_10 = T_org_10, time = T_org_00) %>%
  filter(T_org_10 == Time_UTC) %>%
  select(-T_org_10)

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

wind15 <- tbl_tot %>%
  filter(u15.MAST!=999.99, .preserve = TRUE) %>%
  mutate(rng10.ERA5 = case_when(
    (u10.ERA5<1) ~ "~1",
    (u10.ERA5>=1 & u10.ERA5<2) ~ "~2",
    (u10.ERA5>=2 & u10.ERA5<3) ~ "~3",
    (u10.ERA5>=3 & u10.ERA5<4) ~ "~4",
    (u10.ERA5>=4 & u10.ERA5<5) ~ "~5",
    (u10.ERA5>=5 & u10.ERA5<6) ~ "~6",
    (u10.ERA5>=6 & u10.ERA5<7) ~ "~7",
    (u10.ERA5>=7 & u10.ERA5<8) ~ "~8",
    (u10.ERA5>=8 & u10.ERA5<9) ~ "~9",
    (u10.ERA5>=9 & u10.ERA5<10) ~ "~10",
    (u10.ERA5>=10 & u10.ERA5<11) ~ "~11",
    (u10.ERA5>=11 & u10.ERA5<12) ~ "~12",
    (u10.ERA5>=12 & u10.ERA5<13) ~ "~13",
    (u10.ERA5>=13 & u10.ERA5<14) ~ "~14",
    (u10.ERA5>=14 & u10.ERA5<15) ~ "~15",
    (u10.ERA5>=15 & u10.ERA5<16) ~ "~16",
    (u10.ERA5>=16 & u10.ERA5<17) ~ "~17",
    (u10.ERA5>=17 & u10.ERA5<18) ~ "~18",
    (u10.ERA5>=18 & u10.ERA5<19) ~ "~19",
    (u10.ERA5>=19 & u10.ERA5<20) ~ "~20",
    TRUE ~ "20~")) %>%
  mutate(rng15.MAST = case_when(
    (u15.MAST<1) ~ "~1",
    (u15.MAST>=1 & u15.MAST<2) ~ "~2",
    (u15.MAST>=2 & u15.MAST<3) ~ "~3",
    (u15.MAST>=3 & u15.MAST<4) ~ "~4",
    (u15.MAST>=4 & u15.MAST<5) ~ "~5",
    (u15.MAST>=5 & u15.MAST<6) ~ "~6",
    (u15.MAST>=6 & u15.MAST<7) ~ "~7",
    (u15.MAST>=7 & u15.MAST<8) ~ "~8",
    (u15.MAST>=8 & u15.MAST<9) ~ "~9",
    (u15.MAST>=9 & u15.MAST<10) ~ "~10",
    (u15.MAST>=10 & u15.MAST<11) ~ "~11",
    (u15.MAST>=11 & u15.MAST<12) ~ "~12",
    (u15.MAST>=12 & u15.MAST<13) ~ "~13",
    (u15.MAST>=13 & u15.MAST<14) ~ "~14",
    (u15.MAST>=14 & u15.MAST<15) ~ "~15",
    (u15.MAST>=15 & u15.MAST<16) ~ "~16",
    (u15.MAST>=16 & u15.MAST<17) ~ "~17",
    (u15.MAST>=17 & u15.MAST<18) ~ "~18",
    (u15.MAST>=18 & u15.MAST<19) ~ "~19",
    (u15.MAST>=19 & u15.MAST<20) ~ "~20",
    TRUE ~ "20~")) %>%
  mutate(rng15.ERA5 = case_when(
    (u15.ERA5<1) ~ "~1",
    (u15.ERA5>=1 & u15.ERA5<2) ~ "~2",
    (u15.ERA5>=2 & u15.ERA5<3) ~ "~3",
    (u15.ERA5>=3 & u15.ERA5<4) ~ "~4",
    (u15.ERA5>=4 & u15.ERA5<5) ~ "~5",
    (u15.ERA5>=5 & u15.ERA5<6) ~ "~6",
    (u15.ERA5>=6 & u15.ERA5<7) ~ "~7",
    (u15.ERA5>=7 & u15.ERA5<8) ~ "~8",
    (u15.ERA5>=8 & u15.ERA5<9) ~ "~9",
    (u15.ERA5>=9 & u15.ERA5<10) ~ "~10",
    (u15.ERA5>=10 & u15.ERA5<11) ~ "~11",
    (u15.ERA5>=11 & u15.ERA5<12) ~ "~12",
    (u15.ERA5>=12 & u15.ERA5<13) ~ "~13",
    (u15.ERA5>=13 & u15.ERA5<14) ~ "~14",
    (u15.ERA5>=14 & u15.ERA5<15) ~ "~15",
    (u15.ERA5>=15 & u15.ERA5<16) ~ "~16",
    (u15.ERA5>=16 & u15.ERA5<17) ~ "~17",
    (u15.ERA5>=17 & u15.ERA5<18) ~ "~18",
    (u15.ERA5>=18 & u15.ERA5<19) ~ "~19",
    (u15.ERA5>=19 & u15.ERA5<20) ~ "~20",
    TRUE ~ "20~"))

ord<-c("~1","~2","~3","~4","~5","~6","~7","~8","~9","~10",
       "~11","~12","~13","~14","~15","~16","~17","~18","~19", "~20", "20~")

wind15$rng10.ERA5 <- factor(wind15$rng10.ERA5, levels=ord)
wind15$rng15.MAST <- factor(wind15$rng15.MAST, levels=ord)
wind15$rng15.ERA5 <- factor(wind15$rng15.ERA5, levels=ord)

hst10.ERA5 <- wind15 %>%
  group_by(rng10.ERA5, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng10.ERA5)

hst15.MAST <- wind15 %>%
  group_by(rng15.MAST, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng15.MAST)

hst15.ERA5 <- wind15 %>%
  group_by(rng15.ERA5, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng15.ERA5)

#-------------------------------------------------------------
wind25 <- tbl_tot %>%
  filter(u25.MAST!=999.99, .preserve = TRUE) %>%
  mutate(rng25.ERA5 = case_when(
    (u25.ERA5<1) ~ "~1",
    (u25.ERA5>=1 & u25.ERA5<2) ~ "~2",
    (u25.ERA5>=2 & u25.ERA5<3) ~ "~3",
    (u25.ERA5>=3 & u25.ERA5<4) ~ "~4",
    (u25.ERA5>=4 & u25.ERA5<5) ~ "~5",
    (u25.ERA5>=5 & u25.ERA5<6) ~ "~6",
    (u25.ERA5>=6 & u25.ERA5<7) ~ "~7",
    (u25.ERA5>=7 & u25.ERA5<8) ~ "~8",
    (u25.ERA5>=8 & u25.ERA5<9) ~ "~9",
    (u25.ERA5>=9 & u25.ERA5<10) ~ "~10",
    (u25.ERA5>=10 & u25.ERA5<11) ~ "~11",
    (u25.ERA5>=11 & u25.ERA5<12) ~ "~12",
    (u25.ERA5>=12 & u25.ERA5<13) ~ "~13",
    (u25.ERA5>=13 & u25.ERA5<14) ~ "~14",
    (u25.ERA5>=14 & u25.ERA5<15) ~ "~15",
    (u25.ERA5>=15 & u25.ERA5<16) ~ "~16",
    (u25.ERA5>=16 & u25.ERA5<17) ~ "~17",
    (u25.ERA5>=17 & u25.ERA5<18) ~ "~18",
    (u25.ERA5>=18 & u25.ERA5<19) ~ "~19",
    (u25.ERA5>=19 & u25.ERA5<20) ~ "~20",
    TRUE ~ "20~")) %>%
  mutate(rng25.MAST = case_when(
    (u25.MAST<1) ~ "~1",
    (u25.MAST>=1 & u25.MAST<2) ~ "~2",
    (u25.MAST>=2 & u25.MAST<3) ~ "~3",
    (u25.MAST>=3 & u25.MAST<4) ~ "~4",
    (u25.MAST>=4 & u25.MAST<5) ~ "~5",
    (u25.MAST>=5 & u25.MAST<6) ~ "~6",
    (u25.MAST>=6 & u25.MAST<7) ~ "~7",
    (u25.MAST>=7 & u25.MAST<8) ~ "~8",
    (u25.MAST>=8 & u25.MAST<9) ~ "~9",
    (u25.MAST>=9 & u25.MAST<10) ~ "~10",
    (u25.MAST>=10 & u25.MAST<11) ~ "~11",
    (u25.MAST>=11 & u25.MAST<12) ~ "~12",
    (u25.MAST>=12 & u25.MAST<13) ~ "~13",
    (u25.MAST>=13 & u25.MAST<14) ~ "~14",
    (u25.MAST>=14 & u25.MAST<15) ~ "~15",
    (u25.MAST>=15 & u25.MAST<16) ~ "~16",
    (u25.MAST>=16 & u25.MAST<17) ~ "~17",
    (u25.MAST>=17 & u25.MAST<18) ~ "~18",
    (u25.MAST>=18 & u25.MAST<19) ~ "~19",
    (u25.MAST>=19 & u25.MAST<20) ~ "~20",
    TRUE ~ "20~")) 

wind25$rng25.MAST <- factor(wind25$rng25.MAST, levels=ord)
wind25$rng25.ERA5 <- factor(wind25$rng25.ERA5, levels=ord)

hst25.MAST <- wind25 %>%
  group_by(rng25.MAST, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng25.MAST)

hst25.ERA5 <- wind25 %>%
  group_by(rng25.ERA5, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng25.ERA5)

#-------------------------------------------------------------
wind50 <- tbl_tot %>%
  filter(u50.MAST!=999.99, .preserve = TRUE) %>%
  mutate(rng50.ERA5 = case_when(
    (u50.ERA5<1) ~ "~1",
    (u50.ERA5>=1 & u50.ERA5<2) ~ "~2",
    (u50.ERA5>=2 & u50.ERA5<3) ~ "~3",
    (u50.ERA5>=3 & u50.ERA5<4) ~ "~4",
    (u50.ERA5>=4 & u50.ERA5<5) ~ "~5",
    (u50.ERA5>=5 & u50.ERA5<6) ~ "~6",
    (u50.ERA5>=6 & u50.ERA5<7) ~ "~7",
    (u50.ERA5>=7 & u50.ERA5<8) ~ "~8",
    (u50.ERA5>=8 & u50.ERA5<9) ~ "~9",
    (u50.ERA5>=9 & u50.ERA5<10) ~ "~10",
    (u50.ERA5>=10 & u50.ERA5<11) ~ "~11",
    (u50.ERA5>=11 & u50.ERA5<12) ~ "~12",
    (u50.ERA5>=12 & u50.ERA5<13) ~ "~13",
    (u50.ERA5>=13 & u50.ERA5<14) ~ "~14",
    (u50.ERA5>=14 & u50.ERA5<15) ~ "~15",
    (u50.ERA5>=15 & u50.ERA5<16) ~ "~16",
    (u50.ERA5>=16 & u50.ERA5<17) ~ "~17",
    (u50.ERA5>=17 & u50.ERA5<18) ~ "~18",
    (u50.ERA5>=18 & u50.ERA5<19) ~ "~19",
    (u50.ERA5>=19 & u50.ERA5<20) ~ "~20",
    TRUE ~ "20~")) %>%
  mutate(rng50.MAST = case_when(
    (u50.MAST<1) ~ "~1",
    (u50.MAST>=1 & u50.MAST<2) ~ "~2",
    (u50.MAST>=2 & u50.MAST<3) ~ "~3",
    (u50.MAST>=3 & u50.MAST<4) ~ "~4",
    (u50.MAST>=4 & u50.MAST<5) ~ "~5",
    (u50.MAST>=5 & u50.MAST<6) ~ "~6",
    (u50.MAST>=6 & u50.MAST<7) ~ "~7",
    (u50.MAST>=7 & u50.MAST<8) ~ "~8",
    (u50.MAST>=8 & u50.MAST<9) ~ "~9",
    (u50.MAST>=9 & u50.MAST<10) ~ "~10",
    (u50.MAST>=10 & u50.MAST<11) ~ "~11",
    (u50.MAST>=11 & u50.MAST<12) ~ "~12",
    (u50.MAST>=12 & u50.MAST<13) ~ "~13",
    (u50.MAST>=13 & u50.MAST<14) ~ "~14",
    (u50.MAST>=14 & u50.MAST<15) ~ "~15",
    (u50.MAST>=15 & u50.MAST<16) ~ "~16",
    (u50.MAST>=16 & u50.MAST<17) ~ "~17",
    (u50.MAST>=17 & u50.MAST<18) ~ "~18",
    (u50.MAST>=18 & u50.MAST<19) ~ "~19",
    (u50.MAST>=19 & u50.MAST<20) ~ "~20",
    TRUE ~ "20~")) 

wind50$rng50.MAST <- factor(wind50$rng50.MAST, levels=ord)
wind50$rng50.ERA5 <- factor(wind50$rng50.ERA5, levels=ord)

hst50.MAST <- wind50 %>%
  group_by(rng50.MAST, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng50.MAST)

hst50.ERA5 <- wind50 %>%
  group_by(rng50.ERA5, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng50.ERA5)

#-------------------------------------------------------------
wind63 <- tbl_tot %>%
  filter(u63.MAST!=999.99, .preserve = TRUE) %>%
  mutate(rng63.ERA5 = case_when(
    (u63.ERA5<1) ~ "~1",
    (u63.ERA5>=1 & u63.ERA5<2) ~ "~2",
    (u63.ERA5>=2 & u63.ERA5<3) ~ "~3",
    (u63.ERA5>=3 & u63.ERA5<4) ~ "~4",
    (u63.ERA5>=4 & u63.ERA5<5) ~ "~5",
    (u63.ERA5>=5 & u63.ERA5<6) ~ "~6",
    (u63.ERA5>=6 & u63.ERA5<7) ~ "~7",
    (u63.ERA5>=7 & u63.ERA5<8) ~ "~8",
    (u63.ERA5>=8 & u63.ERA5<9) ~ "~9",
    (u63.ERA5>=9 & u63.ERA5<10) ~ "~10",
    (u63.ERA5>=10 & u63.ERA5<11) ~ "~11",
    (u63.ERA5>=11 & u63.ERA5<12) ~ "~12",
    (u63.ERA5>=12 & u63.ERA5<13) ~ "~13",
    (u63.ERA5>=13 & u63.ERA5<14) ~ "~14",
    (u63.ERA5>=14 & u63.ERA5<15) ~ "~15",
    (u63.ERA5>=15 & u63.ERA5<16) ~ "~16",
    (u63.ERA5>=16 & u63.ERA5<17) ~ "~17",
    (u63.ERA5>=17 & u63.ERA5<18) ~ "~18",
    (u63.ERA5>=18 & u63.ERA5<19) ~ "~19",
    (u63.ERA5>=19 & u63.ERA5<20) ~ "~20",
    TRUE ~ "20~")) %>%
  mutate(rng63.MAST = case_when(
    (u63.MAST<1) ~ "~1",
    (u63.MAST>=1 & u63.MAST<2) ~ "~2",
    (u63.MAST>=2 & u63.MAST<3) ~ "~3",
    (u63.MAST>=3 & u63.MAST<4) ~ "~4",
    (u63.MAST>=4 & u63.MAST<5) ~ "~5",
    (u63.MAST>=5 & u63.MAST<6) ~ "~6",
    (u63.MAST>=6 & u63.MAST<7) ~ "~7",
    (u63.MAST>=7 & u63.MAST<8) ~ "~8",
    (u63.MAST>=8 & u63.MAST<9) ~ "~9",
    (u63.MAST>=9 & u63.MAST<10) ~ "~10",
    (u63.MAST>=10 & u63.MAST<11) ~ "~11",
    (u63.MAST>=11 & u63.MAST<12) ~ "~12",
    (u63.MAST>=12 & u63.MAST<13) ~ "~13",
    (u63.MAST>=13 & u63.MAST<14) ~ "~14",
    (u63.MAST>=14 & u63.MAST<15) ~ "~15",
    (u63.MAST>=15 & u63.MAST<16) ~ "~16",
    (u63.MAST>=16 & u63.MAST<17) ~ "~17",
    (u63.MAST>=17 & u63.MAST<18) ~ "~18",
    (u63.MAST>=18 & u63.MAST<19) ~ "~19",
    (u63.MAST>=19 & u63.MAST<20) ~ "~20",
    TRUE ~ "20~")) 

wind63$rng63.MAST <- factor(wind63$rng63.MAST, levels=ord)
wind63$rng63.ERA5 <- factor(wind63$rng63.ERA5, levels=ord)

hst63.MAST <- wind63 %>%
  group_by(rng63.MAST, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng63.MAST)

hst63.ERA5 <- wind63 %>%
  group_by(rng63.ERA5, .drop = FALSE) %>%
  summarise(n=n()) %>%
  rename(class=rng63.ERA5)

#-----------------
# Bind the results

hst <- left_join(hst10.ERA5, hst15.MAST, by=c("class")) %>%
  left_join(hst15.ERA5, by=c("class")) %>%
  left_join(hst25.MAST, by=c("class")) %>%
  left_join(hst25.ERA5, by=c("class")) %>%
  left_join(hst50.MAST, by=c("class")) %>%
  left_join(hst50.ERA5, by=c("class")) %>%
  left_join(hst63.MAST, by=c("class")) %>%
  left_join(hst63.ERA5, by=c("class")) %>%
  rename(SAR.10 = n.x, MAST.15=n.y, SAR.15 = n.x.x, 
         MAST.25=n.y.y, SAR.25 = n.x.x.x,
         MAST.50=n.y.y.y, SAR.50 = n.x.x.x.x,
         MAST.63=n.y.y.y.y, SAR.63 = n)

hst_tot <- hst %>%
  rename(range=class, SAR_10m=SAR.10, 
         MAST_15m=MAST.15, SAR_15m=SAR.15, MAST_25m=MAST.25, 
         SAR_25m=SAR.25, MAST_50m=MAST.50, SAR_50m=SAR.50, 
         MAST_63m=MAST.63, SAR_63m=SAR.63)

hst_lowest <- hst %>%
  pivot_longer(cols = c("SAR.10", "MAST.15", "SAR.15"),
               names_to = "variable", values_to = "value")

hst_25m <- hst %>%
  pivot_longer(cols = c("SAR.25", "MAST.25"),
               names_to = "variable", values_to = "value")

hst_50m <- hst %>%
  pivot_longer(cols = c("SAR.50", "MAST.50"),
               names_to = "variable", values_to = "value")

hst_63m <- hst %>%
  pivot_longer(cols = c("SAR.63", "MAST.63"),
               names_to = "variable", values_to = "value")

write.csv(hst_tot, "hst_by_wind_spd.csv", quote=F, row.names=F)

#----------------------------------------------------
setwd("./figures/")
g1 <- ggplot(hst_lowest, aes(x = value, fill = variable)) + 
  geom_histogram(position = "dodge", binwidth = 1) +
  labs(x="Wind speed [m/s]",
       title="Frequency classified by wind speed",
       subtitle="Height: 10 & 15 m") +
  xlim(0, 20) + ylim(0, 10)
plot(g1)
ggsave("hist_lowest.png", width = 5, height = 5)


g2 <- ggplot(hst_25m, aes(x = value, fill = variable)) + 
  geom_histogram(position = "dodge", binwidth = 1) +
  labs(x="Wind speed [m/s]",
       title="Frequency classified by wind speed",
       subtitle="Height: 25 m") +
  xlim(0, 20) + ylim(0, 10)
plot(g2)
ggsave("hist_25m.png", width = 5, height = 5)

g3 <- ggplot(hst_50m, aes(x = value, fill = variable)) + 
  geom_histogram(position = "dodge", binwidth = 1) +
  labs(x="Wind speed [m/s]",
       title="Frequency classified by wind speed",
       subtitle="Height: 50 m") +
  xlim(0, 20) + ylim(0, 10)
plot(g3)
ggsave("hist_50m.png", width = 5, height = 5)


g4 <- ggplot(hst_63m, aes(x = value, fill = variable)) + 
  geom_histogram(position = "dodge", binwidth = 1) +
  labs(x="Wind speed [m/s]",
       title="Frequency classified by wind speed",
       subtitle="Height: 63 m") +
  xlim(0, 20) + ylim(0, 10)
plot(g4)
ggsave("hist_63m.png", width = 5, height = 5)

