rm(list = ls(all.names = TRUE))
setwd("~/RWorkspace/")

#load libraries 
library(tidyverse)
library(lubridate)
library(tdr)
library(lattice)
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
  #filter(flag_SAT==0, .preserve = TRUE) %>%
  #select(Time_org, Time_UTC, SAT_15m_WDave, SAT_15m_WSave, 
  #       CUP_50m_WSave, VANE_50m_WDave, orbit) %>%
  rename(u15=SAT_15m_WSave, u25=CUP_25m_WSave, d15=SAT_15m_WDave, 
         u50=CUP_50m_WSave, d50=VANE_50m_WDave, d61=VANE_61m_WDave,
         u63=CUP_63m_WSave) %>%
  mutate(dir15 = case_when(
    d15>360 ~ "invalid",
    (d15<11.25 | d15>348.75) ~ "N",
    (d15>11.25 & d15<33.75) ~ "NNE",
    (d15>33.75 & d15<56.25) ~ "NE",
    (d15>56.25 & d15<78.75) ~ "ENE",
    (d15>78.75 & d15<101.25) ~ "E",
    (d15>101.25 & d15<123.75) ~ "ESE",
    (d15>123.75 & d15<146.25) ~ "SE",
    (d15>146.25 & d15<168.75) ~ "SSE",
    (d15>168.75 & d15<191.25) ~ "S",
    (d15>191.25 & d15<213.75) ~ "SSW",
    (d15>213.75 & d15<236.25) ~ "SW",
    (d15>236.25 & d15<258.75) ~ "WSW",
    (d15>258.75 & d15<281.25) ~ "W",
    (d15>281.25 & d15<303.75) ~ "WNW",
    (d15>303.75 & d15<326.25) ~ "NW",
    TRUE ~ "NNW")) %>%
  mutate(dir50 = case_when(
    d50>360 ~ "invalid",
    (d50<11.25 | d50>348.75) ~ "N",
    (d50>11.25 & d50<33.75) ~ "NNE",
    (d50>33.75 & d50<56.25) ~ "NE",
    (d50>56.25 & d50<78.75) ~ "ENE",
    (d50>78.75 & d50<101.25) ~ "E",
    (d50>101.25 & d50<123.75) ~ "ESE",
    (d50>123.75 & d50<146.25) ~ "SE",
    (d50>146.25 & d50<168.75) ~ "SSE",
    (d50>168.75 & d50<191.25) ~ "S",
    (d50>191.25 & d50<213.75) ~ "SSW",
    (d50>213.75 & d50<236.25) ~ "SW",
    (d50>236.25 & d50<258.75) ~ "WSW",
    (d50>258.75 & d50<281.25) ~ "W",
    (d50>281.25 & d50<303.75) ~ "WNW",
    (d50>303.75 & d50<326.25) ~ "NW",
    TRUE ~ "NNW")) %>%
  mutate(dir61 = case_when(
    d61>360 ~ "invalid",
    (d61<11.25 | d61>348.75) ~ "N",
    (d61>11.25 & d61<33.75) ~ "NNE",
    (d61>33.75 & d61<56.25) ~ "NE",
    (d61>56.25 & d61<78.75) ~ "ENE",
    (d61>78.75 & d61<101.25) ~ "E",
    (d61>101.25 & d61<123.75) ~ "ESE",
    (d61>123.75 & d61<146.25) ~ "SE",
    (d61>146.25 & d61<168.75) ~ "SSE",
    (d61>168.75 & d61<191.25) ~ "S",
    (d61>191.25 & d61<213.75) ~ "SSW",
    (d61>213.75 & d61<236.25) ~ "SW",
    (d61>236.25 & d61<258.75) ~ "WSW",
    (d61>258.75 & d61<281.25) ~ "W",
    (d61>281.25 & d61<303.75) ~ "WNW",
    (d61>303.75 & d61<326.25) ~ "NW",
    TRUE ~ "NNW")) 

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

#-------------------------------------------
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
         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    

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
wind15 <- tbl_tot %>%
  filter(u15.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir15) #%>%
#  summarise(n=n(), u15ave.MAST=mean(u15.MAST), 
#           u10ave.ERA5=mean(u10.ERA5), u15ave.ERA5=mean(u15.ERA5),
#           u15med.MAST=median(u15.MAST), 
#           u10med.ERA5=median(u10.ERA5), u15med.ERA5=median(u15.ERA5))

wind50 <- tbl_tot %>%
  filter(u50.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir50) #%>%
#  summarise(n=n(), u50ave.MAST=mean(u50.MAST),u50ave.ERA5=mean(u50.ERA5),
#            u50med.MAST=median(u50.MAST), u50med.ERA5=median(u50.ERA5))

wind63 <- tbl_tot %>%
  filter(u63.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir61) #%>%
#  summarise(n=n(), u63ave.MAST=mean(u63.MAST),u63ave.ERA5=mean(u63.ERA5),
#            u63med.MAST=median(u63.MAST), u63med.ERA5=median(u63.ERA5))

#-----------------
ord<-c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW")

wind15$dir15<- factor(wind15$dir15, levels=ord)
wind50$dir50<- factor(wind50$dir50, levels=ord)
wind63$dir61<- factor(wind63$dir61, levels=ord)

#-----------------
w15mast <- wind15 %>%
  select(dir15, u15.MAST) %>%
  mutate(cut15 = cut(u15.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))
#w15mast$cut<-cut(w15mast$u15.MAST,breaks=c(0,3,6,10,20,100),labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  ") ,right =TRUE, include.lowest = TRUE) 

# Wind speed distribution
kable(table(w15mast$cut15,w15mast$dir15))

# Frequency distribution
kable(t(apply(table(w15mast$cut15,w15mast$dir15),2,sum)))

# Plotting
setwd("./figures/")

ggplot(w15mast, aes(x =factor(dir15),fill = cut15)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="Mast obs. at 15 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_mast_15m.png", width = 5, height = 5)

#---------------------
w10era5 <- wind15 %>%
  select(dir15, u10.ERA5) %>%
  mutate(cut10 = cut(u10.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w10era5$cut10,w10era5$dir15))

# Frequency distribution
kable(t(apply(table(w10era5$cut10,w10era5$dir15),2,sum)))

# Plotting
ggplot(w10era5, aes(x =factor(dir15),fill = cut10)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 10 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_era5_10m.png", width = 5, height = 5)

#---------------------
w15era5 <- wind15 %>%
  select(dir15, u15.ERA5) %>%
  mutate(cut10 = cut(u15.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w15era5$cut10,w15era5$dir15))

# Frequency distribution
kable(t(apply(table(w15era5$cut10,w15era5$dir15),2,sum)))

# Plotting
ggplot(w15era5, aes(x =factor(dir15),fill = cut10)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 15 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_era5_15m.png", width = 5, height = 5)

#---------------------
w50mast <- wind50 %>%
  select(dir50, u50.MAST) %>%
  mutate(cut50 = cut(u50.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w50mast$cut50,w50mast$dir50))

# Frequency distribution
kable(t(apply(table(w50mast$cut50,w50mast$dir50),2,sum)))

# Plotting
ggplot(w50mast, aes(x =factor(dir50),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="Mast obs. at 50 m",
       subtitle = "Based on the mast wind direction at 50 m")

ggsave("wind_polar_mast_50m.png", width = 5, height = 5)

#---------------------
w50era5 <- wind50 %>%
  select(dir50, u50.ERA5) %>%
  mutate(cut50 = cut(u50.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w50era5$cut50,w50era5$dir50))

# Frequency distribution
kable(t(apply(table(w50era5$cut50,w50era5$dir50),2,sum)))

# Plotting
ggplot(w50era5, aes(x =factor(dir50),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 50 m",
       subtitle = "Based on the mast wind direction at 50 m")

ggsave("wind_polar_era5_50m.png", width = 5, height = 5)

#---------------------
w63mast <- wind63 %>%
  select(dir61, u63.MAST) %>%
  mutate(cut50 = cut(u63.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w63mast$cut50,w63mast$dir61))

# Frequency distribution
kable(t(apply(table(w63mast$cut50,w63mast$dir61),2,sum)))

# Plotting
ggplot(w63mast, aes(x =factor(dir61),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title = "Mast obs. at 63 m",
       subtitle = "Based on the mast wind direction at 61 m")

ggsave("wind_polar_mast_63m.png", width = 5, height = 5)

#---------------------
w63era5 <- wind63 %>%
  select(dir61, u63.ERA5) %>%
  mutate(cut50 = cut(u63.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w63era5$cut50,w63era5$dir61))

# Frequency distribution
kable(t(apply(table(w63era5$cut50,w63era5$dir61),2,sum)))

# Plotting
ggplot(w63era5, aes(x =factor(dir61),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 63 m",
       subtitle = "Based on the mast wind direction at 61 m")

ggsave("wind_polar_era5_63m.png", width = 5, height = 5)

rm(list = ls(all.names = TRUE))
setwd("~/RWorkspace/")

#load libraries 
library(tidyverse)
library(lubridate)
#library(geosphere)
library(tdr)
library(lattice)
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
  #filter(flag_SAT==0, .preserve = TRUE) %>%
  #select(Time_org, Time_UTC, SAT_15m_WDave, SAT_15m_WSave, 
  #       CUP_50m_WSave, VANE_50m_WDave, orbit) %>%
  rename(u15=SAT_15m_WSave, u25=CUP_25m_WSave, d15=SAT_15m_WDave, 
         u50=CUP_50m_WSave, d50=VANE_50m_WDave, d61=VANE_61m_WDave,
         u63=CUP_63m_WSave) %>%
  mutate(dir15 = case_when(
    d15>360 ~ "invalid",
    (d15<11.25 | d15>348.75) ~ "N",
    (d15>11.25 & d15<33.75) ~ "NNE",
    (d15>33.75 & d15<56.25) ~ "NE",
    (d15>56.25 & d15<78.75) ~ "ENE",
    (d15>78.75 & d15<101.25) ~ "E",
    (d15>101.25 & d15<123.75) ~ "ESE",
    (d15>123.75 & d15<146.25) ~ "SE",
    (d15>146.25 & d15<168.75) ~ "SSE",
    (d15>168.75 & d15<191.25) ~ "S",
    (d15>191.25 & d15<213.75) ~ "SSW",
    (d15>213.75 & d15<236.25) ~ "SW",
    (d15>236.25 & d15<258.75) ~ "WSW",
    (d15>258.75 & d15<281.25) ~ "W",
    (d15>281.25 & d15<303.75) ~ "WNW",
    (d15>303.75 & d15<326.25) ~ "NW",
    TRUE ~ "NNW")) %>%
  mutate(dir50 = case_when(
    d50>360 ~ "invalid",
    (d50<11.25 | d50>348.75) ~ "N",
    (d50>11.25 & d50<33.75) ~ "NNE",
    (d50>33.75 & d50<56.25) ~ "NE",
    (d50>56.25 & d50<78.75) ~ "ENE",
    (d50>78.75 & d50<101.25) ~ "E",
    (d50>101.25 & d50<123.75) ~ "ESE",
    (d50>123.75 & d50<146.25) ~ "SE",
    (d50>146.25 & d50<168.75) ~ "SSE",
    (d50>168.75 & d50<191.25) ~ "S",
    (d50>191.25 & d50<213.75) ~ "SSW",
    (d50>213.75 & d50<236.25) ~ "SW",
    (d50>236.25 & d50<258.75) ~ "WSW",
    (d50>258.75 & d50<281.25) ~ "W",
    (d50>281.25 & d50<303.75) ~ "WNW",
    (d50>303.75 & d50<326.25) ~ "NW",
    TRUE ~ "NNW")) %>%
  mutate(dir61 = case_when(
    d61>360 ~ "invalid",
    (d61<11.25 | d61>348.75) ~ "N",
    (d61>11.25 & d61<33.75) ~ "NNE",
    (d61>33.75 & d61<56.25) ~ "NE",
    (d61>56.25 & d61<78.75) ~ "ENE",
    (d61>78.75 & d61<101.25) ~ "E",
    (d61>101.25 & d61<123.75) ~ "ESE",
    (d61>123.75 & d61<146.25) ~ "SE",
    (d61>146.25 & d61<168.75) ~ "SSE",
    (d61>168.75 & d61<191.25) ~ "S",
    (d61>191.25 & d61<213.75) ~ "SSW",
    (d61>213.75 & d61<236.25) ~ "SW",
    (d61>236.25 & d61<258.75) ~ "WSW",
    (d61>258.75 & d61<281.25) ~ "W",
    (d61>281.25 & d61<303.75) ~ "WNW",
    (d61>303.75 & d61<326.25) ~ "NW",
    TRUE ~ "NNW")) 

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

#-------------------------------------------
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
         u63=u10*(log((z63-d)/z0)/log((h-d)/z0)))    

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
wind15 <- tbl_tot %>%
  filter(u15.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir15) #%>%
#  summarise(n=n(), u15ave.MAST=mean(u15.MAST), 
#           u10ave.ERA5=mean(u10.ERA5), u15ave.ERA5=mean(u15.ERA5),
#           u15med.MAST=median(u15.MAST), 
#           u10med.ERA5=median(u10.ERA5), u15med.ERA5=median(u15.ERA5))

wind50 <- tbl_tot %>%
  filter(u50.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir50) #%>%
#  summarise(n=n(), u50ave.MAST=mean(u50.MAST),u50ave.ERA5=mean(u50.ERA5),
#            u50med.MAST=median(u50.MAST), u50med.ERA5=median(u50.ERA5))

wind63 <- tbl_tot %>%
  filter(u63.MAST!=999.99, .preserve = TRUE) %>%
  group_by(dir61) #%>%
#  summarise(n=n(), u63ave.MAST=mean(u63.MAST),u63ave.ERA5=mean(u63.ERA5),
#            u63med.MAST=median(u63.MAST), u63med.ERA5=median(u63.ERA5))

#-----------------
ord<-c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW")

wind15$dir15<- factor(wind15$dir15, levels=ord)
wind50$dir50<- factor(wind50$dir50, levels=ord)
wind63$dir61<- factor(wind63$dir61, levels=ord)

#-----------------
w15mast <- wind15 %>%
  select(dir15, u15.MAST) %>%
  mutate(cut15 = cut(u15.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))
#w15mast$cut<-cut(w15mast$u15.MAST,breaks=c(0,3,6,10,20,100),labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  ") ,right =TRUE, include.lowest = TRUE) 

# Wind speed distribution
kable(table(w15mast$cut15,w15mast$dir15))

# Frequency distribution
kable(t(apply(table(w15mast$cut15,w15mast$dir15),2,sum)))

# Plotting
setwd("./figures/")

ggplot(w15mast, aes(x =factor(dir15),fill = cut15)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="Mast obs. at 15 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_mast_15m.png", width = 5, height = 5)

#---------------------
w10era5 <- wind15 %>%
  select(dir15, u10.ERA5) %>%
  mutate(cut10 = cut(u10.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w10era5$cut10,w10era5$dir15))

# Frequency distribution
kable(t(apply(table(w10era5$cut10,w10era5$dir15),2,sum)))

# Plotting
ggplot(w10era5, aes(x =factor(dir15),fill = cut10)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 10 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_era5_10m.png", width = 5, height = 5)

#---------------------
w15era5 <- wind15 %>%
  select(dir15, u15.ERA5) %>%
  mutate(cut10 = cut(u15.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w15era5$cut10,w15era5$dir15))

# Frequency distribution
kable(t(apply(table(w15era5$cut10,w15era5$dir15),2,sum)))

# Plotting
ggplot(w15era5, aes(x =factor(dir15),fill = cut10)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 15 m",
       subtitle = "Based on the mast wind direction at 15 m")

ggsave("wind_polar_era5_15m.png", width = 5, height = 5)

#---------------------
w50mast <- wind50 %>%
  select(dir50, u50.MAST) %>%
  mutate(cut50 = cut(u50.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w50mast$cut50,w50mast$dir50))

# Frequency distribution
kable(t(apply(table(w50mast$cut50,w50mast$dir50),2,sum)))

# Plotting
ggplot(w50mast, aes(x =factor(dir50),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="Mast obs. at 50 m",
       subtitle = "Based on the mast wind direction at 50 m")

ggsave("wind_polar_mast_50m.png", width = 5, height = 5)

#---------------------
w50era5 <- wind50 %>%
  select(dir50, u50.ERA5) %>%
  mutate(cut50 = cut(u50.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w50era5$cut50,w50era5$dir50))

# Frequency distribution
kable(t(apply(table(w50era5$cut50,w50era5$dir50),2,sum)))

# Plotting
ggplot(w50era5, aes(x =factor(dir50),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 50 m",
       subtitle = "Based on the mast wind direction at 50 m")

ggsave("wind_polar_era5_50m.png", width = 5, height = 5)

#---------------------
w63mast <- wind63 %>%
  select(dir61, u63.MAST) %>%
  mutate(cut50 = cut(u63.MAST,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w63mast$cut50,w63mast$dir61))

# Frequency distribution
kable(t(apply(table(w63mast$cut50,w63mast$dir61),2,sum)))

# Plotting
ggplot(w63mast, aes(x =factor(dir61),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title = "Mast obs. at 63 m",
       subtitle = "Based on the mast wind direction at 61 m")

ggsave("wind_polar_mast_63m.png", width = 5, height = 5)

#---------------------
w63era5 <- wind63 %>%
  select(dir61, u63.ERA5) %>%
  mutate(cut50 = cut(u63.ERA5,breaks=c(0,3,6,10,20,100),
                     labels=c("[0,3]","(3,6]","(6,10]","(10,20]" ,"(20,  "),
                     right =TRUE, include.lowest = TRUE))

# Wind speed distribution
kable(table(w63era5$cut50,w63era5$dir61))

# Frequency distribution
kable(t(apply(table(w63era5$cut50,w63era5$dir61),2,sum)))

# Plotting
ggplot(w63era5, aes(x =factor(dir61),fill = cut50)) + 
  geom_bar(width =0.8) + 
  coord_polar(theta = "x",start = -0.2) + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom")+
  labs(fill = "Wind speed [m/s]", x="Direction", y="",
       title="SAR wind at 63 m",
       subtitle = "Based on the mast wind direction at 61 m")

ggsave("wind_polar_era5_63m.png", width = 5, height = 5)

