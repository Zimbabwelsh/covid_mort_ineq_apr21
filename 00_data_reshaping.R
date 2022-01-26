## Recoding MSOA death counts

## while(!try(require(tidyverse), silent = TRUE)) install.packages("tidyverse")
## while(!try(require(readxl), silent = TRUE)) install.packages("readxl")  

suppressMessages(suppressPackageStartupMessages({
library(tidyverse)
library(readxl)
library(reshape2)
library(sf)
}))

setwd('path/to/data')

#Read in mortality counts data - available from https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsduetocovid19bylocalareaanddeprivation
df <- read_xlsx("MSOAmort_toApr2021.xlsx", sheet=10, skip=5)

#Specify months of observation
obs_months <- 14
#Enumerate months from Jan 2020
months <- c((1:(obs_months))+2)

#Drop weird descriptive final ONS rows
df <- df[c(1:7201),]

#list weird empty ONS columns
dropcols <- c(4,(4+obs_months+2),(4+2*(obs_months+2)))
df <- df[-dropcols]

#Drop total mortality columns - we can calculate from total of COV/NONCOV if we need later.
totcols <- c((1:obs_months)+3)
df <- df[-totcols]

#Drop summed death columns
sumcols <- c(0,obs_months+1, (obs_months+1)*2)+4
df <- df[-sumcols]

#generate colnames
covs <- c()
ncovs <- c()

for(i in months){
  covs[i-2] <- paste0("cov", i)
  ncovs[i-2] <- paste0("ncov", i)
}

names <- c("cd", "ons", "nm", covs, ncovs)
colnames(df) <- names

#ditch cols not needed for this analysis
df <- df[,-c(2,3)]

#melt data for single row per observation, can group by variables otherswise
mdf <- melt(df, id="cd")

## Read age data
age <- read_xlsx("2018_pop_ests.xlsx", sheet =4, skip=4)
# Create age-band proportions <25, 25-44, 45-65, 65-75, 75+
age$`<25` <- rowSums(age[,4:8])
age$`25-44` <- rowSums(age[,9:12])
age$`45-64` <- rowSums(age[,13:16])
age$`65-74` <- rowSums(age[,17:18])
age$`75+` <- rowSums(age[19:22])

# Generate proportion variables
age$`<25prop` <- age$`<25`/age$`All Ages`
age$`25-44prop` <- age$`25-44`/age$`All Ages`
age$`45-64prop` <- age$`45-64`/age$`All Ages`
age$`65-74prop` <- age$`65-74`/age$`All Ages`
age$`75+prop` <- age$`75+`/age$`All Ages`

# drop constituent cols
age <- age[,-c(2,4:27)]

# Copy age profile for each MSOA 
mort <- left_join(mdf, age, by=c("cd"="Area Codes"))

mort <- mort %>% group_by(variable) %>% 
  ### first find total population
  mutate(total_pop = sum(`All Ages`),
         ### then total deaths
         total_death = sum(as.numeric(value)),
         ### expected value is total deaths divided by total population times population of msoa
         expected = (total_death/total_pop)*`All Ages`,
         ### offset is log of this value
         offset = log(expected))

#######
# Integrate area characteristics with shapefiles.
#######

### read in file linking lsoa to stp
# available from https://geoportal.statistics.gov.uk/datasets/520e9cd294c84dfaaf97cc91494237ac_0
lsoa_ccg <- read_csv("LSOA_CCG.csv") %>% 
  select(LSOA11CD, LSOA11NM, LAD19CD, LAD19NM, STP19CD, STP19NM)
### stps are for england for wales there are health boards
# available from https://data.england.nhs.uk/dataset/ods-welsh-local-health-boards-and-sites/resource/464881dd-81fc-4cbe-9681-da24971fafda
wales_healthboard <- read_csv("wales_lad_healthboard.csv") %>% 
  rename("STP19CD" = "LHB19CD",
         "STP19NM" = "LHB19NM",
         "LAD17CD" = "UA19CD",
         "LAD17NM" = "UA19NM") %>% 
  select(-FID)
area_lookups <- read_csv("OAtoLSOAtoMSOAtoLAD.csv") %>% 
  distinct(LSOA11CD, .keep_all = T) %>% 
  select(LSOA11CD, LSOA11NM, LAD17CD, LAD17NM)

wales_healthboard <- inner_join(wales_healthboard, area_lookups, by = c("LAD17CD", "LAD17NM")) %>% 
  rename("LAD19CD" = "LAD17CD",
         "LAD19NM" = "LAD17NM")
lsoa_ccg <- bind_rows(lsoa_ccg,wales_healthboard)

#### read in the IMD data and join to stp data
# available from https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/833970/File_1_-_IMD2019_Index_of_Multiple_Deprivation.xlsx 
IMD <- read_csv("IMD.csv")
lsoa_ccg <- left_join(lsoa_ccg, IMD, by = c("LSOA11CD" = "LSOA code (2011)"))

### read in UK IMD and also join
# available from https://data.bris.ac.uk/data/dataset/1ef3q32gybk001v77c1ifmty7x
UKIMD <- read_csv("UK_IMD_E.csv") %>%
  select(lsoa, UK_IMD_E_score)
lsoa_ccg <- left_join(lsoa_ccg, UKIMD, by = c("LSOA11CD" = "lsoa"))

### create stp level IMD
lsoa_ccg <- lsoa_ccg %>% 
  group_by(STP19CD) %>% 
  mutate(avgIMDstp = mean(`Index of Multiple Deprivation (IMD) Rank`),
         avgUKIMDstp = mean(UK_IMD_E_score))


### a file of imd data just for stp
stp_imd <- lsoa_ccg %>% 
  select(STP19CD, avgIMDstp, avgUKIMDstp) %>% 
  distinct(STP19CD, .keep_all = TRUE)


### Generate key MSOA file to link other area characteristics to
### read in file linking msoa to lsoa
# Master linking file available at https://geoportal.statistics.gov.uk/datasets/ons::oa-to-lad-to-lsoa-to-msoa-to-lep-overlapping-part-april-2020-lookup-in-england-partial-coverage/about
msoa_lsoa <- read_csv("OAtoLSOAtoMSOAtoLAD.csv") %>% 
  distinct(LSOA11CD, .keep_all = T) %>% 
  select(LSOA11CD,MSOA11CD, MSOA11NM,LAD17CD, LAD17NM, RGN11CD,RGN11NM)


### join ttwas
ttwas <- read_csv("LSOA2011toTTWA2011.csv")
ttwas <- ttwas %>% select(LSOA11CD, TTWA11CD, TTWA11NM)
msoa_lsoa <- left_join(msoa_lsoa, ttwas, by = "LSOA11CD")

### join with deprivation data
msoa_lsoa <- left_join(msoa_lsoa, IMD, by = c("LSOA11CD" = "LSOA code (2011)"))
msoa_lsoa <- left_join(msoa_lsoa, UKIMD, by = c("LSOA11CD" = "lsoa"))

### decompose the deprivation data into structural mean for between-within formulation
msoa_lsoa <- msoa_lsoa %>% 
  group_by(RGN11CD) %>% 
  mutate(avgIMDreg = mean(`Index of Multiple Deprivation (IMD) Rank`),
         avgUKIMDreg = mean(UK_IMD_E_score)) %>% 
  ungroup() %>% 
  group_by(LAD17CD) %>% 
  mutate(avgIMDlad = mean(`Index of Multiple Deprivation (IMD) Rank`),
         avgUKIMDlad = mean(UK_IMD_E_score)) %>%
  ungroup() %>% 
  group_by(TTWA11CD) %>% 
  mutate(avgIMDttwa = mean(`Index of Multiple Deprivation (IMD) Rank`, na.rm =T),
         avgUKIMDttwa = mean(UK_IMD_E_score)) %>%
  ungroup() %>% 
  group_by(MSOA11CD) %>% 
  mutate(avgIMDmsoa = mean(`Index of Multiple Deprivation (IMD) Rank`),
         avgUKIMDmsoa = mean(UK_IMD_E_score))


### join two together usin lsoa as key
msoa_stp <- left_join(msoa_lsoa, lsoa_ccg, by = "LSOA11CD")

### there is a slight problem here because 17 MSOAs are in different STPs
### for the purpose of this analysis we will assign each MSOA to the STP/
### within within which the majority of the lsoas within that msoa are part of
### this code does this and matches each MSOA with exactly one STP
msoa_stp <- msoa_stp %>% 
  group_by(MSOA11CD) %>% 
  count(STP19CD, STP19NM) %>% 
  ungroup() %>% 
  arrange(MSOA11CD,-n) %>% 
  distinct(MSOA11CD, .keep_all = T)

msoa_stp <- left_join(msoa_stp, stp_imd, by = "STP19CD")

### Rejoin STP information with MSOA data to generate overall area-information list
areas <- left_join(msoa_lsoa,msoa_stp, by = "MSOA11CD") %>% 
  distinct(MSOA11CD, .keep_all = T) %>% 
  select(MSOA11CD, MSOA11NM, RGN11CD, RGN11NM,STP19CD, STP19NM, LAD17CD, LAD17NM, TTWA11CD, TTWA11NM, 
         avgIMDreg, avgIMDlad, avgIMDmsoa, avgIMDstp, avgIMDttwa, avgUKIMDreg, avgUKIMDlad, avgUKIMDmsoa, avgUKIMDstp, avgUKIMDttwa) %>% 
  filter(RGN11NM != "Scotland")


### Merge with MSOA death counts data
# mort <- read_csv("MSOAmort.csv")
mortareas <- left_join(mort,areas, by = c("cd"  = "MSOA11CD"))



###############################
# Care Homes spatial matching #
###############################

#### add care home data (available from https://covid19.esriuk.com/datasets/e4ffa672880a4facaab717dea3cdc404_0)
care_homes <- st_read("Geolytix_Carehomes/Geolytix_UK_Care_Homes_2020.shp")
#### transform to british national grid to match with msoa data
care_homes <- st_transform(care_homes,27700)

### read in msoa shapefile (available from https://geoportal.statistics.gov.uk/datasets/826dc85fb600440889480f4d9dbb1a24_0)
msoa <- st_read("MSOA_shapefiles/Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries_Generalised_Clipped_(BGC)_EW_V3.shp")
### join each care home with an msoa
data <- st_join(care_homes, msoa, join = st_within)

### count how many care home are in each msoa
care_home_msoa <- data %>% 
  as_tibble() %>% 
  group_by(MSOA11CD) %>% 
  tally()

#### join count in each msoa back to msoa data
care_home_msoa <- left_join(msoa, care_home_msoa, by = "MSOA11CD") %>% 
  as_tibble() %>% 
  select(MSOA11CD,n) %>% 
  rename(care_homes = n) %>% 
  ### msoas with no care homes currently read NA so replace this with 0
  replace(is.na(.), 0)

#### join this with the rest of the data
mortareas <- left_join(mortareas, care_home_msoa, by = c("cd" = "MSOA11CD"))


##########################################
# Write out final data file for analyses #
##########################################

### write to csv
write_csv(mortareas, "model_data_toapril.csv")

