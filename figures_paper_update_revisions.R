suppressMessages(suppressPackageStartupMessages({
library(data.table)
library(ggplot2)
library(tidyverse)
library(corrplot)
library(gridExtra)
library(viridis)
library(lubridate)
}))

setwd("path/to/data")

# Set names for random coefficients, here months post March 2020
rcoefnames <- c("Mar20", "Apr20", "May20", "Jun20", "July20", "Aug20", "Sep20", "Oct20", "Nov20", "Dec20", "Jan21", "Feb21", "Mar21", "Apr21")
# Set names for structural levels in the model
levnames <- c("MSOA","TTWA", "Region")

## Read in using baseline model output from MLwiN (extracted to csv from c1090) - uses bespoke "data_process" 
## function, written out at the end
# Produces MRRs, Correlations, and estimates of Total Variance
COV_noIMD <- data_process("cov_noimd_500k.csv", 500000, 3, c("MSOA", "TTWA", "Region"), 32, 14, rcoefnames, "COVID-19")
covMRRs <- COV_noIMD[[1]]
covcorrs <- COV_noIMD[[2]]
totalvars <- COV_noIMD[[3]]

NCOV_noIMD <- data_process("ncov_noimd_500k.csv", 500000, 3, c("MSOA", "TTWA", "Region"), 32, 14, rcoefnames, "Non COVID-19")
ncovMRRS <- NCOV_noIMD[[1]]
ncovcorrs <- NCOV_noIMD[[2]]

#Create MRR dataset for both COVID and non-COVID mortality
MRR_noIMD <- bind_rows(list(covMRRs, ncovMRRS))

# Convert month indicators from model into dates (mid-month)
date <- dmy("11/12/2019")
MRR_noIMD$date <- ymd(date %m+% months(as.numeric(MRR_noIMD$Month)))
#Create outcome indicator
MRR_noIMD$outcome <- factor(MRR_noIMD$outcome, levels=c("COVID-19", "Non COVID-19"))

#Create correlations data-frame
corr_noIMD <- bind_rows(list(covcorrs, ncovcorrs))


#Repeat for models adjusted for IMD and population density
COV_IMD <- data_process("cov_ukimd_popden_500k.csv", 500000, 3, c("MSOA", "TTWA", "Region"), 74, 14, rcoefnames, "COVID-19")
covMRRs <- COV_IMD[[1]]
covcorrs <- COV_IMD[[2]]

NCOV_IMD <- data_process("ncov_ukimd_popden_500k.csv", 500000, 3, c("MSOA", "TTWA", "Region"), 74, 14, rcoefnames, "Non COVID-19")
ncovMRRs <- NCOV_IMD[[1]]
ncovcorrs <- NCOV_IMD[[2]]

# Create MRR dataset for COVID and nonCOVID mortality, adjusted for IMD and pop-density
MRR_UKIMD <- bind_rows(list(covMRRs, ncovMRRs))
MRR_UKIMD$date <- ymd(date %m+% months(as.numeric(MRR_UKIMD$Month)))
corr_UKIMD <- bind_rows(list(covcorrs, ncovcorrs))
#Create outcome indicator
MRR_UKIMD$outcome <- factor(MRR_UKIMD$outcome, levels=c("COVID-19", "Non COVID-19"))


#Read study dates

dates <- read.csv("../../studyDates.csv")
dates$StartDate <- dates$StartDate %>% dmy()
dates$EndDate <- dates$EndDate %>% dmy()
colnames <- c("Period", "Start", "End")
colnames(dates) <- colnames
dates$Period <- factor(dates$Period, levels=c(1,2,3,4))

#########
# PLOTS #
#########

########################
# Deaths timeline plot #
########################

#Read in daily deaths - accessible at https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/articles/coronaviruscovid19latestinsights/deaths
daily <- read.csv("../../../nation_daily_deaths_toOct21.csv")
#Create lubridate interval for study period
mort_period <- interval(ymd("2020-03-01"), ymd("2021-04-30")) 
# Select relevant columns, dates within interval, and only deaths in England and Wales
daily$date <- daily$date %>% dmy()
daily <- daily %>% select(areaName, date, newDeaths28DaysByDeathDate) %>% filter(areaName=="England"|areaName=="Wales")
daily <- daily %>% filter(ymd(date) %within% mort_period)
daily <- daily %>% pivot_wider(names_from = areaName, values_from = newDeaths28DaysByDeathDate)
daily$Wales <- replace_na(daily$Wales, 0)
daily <- daily %>% mutate(total=England+Wales) %>% select(date, total)


# Add study dates for annotation
dates <- read.csv("../../studyDates.csv")
dates$StartDate <- dates$StartDate %>% dmy()
dates$EndDate <- dates$EndDate %>% dmy()
colnames <- c("Period", "Start", "End")
colnames(dates) <- colnames
dates$Period <- factor(dates$Period, levels=c(1,2,3,4))



#####################
# Generate Figure 1 #
#####################

deathplot <- ggplot(daily,
       aes(x=date, y=total))+
  geom_ribbon(aes(ymin=0, ymax=total),  fill="dark grey")+
  ylab("Daily COVID-19 deaths in England and Wales")+
  xlab("Date")+
  theme(legend.position = "none")
  
# Annotation text
ann.text <- data.frame(Start = ymd(c("2020-10-20", "2020-11-03")), 
                       y=c(-112.5, -40),
                       lab= c("Tiered Lockdowns", "National")) 
# Annotation boxes
annot_death_plot <- deathplot+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=-150, ymax=0, alpha=0.6, fill="#7570b3")+
  annotate("rect", xmin=ymd("2020-11-05"), xmax=ymd("2020-12-02"), ymin=-75, ymax=0, alpha=0.5, fill="#1d449f")+
  annotate("rect", xmin=dates$Start[3], xmax=ymd("2020-11-05"), ymin=0, ymax=Inf, alpha =0.3, fill="#7570b3")+
  annotate("rect", xmin=ymd("2020-12-02"), xmax=dates$End[3], ymin=0, ymax=Inf, alpha =0.3, fill="#7570b3")+
  annotate("rect", xmin=ymd("2020-11-05"), xmax=ymd("2020-12-02"), ymin=0, ymax=Inf, alpha=0.4, fill = "#1d449f")+
  geom_text(data=ann.text,
            mapping = aes(x=Start, y=y, label = lab),
            size=2.7,
            hjust=-0.05,
            vjust=0.05,
            colour="white")+
  geom_vline(xintercept = c(dates$End[2:3]), linetype=2, alpha=0.6)
annot_death_plot

ggsave("dateDeathPlot.png", dpi=350)




##########################
### Combined MRR plots ###
##########################

# Produce ggplot graphic
NOIMD_MRRplot <- ggplot(MRR_noIMD, aes(x=date, y=Median, colour=level))+
  geom_line(size=0.8)+
  scale_colour_manual(values=c("#1b9e77","#d95f02", "#7570b3"))+
  ylab("Median Rate Ratio (MRR)")+
  scale_y_continuous(trans='log2')+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8)) +
  scale_x_date(date_labels = "%b %y")+
  geom_ribbon(aes(ymin=Lower, ymax=Upper), linetype=3, alpha=0.2) +
  geom_vline(xintercept=13, linetype=2, alpha=0.5, col="black")+
  facet_grid(outcome~level)+
  #theme_bw()+
  theme(strip.background = element_blank())+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0, ymax=Inf, alpha=0.15, fill="#7570b3")

NOIMD_MRRplot

ggsave("NOIMD_MRRplot.png", dpi=350)

UKIMD_MRRplot <- ggplot(MRR_UKIMD, aes(x=date, y=Median, colour=level))+
  geom_line(size=0.8)+
  scale_colour_manual(values=c("#1b9e77","#d95f02", "#7570b3"))+
  ylab("Median Rate Ratio (MRR)")+
  scale_y_continuous(breaks=c(1,2,4,8,16),trans='log')+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8)) +
  scale_x_date(date_labels = "%b %y")+
  geom_ribbon(aes(ymin=Lower, ymax=Upper), linetype=3, alpha=0.2) +
  geom_vline(xintercept=13, linetype=2, alpha=0.5, col="black")+
  facet_grid(outcome~level)+
  #theme_bw()+
  theme(strip.background = element_blank())+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0, ymax=Inf, alpha=0.15, fill="#7570b3")

UKIMD_MRRplot

ggsave("UKIMD_MRRplot.png", dpi=350)

######################
# Combined VPC plots #
######################
NOIMD_VPCplot <- ggplot(data=MRR_noIMD, aes(x=date, y=VPC, colour=level)) +
  geom_line(size=1)+
  ylim(0,1) +
  ylab("Variance Partitioning Coefficient") +
  facet_grid(outcome~level)+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8))+
  scale_x_date(date_labels = "%b %y")+
  theme(strip.background = element_blank())+
  scale_colour_manual(values=c("#1b9e77","#d95f02", "#7570b3"))+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0, ymax=Inf, alpha=0.15, fill="#7570b3")

NOIMD_VPCplot

ggsave("NOIMD_VPCplot.png", dpi=350)

UKIMD_VPCplot <- ggplot(data=MRR_UKIMD, aes(x=date, y=VPC, colour=level)) +
  geom_line(size=1)+
  ylim(0,1) +
  ylab("Variance Partitioning Coefficient") +
  facet_grid(outcome~level)+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8))+
  scale_x_date(date_labels = "%b %y")+
  theme(strip.background = element_blank())+
  scale_colour_manual(values=c("#1b9e77","#d95f02", "#7570b3"))+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0, ymax=Inf, alpha=0.15, fill="#7570b3")

UKIMD_VPCplot

ggsave("UKIMD_VPCplot.png", dpi=350)

###########################
### Combined Corr plots ###
###########################

### combine covid and other deaths together
correstsMod2 <-  bind_rows(list(covcorrs, ncovcorrs))

### make the plot
ggplot(data = corr_noIMD,
       aes(month1,
           fct_rev(month2), 
           fill=Median, 
           label=sprintf("%.2f", round(Median,2)))) +
  geom_tile() +
  labs(x = NULL,
       y = NULL,
       fill = "Pearson's\nCorrelation", 
       title="Multiscale monthly mortality correlations") + 
  ### add colour scale
  scale_fill_gradient2(mid="#FBFEF9",low="#A63446",high="#0C6291", limits=c(-1,1)) +
  ### add in the numbers on the tiles
  # geom_text() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_minimal() +
  ### change text sizes
  theme(axis.text=element_text(size=8, colour = "black"),
        axis.text.x=element_text(angle=90),
        strip.text = element_text(size=12),
        legend.title=element_text(size=12),
        plot.title = element_text(size = 16)) +
  ### create 2 by 4 plot with each level as row starting with region (fct_rev)
  facet_grid(fct_rev(level)~outcome) 

ggplot(data = corr_UKIMD,
       aes(month1,
           fct_rev(month2), 
           fill=Median, 
           label=sprintf("%.2f", round(Median,2)))) +
  geom_tile() +
  labs(x = NULL,
       y = NULL,
       fill = "Pearson's\nCorrelation", 
       title="Multiscale monthly mortality correlations") + 
  ### add colour scale
  scale_fill_gradient2(mid="#FBFEF9",low="#A63446",high="#0C6291", limits=c(-1,1)) +
  ### add in the numbers on the tiles
  # geom_text() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_minimal() +
  ### change text sizes
  theme(axis.text=element_text(size=8, colour = "black"),
        axis.text.x=element_text(angle=90),
        strip.text = element_text(size=12),
        legend.title=element_text(size=12),
        plot.title = element_text(size = 16)) +
  ### create 2 by 4 plot with each level as row starting with region (fct_rev)
  facet_grid(fct_rev(level)~outcome) 



ggplot(corr_noIMD,
       aes(month1, month2, fill=Median))+
  geom_raster(hjust = 0, vjust = 0)+
  facet_grid(fct_rev(level)~outcome) +
  scale_color_viridis()+
  theme(axis.text.x=element_text(size=10, angle=90))

ggsave("Mod2ENGIMDMortCorrs.tiff")

########################
# UK COVID deaths data #
########################


#######################
# Fixed effects plots #
#######################
### Data reading
MCMC <- fread("cov_ukimd_popden_500k.csv")

# Read in & Reshape data
df <- as.data.frame(matrix(MCMC$V1, 500000, (nrow(MCMC)/500000), byrow=TRUE))
# Drop final col as it only contains level 1 (constant)
df <- df[,-(ncol(df))]
raw_fixed_ests <- df[,15:74]

cov_density <- df[,75:88]

# Quants for percentiles
quants <- c(0.025,0.05,0.50,0.95,0.975)
# Calculate percentiles of fixed effects and reshape, rename
fixed_ests <- as.data.frame(apply(raw_fixed_ests, 2, quantile, probs = quants))
cov_density <- as.data.frame(apply(density, 2, quantile, probs=quants))
colnames <- c("25-44", "45-64", "65-74", "75+",
              paste("CH", rcoefnames, sep=""),
              paste("RegIMD", rcoefnames, sep=""),
              paste("M-TIMD", rcoefnames, sep=""),
              paste("T_RIMD", rcoefnames, sep=""))
colnames(fixed_ests) <- colnames

fixed_ests <-  fixed_ests %>% select(contains("IMD"))

fixed_ests <- fixed_ests %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)
fixed_ests$month <- rep(3:16, 3)
fixed_ests$level <- c(rep("Reg", 14), rep("MSOA", 14), rep("TTWA", 14))
fixed_ests$level <- factor(fixed_ests$level, levels=c("Reg", "TTWA", "MSOA"))


MCMC <- fread("ncov_ukimd_popden_500k.csv")

# Read in & Reshape data
df <- as.data.frame(matrix(MCMC$V1, 500000, (nrow(MCMC)/500000), byrow=TRUE))
# Drop final col as it only contains level 1 (constant)
df <- df[,-(ncol(df))]
raw_fixed_ests <- df[,15:74]
ncov_density <- df[,75:88]

# Quants for percentiles
quants <- c(0.025,0.05,0.50,0.95,0.975)
# Calculate percentiles of fixed effects and reshape, rename
ncov_fixed_ests <- as.data.frame(apply(raw_fixed_ests, 2, quantile, probs = quants))

colnames <- c("25-44", "45-64", "65-74", "75+",
              paste("CH", rcoefnames, sep=""),
              paste("RegIMD", rcoefnames, sep=""),
              paste("T_RIMD", rcoefnames, sep=""),
              paste("M-TIMD", rcoefnames, sep=""))
colnames(ncov_fixed_ests) <- colnames

ncov_fixed_ests <-  ncov_fixed_ests %>% select(contains("IMD"))

ncov_fixed_ests <- ncov_fixed_ests %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)
ncov_fixed_ests$month <- rep(3:16, 3)
ncov_fixed_ests$level <- c(rep("Reg", 14), rep("TTWA", 14), rep("MSOA", 14))
ncov_fixed_ests$level <- factor(ncov_fixed_ests$level, levels=c("Reg", "TTWA", "MSOA"))

#######################
# Combined data frame #
#######################

density_tot <- bind_cols(cov_density, ncov_density)
colnames(density_tot) <- c(paste("dens_cov_", rcoefnames, sep=""),
                           paste("dens_ncov_", rcoefnames, sep=""))

density_tot <- as.data.frame(apply(density_tot, 2, quantile, probs=quants))

density_tot <- density_tot %>% rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname,  values_from=value)
density_tot$month <- rep(3:16, 2)
density_tot$date <- date %m+% months(as.numeric(density_tot$month))
density_tot$outcome <- c(rep("COVID-19", 14), rep("Non COVID-19", 14))
density_tot <- density_tot %>% mutate(densityOR = exp(density_tot$`50%`))

ggplot(density_tot, 
       aes(x=date, y=densityOR, 
           group=outcome, 
           colour=outcome))+
  geom_line()+
  ylim(0.90, 2.5)



fixed_tot <- bind_rows(fixed_ests, ncov_fixed_ests, .id="id")
fixed_tot$type <- fixed_tot$id %>% recode(`1`="COVID-19", `2`="Non COVID-19")

date <- dmy("15/12/2019")
fixed_tot$date <- date %m+% months(as.numeric(fixed_tot$month))

### Update to convert units to SDs for each UKIMD metric

fixed_tot <- fixed_tot %>% mutate(
  SD = case_when(
    level == "Reg" ~ 4.446,
    level == "TTWA" ~ 6.171,
    level == "MSOA" ~ 13.176))

fixed_tot <- fixed_tot %>% mutate(SDOR = exp(fixed_tot$`50%`*fixed_tot$SD),
                                  SD2.5=exp(fixed_tot$`2.5%`*fixed_tot$SD),
                                  SD97.5=exp(fixed_tot$`97.5%`*fixed_tot$SD),
                                  OR = exp(fixed_tot$`50%`))


##################
# Total Variance #
##################
MRR_noIMD$variance <- (log(MRR_noIMD$Median)/0.6745)^2/2

## Create new dataframe & sum over months, collapse rows
covTotVar <- MRR_noIMD %>%
  filter(outcome=="COVID-19") %>% 
  select("outcome", "date","variance") %>%
  group_by(date) %>%
  summarise(sum(variance))

# Convert variance to a proportion
covTotVar$propmax <- covTotVar$`sum(variance)`/max(covTotVar$`sum(variance)`)
covTotVar$Measure <- "Total Variance"

# Convert deaths to a proportion
daily$propmax <- daily$total/max(daily$total)
daily$Measure <- "Total Mortality"

dailyprop <- daily %>% select(date, propmax, Measure)
varprop <- covTotVar %>% select(date, propmax, Measure)
ineq_data <- bind_rows(dailyprop, varprop)

inequality_fig <- ggplot(data=ineq_data, 
                         aes(x=date, y=propmax, group=Measure, colour=Measure))+
  geom_line(size=0.8)+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0, ymax=1, alpha=0.15, fill="#7570b3")+
  xlab("Date")+
  ylab("Proportion of Maximum")+
  scale_colour_manual(values=c("#7570b3", "#d95f02"))
  
inequality_fig

ggsave("ineq_plot.png", dpi=350)


######################
# Fixed effect Plots #
######################

IMDplot <- ggplot(data=fixed_tot, 
                  aes(x=date, y=OR,
                      group=level,colour=level))+
  geom_line(size=1.2)+
  geom_hline(yintercept = 1, colour="black", linetype=2)+
  #ylim(0.9, 1.25)+
  xlab("Date")+
  ylab("UKIMD Relative Risk Ratio")+
  # + geom_ribbon(aes(ymin=exp(`2.5%`), ymax=exp(`97.5%`)),
  # linetype=3, alpha=0.2)+
  facet_grid(vars(type))+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8))+
  scale_x_date(date_labels = "%b %y")+
  theme(strip.background = element_blank())+
  scale_colour_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0.9, ymax=1.25, alpha=0.15, fill="#7570b3")

IMDplot

ggsave('imdplot.png', dpi=300)


unstdIMDplot <- ggplot(data=fixed_tot, 
                  aes(x=date, y=SDOR,
                      group=level,colour=level))+
  geom_line(size=1.2)+
  geom_hline(yintercept = 1, colour="black", linetype=2)+
  #ylim(0.9, 1.25)+
  xlab("Date")+
  ylab("UKIMD Relative Risk Ratio")+
  # + geom_ribbon(aes(ymin=exp(`2.5%`), ymax=exp(`97.5%`)),
  # linetype=3, alpha=0.2)+
  facet_grid(vars(type))+
  theme(axis.text.x=element_text(color = "black",
                                 size=7.5, angle=30, hjust=0.8))+
  scale_x_date(date_labels = "%b %y")+
  theme(strip.background = element_blank())+
  scale_colour_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  annotate("rect", xmin=dates$Start[3], xmax=dates$End[3], ymin=0.8, ymax=2.5, alpha=0.15, fill="#7570b3")

unstdIMDplot

ggsave("IMD_combined_mort.png", dpi=350)



############################
# Data Processing Function #
############################

data_process <- function(data, iter, levels, levnames, fcoefs, rcoefs, rcoefnames, type){
  
  ##########################
  # Global data processing #
  ##########################
  
  MCMC <- fread(data)
  
  # Total coefs
  coefs <- rcoefs+fcoefs
  # Calc number of parameter estimates in iteration
  ests <- nrow(MCMC)/iter
  # Calc size of variance matrix minus diagonal
  size <- ((rcoefs*(rcoefs+1))/2)-rcoefs
  # Calc N of random estimates
  rests <- size*levels
  
  # Read in & Reshape data
  df <- as.data.frame(matrix(MCMC$V1, iter, ests, byrow=TRUE))
  # Drop final col as it only contains level 1 (constant)
  df <- df[,-(ncol(df))]
  
  ####################
  # Calculating MRRs #
  ####################
  
  ## Create variance list for omission/inclusion 
  varlist <- seq(ncol(df))
  # Mask based on MLwiN output
  # Ignore beta estimates
  varlist <- varlist[-c(seq(coefs))] 
  # Create mask for random effect variances
  variances <- vector()
  maxval=ncol(df)
  variances <- trino(rcoefs)
  basevariances <- variances
  
  # Repeat variance sequence n=level times
  repeat{  
    variances <- append(variances, basevariances+max(variances))
    if (length(variances)==(rcoefs*levels)){
      break
    }
  }
  
  # Take just variance estimates, and calculate quantiles
  quants <- c(0.025,0.05,0.50,0.95,0.975)
  
  # Dataset of solely variance iteration estimates
  rawvars <- df[,varlist[variances]]
  
  # Take quantiles from "quants" and apply to rawvars
  varests <- as.data.frame(apply(rawvars, 2, quantile, probs = quants))
  varests <- varests %>%
    tibble::rownames_to_column() %>%  
    pivot_longer(-rowname) %>% 
    pivot_wider(names_from=rowname, values_from=value)
  varests <- varests[,-1]
  # Rename Cols to something sensible
  colnames(varests) <- c("Lower", "Lower 0.05", "Median", "Upper 0.95", "Upper")
  # 
  MRRests <- exp(sqrt(2*varests)*0.6745)
  # Generate level-names variable and ensure factor to maintain order for plotting
  MRRests$level <- rev(rep(levnames, each=rcoefs))
  MRRests$level <- factor(MRRests$level, levels = levnames)
  # Generate months as continuous variable
  MRRests$Month <- rep((1:rcoefs)+2, levels)
  
  # Create Indicator for outcome-type
  MRRests$outcome <- type
  varests$level <- rev(rep(levnames, each=rcoefs))
  varests$level <- factor(varests$level, levels = levnames)
  varests$Month <- rep((1:rcoefs)+2, levels)
  
  ####################
  # Calculating VPCs #
  ####################
  
  # Take Fixed effects
  raw_fixed_ests <- df[,1:rcoefs]
  
  # Calculate percentiles of fixed effects and reshape, rename
  fixed_ests <- as.data.frame(apply(raw_fixed_ests, 2, quantile, probs = quants))
  fixed_ests <- fixed_ests %>%
    tibble::rownames_to_column() %>%  
    pivot_longer(-rowname) %>% 
    pivot_wider(names_from=rowname, values_from=value)
  fixed_ests <- fixed_ests[,-1]
  colnames(fixed_ests) <- c("Lower", "Lower 0.05", "Median", "Upper 0.95", "Upper")
  # Generate level, and select only median ests
  fixed_ests$level <- (rcoefnames)
  fixed_slim <- fixed_ests %>% select(c("Median"))
  lv3_slim <- varests[1:rcoefs,] %>% select(c("Median"))
  lv2_slim <- varests[(rcoefs+1):(2*rcoefs),] %>% select("Median")
  lv1_slim <- varests[((2*rcoefs)+1):(3*rcoefs),] %>% select("Median")
  # Take only median data for "slim" data
  VPC_data <- cbind(fixed_slim, lv3_slim, lv2_slim, lv1_slim)
  VPC_data$level <- rcoefnames
  # Update Column names
  VPC_colnames <- c("expected", rev(levnames), "month")
  colnames(VPC_data) <- VPC_colnames
  
  #Construct Marginal expectation (in average context) - taken from Austin et. al paper
  VPC_data$mu_1 <- exp(VPC_data$expected + (VPC_data$Region/2) + (VPC_data$TTWA/2) + (VPC_data$MSOA/2))
  #Construct Level-3 component of total variance
  VPC_data$lv3var <- ((VPC_data$mu_1)^2)*(exp(VPC_data$Region)-1)
  #Construct Level-2 component of total variance
  VPC_data$lv2var <- ((VPC_data$mu_1)^2)*(exp(VPC_data$Region))*(exp(VPC_data$TTWA)-1)
  #Construct level-1 component of total variance
  VPC_data$lv1var <- (VPC_data$mu_1)+((VPC_data$mu_1)^2)*exp(VPC_data$Region + VPC_data$TTWA)*(exp(VPC_data$MSOA)-1)
  # Relative contribution of Level 3 variance 
  VPC_data$VPC3 <- (VPC_data$lv3var)/(VPC_data$lv3var+VPC_data$lv2var+VPC_data$lv1var)
  #Relative contribution of Level 2 variance
  VPC_data$VPC2 <- (VPC_data$lv2var)/(VPC_data$lv3var+VPC_data$lv2var+VPC_data$lv1var)
  # Relative contribution of Level 1 variance
  VPC_data$VPC1 <- (VPC_data$lv1var)/(VPC_data$lv3var+VPC_data$lv2var+VPC_data$lv1var)
  # Total variance across all levels
  VPC_data$totalvar <- (VPC_data$mu_1)+(((VPC_data$mu_1)^2)*(exp(VPC_data$Region + VPC_data$TTWA + VPC_data$MSOA)-1))
  
  #Append to MRRests dataframe
  MRRests$VPC <- c(VPC_data$VPC3, VPC_data$VPC2, VPC_data$VPC1)
  
  #####################
  # Calculating Corrs #
  #####################
  
  ### Producing and plotting correlations 
  # Generate list of columns containing covariances
  covarlist <- varlist[-variances]
  rawcovars <- df[covarlist]
  
  # Create indices to extract from covarlist
  diag <- trino(rcoefs)
  # Covarseq1 gives the first variance location for each covariance 
  covarseq1 <- diag[cvseq(rcoefs-1)]
  
  # Add total size of v/cv matrix and repeat until list is 40 elements 
  repeat{
    covarseq1 <- append(covarseq1, tail(covarseq1, size)+(size+rcoefs))
    if (length(covarseq1)==rests){
      break
    }
  }
  
  covarseq1 <- covarseq1+coefs
  
  # Covarseq2 gives second variance location for correlation calculation
  covarseq2 <- diag[(rep(1:(rcoefs-1), 1:(rcoefs-1))+1)]
  
  # Extend using same approach as above
  repeat{
    covarseq2 <- append(covarseq2, tail(covarseq2, size)+(size+rcoefs))
    if (length(covarseq2)==rests){
      break
    }
  }
  
  # Include offset of initial fixed effects estimates
  covarseq2 <- covarseq2+coefs
  
  # Combine lists to calculate correlations for each covariance estimate and store in rawcorrs
  
  # Create empty matrix
  rawcorrs <- as.data.frame(matrix(0, iter, rests))
  
  # Apply correlation calculation (dividing by sum of relevant variances) for each covariance in matrix
  for (i in 1:length(covarlist)){
    rawcorrs[i] <- df[,covarlist[i]]/(sqrt(df[,covarseq1[i]])*sqrt(df[,covarseq2[i]]))
  }
  
  # Create correlation ests from posterior distributions of correlations
  corrests <- as.data.frame(apply(rawcorrs, 2, quantile, probs = quants))
  corrests <- corrests %>%
    tibble::rownames_to_column() %>%  
    pivot_longer(-rowname) %>% 
    pivot_wider(names_from=rowname, values_from=value)
  corrests <- corrests[,-1]
  
  # Colnames for new table
  colnames(corrests) <- c("Lower", "Lower 0.05", "Median", "Upper 0.95", "Upper")
  # Generate level-names variable and ensure factor to maintain order for plotting
  corrests$level <- rev(rep(levnames, each=size))
  corrests$level <- factor(corrests$level, levels = levnames)
  
  
  # Generate month lists for faceting
  corrests$month1 <- factor(rep(rcoefnames[cvseq(rcoefs-1)],levels), levels=rcoefnames)
  corrests$month2 <- factor(rep(rcoefnames[(rep(1:(rcoefs-1), 1:(rcoefs-1))+1)], levels), levels=rcoefnames)
  
  # create column for type of mortality
  corrests$outcome <- type
  
  ##############
  # Outputting #
  ##############
  out <- list(MRRests, corrests,varests)
  return(out)
}

# Function to create repeating ascending list (1,1,2,1,2,3...)
cvseq <- function(n){
  seq(1, n, 1) %>%
    map(~seq(1, ., 1)) %>%
    unlist()
}

# Function to create repeating ascending list (1,1,2,1,2,3...)
trino <- function(n){
  sapply(1:n, function(x) sum(1:x))
}
