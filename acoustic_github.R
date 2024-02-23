#=======================================

#Start by loading packages, remember you might need to install some if you 
# are using them for the first time

#Load packages
library(arm)
library(car)
library(effects)
library(ggplot2)
library(GGally)
library(lattice)
library(lawstat)
library(outliers)
library(tidyverse)
library(rlang)
library(gridExtra)
library(glmmTMB)
library(tidyverse)
library(DHARMa)
library(AICcmodavg)
library(performance)
library(ggpubr)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

#=======================================

#Import the data
acoustic <- read.csv(file = "acoustic_full_3.csv", 
                     header = TRUE, 
                     dec = ".", 
                     stringsAsFactors = TRUE)
str(acoustic, width = 50, strict.width = "cut")

# to make the analaysis more consistent with the descriptive part and
# "easier" to analyse I created a dataset below. It includes more general categories
# eg. not particular country but Europe, NAmerica, etc or migratory/resident
# without distinguishing into anadromous, catadromous, etc.

names(acoustic)[names(acoustic) == 'X'] <- 'analysed.sample'

# adding year
year <- data.frame(
  citation = c("Smith et al. 2021", "Smith & Smith 2021", "Johnson et al. 2022")
)

year <- year %>%
  mutate(all_years = str_extract_all(citation, "\\b\\d{4}\\b"))
acoustic$year <- str_extract(acoustic$study, "\\b\\d{4}\\b")

# adding category salmonids / non-salmonids => salmonids
acoustic <- acoustic %>%
  mutate(salmonids = case_when(
    family %in% c("Salmonidae") ~ "salmonids",
    TRUE ~ "non-salmonids"  
  ))

# combining locations into 4 categories: Europe, North America, Australia 
# (with New Zealand) and other

acoustic <- acoustic %>%
  mutate(continent = case_when(
    location %in% c("Poland", "Germany", "France", "Spain", "Netherlands",
                    "Portugal", "UK", "Norway", "Sweden", "Estonia",
                    "Scotland", "Czech Republic", "Ireland", "Denmark", 
                    "Belgium", "Romania", "Italy") ~ "Europe",
    location %in% c("USA", "Canada") ~ "North America", #removed Mexico
    location %in% c("Australia", "New Zealand", "New Caledonia", "French Polynesia") ~ "Australia",#added
    TRUE ~ "Other"  
  ))

#
acoustic <- acoustic %>%
  mutate(migration = case_when(
    migratory.status %in% c("anadromous", "catadromous", "amphidromous") ~ "diadromous",
    TRUE ~ "non-diadromous"  
  ))

# There were 3 blanks in IUCN Red list category and I didn't want to loose 
# the observations so I filled them with 
acoustic$Red.list.category[acoustic$Red.list.category == ""] <- "NE"

# removing unimportant columns (must be done after adding columns above)
acoustic <- acoustic[, !(names(acoustic) %in% c("species.common.name", "X.",
                                                "conservation.designation",
                                                "secondary.conservation","location", "fish.size..mm.",
                                                "journal", "study", "migratory.status",
                                                "life.stage", "X.1", "analysed.sample"))]

# Check the resulting dataframe
str(acoustic, width = 50, strict.width = "cut")

# 'data.frame':	578 obs. of  11 variables:
#  $ species.latin.name: Factor w/ 255 levels "","..
#  $ family            : Factor w/ 79 levels "",""..
#  $ system            : Factor w/ 2 levels "fres"..
#  $ sample.size       : int  30 40 27 24 26 5 9 8..
#  $ Red.list.category : Factor w/ 9 levels "","C"..
#  $ trophic.guild     : Factor w/ 12 levels "",""..
#  $ journal.IF        : num  2.8 2 3.2 1.6 1.6 1...
#  $ year              : chr  "2023" "2023" "2023"..
#  $ salmonids         : chr  "non-salmonids" "no"..
#  $ continent         : chr  "Australia" "Austra"..
#  $ migration         : chr  "non-diadromous" "n"..



##########################################

# DATA EXPLORATION

##########################################

# Are there missing values?
colSums(is.na(acoustic))
# 2 missing JIF

# Remove these rows
ac1 <- na.omit(acoustic)
colSums(is.na(ac1))
dim(ac1) #2 data points lost
#=======================================

# 1. OUTLIERS

# Start by defining a preferred figure format, called 'My_theme'
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1),
                  strip.background = element_rect(fill = "white", 
                                                  color = "white", linewidth = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", linewidth = 0.1),
                  panel.grid.minor = element_line(colour = "white", linewidth = 0.1))

# A function for dotplots
multi_dotplot <- function(filename, Xvar, Yvar){
  filename %>%
    ggplot(aes(x = {{Xvar}})) +
    geom_point(aes(y = {{Yvar}})) +
    My_theme +
    coord_flip() +
    labs(x = "Order of Data")}

#=======================================

# Region
box1 <- ac1 %>%
  ggplot(aes(y = sample.size, x = continent)) +
  ylim(0,200) +
  labs(y = "sample size", x = "continent") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

# Taxon
box2 <- ac1 %>%
  ggplot(aes(y = sample.size, x = salmonids)) +
  ylim(0,200) +
  labs(y = "sample size", x = "salmonids") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

# Migration
box3 <- ac1 %>%
  ggplot(aes(y = sample.size, x = migration)) +
  ylim(0,200) +
  labs(y = "sample size", x = "migration") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

# System
box4 <- ac1 %>%
  ggplot(aes(y = sample.size, x = system)) +
  ylim(0,200) +
  labs(y = "sample size", x = "system") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

# Year
box5 <- ac1 %>%
  ggplot(aes(y = sample.size, x = year)) +
  ylim(0,200) +
  labs(y = "sample size", x = "year") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

# guild
box6 <- ac1 %>%
  ggplot(aes(y = sample.size, x = trophic.guild)) +
  ylim(0,200) +
  labs(y = "sample size", x = "guild") +
  geom_boxplot(fill = "gray88", outlier.shape = NA) +
  My_theme

#Plot as a grid
grid.arrange(box1, box2, box3, 
             box4, box5, box6, nrow = 2)
# Trophic guild is not usable, but the others are fine

# A multi-panel dotchart to view all 
# continuous variables simultaneously.

#Order data to plot continuous variables
ac1 <- ac1 %>%
  mutate(order = seq(1:nrow(ac1)))

#Select continuous variables to plot
p1 <- multi_dotplot(ac1, order, sample.size)
p2 <- multi_dotplot(ac1, order, journal.IF)
 
#Plot as a grid
grid.arrange(p1, p2, nrow = 2)


# Deal with sample size first
# Remove sample size >400
ac2 <- ac1[which(ac1$sample.size < 400),]
dim(ac1)
dim(ac2)
range(ac2$sample.size)

p3 <- multi_dotplot(ac2, order, sample.size)
grid.arrange(p1, p3, nrow = 2)


# And IF
p4 <- multi_dotplot(ac2, order, journal.IF)
grid.arrange(p4, nrow = 2)


# Transform JIF with square-root
ac2$sqrt.journal.IF <- sqrt(ac2$journal.IF)

p5 <- multi_dotplot(ac2, order, sqrt.journal.IF)
grid.arrange(p4, p5, nrow = 2)


#=======================================

# 2. NORMALITY AND HOMOGENEITY OF DEPENDENT VARIABLE

# Frequency polygon plot
ac2 %>% ggplot(aes(sample.size)) +
  geom_freqpoly(bins = 7) +
  labs(x = "sample size", y = "Frequency") +
  My_theme
# Positively skewed

#Shapiro-Wilk test for deviation from normality
shapiro.test(ac2$sample.size)

#data:  ac2$sample.size
# W = 0.73257, p-value < 2.2e-16

# p<0.05 indicates non-normality


#=======================================

# 3. BALANCE
# Examine the balance of categorical variables

table(ac2$salmonids)
#non-salmonids  salmonids 
#457            98 
# NOT BALANCED (but variance homogenous)

table(ac2$continent)
#Australia        Europe North America         Other 
#  62           126           283            84 
# NOT BALANCED (but variance homogenous)

table(ac2$system)
#freshwater     marine 
#258        297
# BALANCED (but variance homogenous)

table(ac2$migration)
#diadromous non-diadromous 
# 130            425 
# NOT BALANCED (but variance homogenous)

table(ac2$Red.list.category)
#CR  DD  EN  EX  LC  NE  NT  VU 
#0  30   8  45   1 289  56  28  98 
# NOT BALANCED - can't use this (unless levels are collapsed together)

table(ac2$year)
#2019 2020 2021 2022 2023 
#91  114  174   87   89 
# NOT BALANCED (variance is fine)


#=======================================

# 4. CALCULATE NUMBER OF ZEROS

# What is the percentage of zeros for sample size?

sum(ac2$sample.size == 0) * 100 / nrow(ac2)

#0% zeros

#=======================================

# 5. COLLINEARITY
#
# Fig. 2.5: a summary using the ggpairs command from the GGally library
ac2 %>% 
  ggpairs(columns = c("system", "salmonids", "migration", "sqrt.journal.IF"), 
          aes(alpha = 0.8), 
          lower = list(combo = wrap("facethist", binwidth = 2))) + 
  My_theme


#=======================================
# 6. PLOT RELATIONSHIPS

# Tidy up names and ensure factors properly designated
ac2$fSyst <- as.factor(ac2$system)
ac2$fYear <- as.factor(ac2$year)
ac2$fSal  <- as.factor(ac2$salmonids)
ac2$fCont <- as.factor(ac2$continent)
ac2$fMig  <- as.factor(ac2$migration)

ac2$ss  <- ac2$sample.size 
ac2$jif <- ac2$sqrt.journal.IF  

# System
ggplot(ac2, aes(x = jif, y = (ss))) +
  geom_jitter(shape = 16, size = 5, alpha = 0.5, height = 0.5, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, linewidth = 1.5) +
  My_theme +
  xlab("JIF") + ylab("Sample size") +
  facet_grid(.~fSyst)

# Year
ggplot(ac2, aes(x = jif, y = (ss))) +
  geom_jitter(shape = 16, size = 5, alpha = 0.5, height = 0.5, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, linewidth = 1.5) +
  My_theme +
  xlab("JIF") + ylab("Sample size") +
  facet_grid(.~fYear)

# Salmonids
ggplot(ac2, aes(x = jif, y = (ss))) +
  geom_jitter(shape = 16, size = 5, alpha = 0.5, height = 0.5, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, linewidth = 1.5) +
  My_theme +
  xlab("JIF") + ylab("Sample size") +
  facet_grid(.~fSal)

# Continent
ggplot(ac2, aes(x = jif, y = (ss))) +
  geom_jitter(shape = 16, size = 5, alpha = 0.5, height = 0.5, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, linewidth = 1.5) +
  My_theme +
  xlab("JIF") + ylab("Sample size") +
  facet_grid(.~fCont)

# Migration
ggplot(ac2, aes(x = jif, y = (ss))) +
  geom_jitter(shape = 16, size = 5, alpha = 0.5, height = 0.5, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, linewidth = 1.5) +
  My_theme +
  xlab("JIF") + ylab("Sample size") +
  facet_grid(.~fMig)

#=======================================
# 7. INDEPENDENCE

# data on sample size are independent, however, sample size and analysed sample 
# are dependent as analysed sample derive from sample size.

#=======================================
# 
# The data exploration showed:
#   
# 1.	There are outliers in the response variable sample size (drop these)
# 2.	A non-normally distributed and not homogenous response variable (use Poisson)
# 3.  Most of categorical variables are not balanced (families, location, year,
#     migration and IUCN red list category) (but variances are similar)
# 4.	No zeros in the response variable.
# 5.	There could be some collinearity, however, it was difficult to see on a plot. (not sure...)
# 6.	A potential interactions?
# 7.	Independency of the response variable, however, analysed sample size is 
#     dependant on sample size.
#

#=======================================

# Sample size is a count

# A Poisson distribution is an appropriate starting point
# Model needs to generalise among continents and years, include these
# as random terms.
# Also, salmonids are usually migratory - so also include migration as random  

# Fit Poisson GLM
pois1 <- glmmTMB(ss ~ jif + fSyst + fSal + 
                      (1|fYear) + (1|fCont) + (1|fMig),
                      family = poisson(link = "log"),
                      data = ac2)
check_overdispersion(pois1)


#    dispersion ratio =    79.843
# Pearson's Chi-Squared = 43754.230
#                 p-value =   < 0.001

#=======================================
# Why is the model overdispersed?
# A. Outliers?                  ==> We dealt with this...
# B. Missing variables?         ==> Possible interactions?
# C. Zero inflation?            ==> No
# D. High variance?             ==> Possible
# E. Correlation?               ==> Using a GLMM
# F. Non-linear patterns        ==> GAM(M)?
# G. Wrong link function        ==> Change it? 
#=======================================

pois2  <- glmmTMB(ss ~ jif * fSyst + jif * fSal + 
                      (1|fYear) + (1|fCont) + (1|fMig),
                      family = poisson(link = "log"),
                      data = ac2)

check_overdispersion(pois2)
# Slight improvement

# Try negative binomial
nb1 <- glmmTMB(ss ~ jif * fSyst + jif * fSal + 
                 (1|fYear) + (1|fCont) + (1|fMig),
                    family = nbinom1(link = "log"),
                    data = ac2)

check_overdispersion(nb1)
# dispersion ratio =   1.260
# Pearson's Chi-Squared = 686.569
#                 p-value = < 0.001


# Generalised Poisson
gp1 <- glmmTMB(ss ~ jif * fSyst + jif * fSal + 
                   (1|fYear) + (1|fCont) + (1|fMig),
                   family = genpois(link = "log"),
                   data = ac2)
check_overdispersion(gp1)
# dispersion ratio =   0.793
# Pearson's Chi-Squared = 432.266
#                 p-value = 1


round(AIC(pois1,pois2,nb1,gp1),0)

#       df   AIC
# pois1  7 37183
# pois2  9 36299
# nb1   10  5734
# gp1   10  5697 <- najlepszy

# Can we improve model fit?
drop1(gp1)

#           Df    AIC
# <none>       5697.1
# jif:fSyst  1 5703.1
# jif:fSal   1 5696.1 <- drop this interaction

gp2 <- glmmTMB(ss ~ jif * fSyst + fSal + 
                 (1|fYear) + (1|fCont) + (1|fMig),
               family = genpois(link = "log"),
               data = ac2)

drop1(gp2)

#           Df    AIC
# <none>       5696.1
# fSal       1 5694.9 <- drop salmonids
# jif:fSyst  1 5702.8

gp3 <- glmmTMB(ss ~ jif * fSyst +
                   (1|fYear) + (1|fCont) + (1|fMig),
                   family = genpois(link = "log"),
                   data = ac2)

drop1(gp3)
# Cannot improve the model further

# 1. Plot residuals against fitted values
# Simulate data using model parameters
Res1 <- simulateResiduals(fittedModel = gp3, plot = F)

# Plot standardised residuals
par(mfrow = c(1,1), mar = c(5,5,4,4))
plotResiduals(Res1)


# 2. Plot model residuals against covariates in the model
par(mfrow = c(1,2), mar = c(5,5,4,4))
plotResiduals(Res1, form = ac2$jif)
plotResiduals(Res1, form = ac2$fSyst)


# Check normality of residuals, dispersion and presence of outliers
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Res1, testOutliers = TRUE, testDispersion = TRUE)


# Model summary
tab_model(gp3,
          show.zeroinf = TRUE,
          dv.labels = c("Generalised Poisson GLM (bison)"),
          string.pred = "Coeffcient",
          pred.labels = c("Intercept(fw)", "Impact Factor",
                          "System(marine)", "IF:System(marine)"),
          string.ci = "Conf. Int (95%)",
          collapse.ci = FALSE,
          auto.label = FALSE,
          string.p = "P-value",
          p.style = c("numeric"),
          emph.p = FALSE,
          transform = NULL)
# Most variance explained by random part of the model

# Visualise the results
plot_model(gp1, type = "eff", 
           terms = c("jif", "fSyst"),
           colors = c("blue2","firebrick2"),
           show.data = TRUE,
           title = "", 
           axis.title = c("sqrtJIF",
                          "Sample size"),
           show.legend = T) + My_theme +
scale_y_continuous(limits = c(0, 300))
# Interpretation?

# For marine species, higher impact papers have larger sample sizes,
# but the opposite is true for freshwater papers; high impact studies have 
# lower sample sizes.

# Are the axes the wrong way here? i.e. should sample size predict IF?
