#--------------------------------------------------------------------------------------------------------------------------
# Okuns law and potential output - recreation of Lancaster and Tulips 2015 RDP
#--------------------------------------------------------------------------------------------------------------------------

#setwd("C:/Users/aelde/OneDrive/Documents/GitHub/Okuns law")

library(dlm)
library(dplyr)
library(tidyverse)
library(readxl)
library(lubridate)
library(stargazer)
source("~/GitHub/TST themes/Chart themes.R")
source("~/GitHub/Packages/tst.package/R/tst.macrodata.R")

#--------------------------------------------------------------------------------------------------------------------------
# import data 
#--------------------------------------------------------------------------------------------------------------------------

# Data series used:
#   Unemployment rate
#   real unit labour costs (authors estimates - COE+TLS/(GDP*Deflator)) - data is spliced with rulc from NA growth
#   Real gdp


tldata <- read_excel("C:/Users/aelde/OneDrive/Documents/GitHub/Okuns law/TLdata2018.xlsx", sheet = "Sheet1", range = "A1:D238")


#--------------------------------------------------------------------------------------------------------------------------
# data preperation
#--------------------------------------------------------------------------------------------------------------------------


# Overwrite and update TL data
colnames(tldata)[1] <- "Date"
tldata$Date <- as.Date(tldata$Date)
lastdate <-  last(m$Date)

#load macrodata
m <- tst.macrodata()

updatetl <- right_join(tldata, m[,c("Date","rate_UNE","GDPE_r","ULC_r")]) %>% 
  mutate(ur = if_else(is.na(ur), rate_UNE, ur),
         gdp = GDPE_r,
         rulc =if_else(is.na(rulc), lag(rulc)*ULC_r/lag(ULC_r), rulc) 
         )


dur <- updatetl$ur-lag(updatetl$ur)                                          # first difference of the unemployment rate
d2lgdp <- 400*(log(updatetl$gdp)-log(lag(updatetl$gdp,2)))/2                 # Two-quarter GDP growth, annualised log changes
d2lrulc <- 400*(log(updatetl$rulc)-log(lag(updatetl$rulc,2)))/2              # Two-quarter rulc growth, annualised log changes

mod_data <- data.frame(
  dur = dur,
  
  d2lgdp = d2lgdp,
  
  d2lrulc = d2lrulc
)

mod_data <- mod_data[-c(1:2),]

mod_data <- mod_data %>% 
  mutate(dur_lag1 = lag(dur),
         d2lrulc_lag2 =lag(d2lrulc,2)) %>% 
  filter(!is.na(d2lrulc_lag2))


#--------------------------------------------------------------------------------------------------------------------------
# constant coefficients model
#-------------------------------------------------------------------------------------------------------------------------

  
const_coef <- nls(formula = dur~ b1*dur_lag1 + b2*(d2lgdp-b0) + b3*d2lrulc_lag2 ,
    start = list(b0 =0.1, b1=0.1, b2=0.1, b3=0.1),
    data = mod_data) 

#--------------------------------------------------------------------------------------------------------------------------
# tvp model full model - dur = alpha*dur(-1)+ beta(dgdp-potential) + gamma*wages
#--------------------------------------------------------------------------------------------------------------------------

beta.start <- as.vector(coef(const_coef)[3])

# Construct DLM

OkunsDLMfm <- dlm(
  
  
  FF = matrix(c(1,1,1,1),ncol = 4, byrow = TRUE),
  
  V = matrix(1),
  
  GG = matrix(c(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1), ncol = 4, byrow = TRUE),
  
  W =  matrix(c(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1), ncol = 4, byrow = TRUE),
  
  JFF = matrix(c(1,2,3,0),ncol = 4, byrow = TRUE),
  
  X = cbind(mod_data$dur_lag1,mod_data$d2lgdp, mod_data$d2lrulc_lag2),
  
  m0 = c(0,0,0,0),
  
  C0 = matrix(c(1e+07,0,0,0,
                0,1e+07,0,0,
                0,0,1e+07,0,
                0,0,0,1e+07), ncol = 4, byrow = TRUE)
  
)


buildOkunsFM <- function(p){
  
  V(OkunsDLMfm)  <- exp(p[1])
  
  GG(OkunsDLMfm)[1,1]  <- 1
  
  GG(OkunsDLMfm)[2,2]  <- 1
  
  GG(OkunsDLMfm)[3,3]  <- 1 
  
  GG(OkunsDLMfm)[4,4]  <- 1
  
  W(OkunsDLMfm)[1,1] <- exp(p[2])
  
  W(OkunsDLMfm)[2,2] <- 0
  
  W(OkunsDLMfm)[3,3] <- 0
  
  W(OkunsDLMfm)[4,4] <- exp(p[3])
  
  m0(OkunsDLMfm) <- c(0,0,0,0)
  
  C0(OkunsDLMfm)[1,1] <- 1
  
  C0(OkunsDLMfm)[4,4] <- 5
  
  
  return(OkunsDLMfm)
  
}



okuns.estfm <-  dlmMLE(y = mod_data$dur, parm = c(-1.4,-6,-5), build = buildOkunsFM)


OkunsDLM1fm <- buildOkunsFM(okuns.estfm$par)


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

filtered.fm <- dlmFilter(y = mod_data$dur, mod = OkunsDLM1fm)

smoothed <- dlmSmooth(y = mod_data$dur, mod = OkunsDLM1fm)

variances <-  dlmSvd2var(filtered.fm$U.C,filtered.fm$D.C) %>%
  lapply(sqrt) %>%
  sapply(diag) %>%
  t() %>% 
  data.frame()

variances$Date <- seq(as.Date("1960-06-01"),as.Date(lastdate), by = "quarter")

Upperb <- filtered.fm$m[-1,2]+ variances[-1,2]#*qnorm(0.05,lower = FALSE)
Lowerb <- filtered.fm$m[-1,2]- variances[-1,2]#*qnorm(0.05,lower = FALSE)


Upperbp <- filtered.fm$m[-1,4]+variances[-1,4]#*qnorm(0.05,lower = FALSE)
Lowerbp <- filtered.fm$m[-1,4]-variances[-1,4]#*qnorm(0.05,lower = FALSE)



chart1 <- data.frame(GDPgrowth = mod_data$d2lgdp,
           Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
           upper = Upperbp/(abs(Upperb)),
           lower = Lowerbp/(abs(Lowerb)),
           Date = seq(as.Date("1960-09-01"),as.Date(lastdate), by = "quarter")) %>% 
  gather(Var, Val, -Date) %>% 
  filter(Var != "GDPgrowth" & Date >= "1979-12-01") %>%
  spread(Var, Val) %>% 
  ggplot()+
  geom_line(aes(Date,Potential), colour =tst_colors[2]  , size = 1)+
  geom_ribbon(aes(x = Date, ymin = lower, ymax = upper) , alpha = 0.2)+
  tst_theme()+
  theme(legend.position = "none")+
  ylim(0.5,8.5)+
  ylab("")+
  xlab("")+
  ggtitle("Potential output growth",subtitle = "annualised % change")+
  annotate("segment", x=ymd("2010-06-01"), xend =ymd("2012-06-01") , y= 6, yend= 4.5, arrow = arrow(), colour = "dark grey" )+
  annotate("text", x=ymd("2010-06-01") , y= 6.5, label = "+/- One S.D", colour = "dark grey")+
  annotate("segment", x=ymd("2012-06-01"), xend =ymd("2017-06-01") , y= 1, yend= 2.5, arrow = arrow(), colour = tst_colors[2] )+
  annotate("text", x=ymd("2012-06-01") , y=0.75 , label = "Potential output growth", colour = tst_colors[2])


chart2 <-  data.frame(GDPgrowth = mod_data$d2lgdp,
           Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
           upper = Upperbp/(abs(Upperb)),
           lower = Lowerbp/(abs(Lowerb)),
           Date = seq(as.Date("1960-09-01"),as.Date(lastdate), by = "quarter")) %>% 
  gather(Var, Val, -Date) %>% 
  filter(!Var %in% c("upper","lower") & Date >= "1979-12-01") %>%
  spread(Var, Val) %>% 
  ggplot()+
  geom_line(aes(Date, GDPgrowth), colour = tst_colors[1], size = 1)+
  geom_line(aes(Date,Potential), colour = tst_colors[2] , size = 1)+
  tst_theme()+
  scale_colour_tst()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  ggtitle("Potential and actual output growth",subtitle = "annualised % change")+
  annotate("text", x=ymd("2011-06-01") , y= 5.5, label = "Real GDP growth", colour = tst_colors[1])+
  annotate("text", x=ymd("2012-06-01") , y=0 , label = "Potential output growth", colour = tst_colors[2])



States <- data.frame(Date = seq(as.Date("1960-09-01"),as.Date(lastdate), by = "quarter"), 
                     GDPgrowth = mod_data$d2lgdp,
                     Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
                     alpha = dropFirst(filtered.fm$m[,1]),
                     beta = dropFirst(filtered.fm$m[,2]),
                     gamma = dropFirst(filtered.fm$m[,3])
                     
                     
                     )


#--------------------------------------------------------------------------------------------------------------------------
# recreating estiamtes in table 1
#--------------------------------------------------------------------------------------------------------------------------

fullsample <-tibble(Model = "TVP - full sample",
                        alpha = last(States$alpha),
                        se_alpha = last(variances$X1),
           
                        beta = last(States$beta),
                        se_beta = last(variances$X2),
           
           `Potential output growth` = last(States$Potential),
           gamma = last(States$gamma),
           se_gamma =last(variances$X2),
           
           `Okun's coef` = 4*last(States$beta)/(1-last(States$alpha))
)%>% 
  gather(Variable, `TVP full sample`, -Model) %>% 
  dplyr::select(-Model)


tlsample <- tibble(Model = "TVP - Estimate at 2015Q1",
                      alpha = States[States$Date =="2015-03-01","alpha"],
                      se_alpha = variances[variances$Date == "2015-03-01","X1"],
                      
                       beta = States[States$Date =="2015-03-01","beta"] ,
                      se_beta = variances[variances$Date == "2015-03-01","X2"],
                       `Potential output growth` = States[States$Date =="2015-03-01","Potential"],
                       gamma = States[States$Date =="2015-03-01","gamma"],
                      se_gamma = variances[variances$Date == "2015-03-01","X3"],
                       `Okun's coef` = 4*States[States$Date =="2015-03-01","beta"]/(1-States[States$Date =="2015-03-01","alpha"])
) %>% 
  gather(Variable, `TVP est at 2015Q1`, -Model) %>% 
  dplyr::select(-Model)



tlorig <- tibble(Model = "TL (2015)",
                     alpha = 0.443,
                     se_alpha = 0.251,
                     
                     beta = -0.049 ,
                     se_beta = 0.006,
                     `Potential output growth` = 2.94,
                     gamma = 0.019,
                     se_gamma = 0.004,
                     `Okun's coef` = -0.35
) %>% 
  gather(Variable, `TL (2015)`, -Model) %>% 
  dplyr::select(-Model)

tabledat <- left_join(fullsample,
      tlsample) %>% left_join(tlorig)

#--------------------------------------------------------------------------------------------------------------------------
# Save workspace
#--------------------------------------------------------------------------------------------------------------------------


save.image("~/GitHub/Okuns law/OkunsLaw.RData")

