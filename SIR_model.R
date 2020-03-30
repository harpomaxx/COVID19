# Code adapted from 
# https://timchurches.github.io/blog/posts/2020-02-18-analysing-covid-19-2019-ncov-outbreak-data-with-r-part-1/#estimating-changes-in-the-effective-reproduction-number

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
require(deSolve)

population<-list("CABA"=2890000,"Tucumán"=1592878,"Santa Fe"=3369000,"Mendoza"=1741610)
province="Santa Fe"
covid19_arg<-read_delim("COVID19_ARG.tsv",delim ='\t')
covid19_arg<-covid19_arg %>% mutate(date=mdy(date))
cumulative_incidence<-covid19_arg[,c(12:32)] %>% replace(is.na(.), 0) %>% 
  tibble::add_column(date=covid19_arg$date) %>%group_by(date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "new_cases") %>% 
  group_by(Province) %>% mutate(cumulative_cases=cumsum(new_cases)) %>% rename(Date=date)


# put the daily cumulative incidence numbers 
sir_start_date <- (cumulative_incidence %>% filter(Province == province & cumulative_cases!=0) %>% 
  select(Date)   %>% top_n(-1))[2]
sir_fitted_date <- sir_start_date + 10

sir_start_date <- sir_start_date[[1]]
sir_fitted_date <- sir_fitted_date[[1]]

Infected <- cumulative_incidence %>% filter(Province == province, 
                                            Date >= ymd(sir_start_date), 
                                            Date <= ymd(sir_fitted_date)) %>% 
  pull(cumulative_cases)

# Create an incrementing Day vector the same length as our
# cases vector
Day <- 1:(length(Infected))
N=population[[province]]
#N=2890000 #Cdad Tucumán
#N=600000 # PCIA_BS_AS
#N=1000000 #MZACity
#N=1741610 #MZA
#N=1592878 #TUCUMAN
#N=3369000 #SANTAFE

# now specify initial values for S, I and R
init <- c(S = N - Infected[1], I = Infected[1], R = 0)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S/N
    dI <- beta * I * S/N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


# define a function to calculate the residual sum of squares
# (RSS), passing in parameters beta and gamma that are to be
# optimised for the best fit to the incidence data
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[, 3]
  sum((Infected - fit)^2)
}

# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0
Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

# time in days for predictions
t <- 1:as.integer(today() - ymd(sir_start_date))
# get the fitted values from our SIR model
fitted_cumulative_incidence <- data.frame(ode(y = init, times = t, 
                                              func = SIR, parms = Opt_par))
# add a Date column and join the observed incidence data
fitted_cumulative_incidence <- fitted_cumulative_incidence %>% 
  mutate(Date = ymd(sir_start_date) + days(t - 1), Province = province) %>% 
  left_join(cumulative_incidence %>% ungroup()  %>% filter(Province == 
                                                            province) %>% 
select(Date, cumulative_cases))

# plot the data
fitted_cumulative_incidence %>% filter(Date >= ymd(sir_start_date)) %>% 
  ggplot(aes(x = Date)) + 
  geom_line(aes(y = I), colour = "red") + 
  annotate(geom="text",x=ymd("2020-03-20")-0.3,y=10,label="Quarantine start",angle=90,size=2.5)+
  geom_vline(xintercept = ymd("2020-03- 20"), 
             color = "blue", size=0.5)+
  annotate(geom="text",x=ymd(sir_fitted_date)-0.3,y=10,label="Fitted Model",angle=90,size=2.5)+
  geom_vline(xintercept = ymd(sir_fitted_date), 
             color = "green", size=0.5)+
  geom_point(aes(y = cumulative_cases), colour = "orange") + 
  scale_y_continuous(trans='log10')+
  scale_x_date(date_breaks = "1 day",date_labels = "%d %b")+
  labs(y = "Cumulative incidence [log10]", title = paste("COVID-19 fitted vs observed cumulative incidence,",province), 
       subtitle = "(red=fitted incidence from SIR model, orange=observed incidence)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename=paste("images/sir",province,".png",sep=""),height = 3, width = 6)
