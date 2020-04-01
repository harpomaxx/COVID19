# Code adapted from 
# https://timchurches.github.io/blog/posts/2020-02-18-analysing-covid-19-2019-ncov-outbreak-data-with-r-part-1/#estimating-changes-in-the-effective-reproduction-number

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
require(deSolve)
library(ggplot2)


jhu_url_confirmed <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
jhu_url_recovered <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
jhu_url_dead <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"


confirmed_long_jhu <- read_csv(jhu_url_confirmed) %>% rename(province = "Province/State", country_region = "Country/Region") %>%  
  filter(country_region=="Argentina") %>% 
  select(-province,-Lat, -Long) %>%
  pivot_longer(-country_region, names_to = "Date", values_to = "infected_cases") 
               
recovered_long_jhu <- read_csv(jhu_url_recovered) %>% rename(province = "Province/State", country_region = "Country/Region") %>%  
  filter(country_region=="Argentina") %>% 
  select(-province,-Lat, -Long) %>%
  pivot_longer(-country_region, names_to = "Date", values_to = "recovered_cases") 

dead_long_jhu <- read_csv(jhu_url_dead) %>% rename(province = "Province/State", country_region = "Country/Region") %>%  
  filter(country_region=="Argentina") %>% 
  select(-province,-Lat, -Long) %>%
  pivot_longer(-country_region, names_to = "Date", values_to = "dead_cases") 


cumulative_incidence<-confirmed_long_jhu %>% select(Date,infected_cases) %>% rename(cumulative_cases=infected_cases)
cumulative_recover<-recovered_long_jhu %>% select(Date,recovered_cases) 

#cumulative_incidence<-inner_join(confirmed_long_jhu,recovered_long_jhu)%>% inner_join(dead_long_jhu)  %>%mutate(cumulative_cases=infected_cases - recovered_cases - dead_cases) %>% select(Date,cumulative_cases)

province="Argentina"
cumulative_incidence$Province=province
cumulative_incidence<-cumulative_incidence %>% mutate(Date=mdy(Date))

cumulative_recover$Province=province
cumulative_recover<-cumulative_recover %>% mutate(Date=mdy(Date))


sir_start_date <- today()-15
sir_fitted_date <- sir_start_date + 10
Infected <- cumulative_incidence %>%filter(Province == province, 
                                            Date >= ymd(sir_start_date), 
                                            Date <= ymd(sir_fitted_date)) %>% 
  pull(cumulative_cases)

Recovered <- cumulative_recover %>%filter(Province == province, 
                                           Date >= ymd(sir_start_date), 
                                           Date <= ymd(sir_fitted_date)) %>% 
  pull(recovered_cases)


Infected <- Infected - Recovered

# Create an incrementing Day vector the same length as our
# cases vector
Day <- 1:(length(Infected))
N=40480000


# now specify initial values for S, I and R
init <- c(S = N - Infected[1], I = Infected[1], R = Recovered[1])

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S/N
    dI <- beta * I * S/N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

SEIR <- function (current_timepoint, state_values, parameters)
{
  #print(state_values)
  # create state variables (local variables)
  #S = state_values [1]        # susceptibles
  #E = state_values [2]        # exposed
  #I = state_values [3]        # infectious
  #R = state_values [4]        # recovered
  
  with ( 
    as.list (c(state_values, parameters)),     # variable names within parameters can be used 
    {
      # compute derivatives
      dS = (-beta * S * I)/N
      dE = (beta * S * I)/N - (delta * E)
      dI = (delta * E) - (gamma * I)
      dR = (gamma * I)
      
      # combine results
      results = c (dS, dE, dI, dR)
      list (results)
    }
  )
}


# define a function to calculate the residual sum of squares
# (RSS), passing in parameters beta and gamma that are to be
# optimised for the best fit to the incidence data
RSS <- function(parameters) {

  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  
  fit1 <- out[, 3]
  fit2 <-out[,4]
  alpha = 0.1
  l2 = sqrt(mean((Recovered-fit2)**2))
  l1= sqrt(mean(sum((Infected - fit1)^2)))
  
  return (alpha * l1 + (1 - alpha) * l2)
  
 }




#beta is the average contact rate in the population
#gamma is the inverse of the mean infectious period (1/t_infectious)
#delta is the inverse of the incubation period (1/t_incubation)


  Opt <- optim(c(0.5,0.5), RSS, method = "L-BFGS-B", lower = c(0.00001,0.00001), 
               upper = c(4, 1))
  
 Opt_par <- setNames(Opt$par, c("beta", "gamma"))
  
  #Opt_par[1]=1.235
  #Opt_par[2]=0.8
  #Opt_par[3]=0.7
  #print(Opt_par)
  # time in days for predictions
  t <- 1:as.integer(today() - ymd(sir_start_date))
  t <- 1:as.integer(ymd(sir_fitted_date) - ymd(sir_start_date)+10)
  # get the fitted values from our SIR model
  fitted_cumulative_incidence <- data.frame(ode(y = init, times = t, 
                                              func = SIR, parms = Opt_par))
  # add a Date column and join the observed incidence data
  fitted_cumulative_incidence <- fitted_cumulative_incidence %>% 
  mutate(Date = ymd(sir_start_date) + days(t - 1), Province = province) %>% 
  left_join(cumulative_incidence %>% ungroup()  %>% filter(Province == 
                                                             province) %>% 
              select(Date, cumulative_cases))
  
  
  fitted_cumulative_incidence <-fitted_cumulative_incidence %>% 
    left_join(recovered_long_jhu %>% 
                select(Date,recovered_cases) %>% 
                mutate(Date=mdy(Date)),by="Date")
  # plot the data
  fitted_cumulative_incidence %>% filter(Date >= ymd(sir_start_date)) %>% 
  ggplot(aes(x = Date)) + 
  geom_line(aes(y = I), colour = "orange") + 
  #geom_line(aes(y = E), colour = "blue") + 
  geom_line(aes(y = R), colour = "black") + 
    
  #annotate(geom="text",x=ymd("2020-03-20")-0.3,y=10,label="Quarantine start",angle=90,size=2.5)+
  #geom_vline(xintercept = ymd("2020-03- 20"), 
  #           color = "blue", size=0.5)+
  annotate(geom="text",x=ymd(sir_fitted_date)-0.3,y=10,label="Fitted Model",angle=90,size=2.5)+
  geom_vline(xintercept = ymd(sir_fitted_date), 
             color = "green", size=0.5)+
  geom_point(aes(y = cumulative_cases-recovered_cases), colour = "orange") + 
  geom_point(aes(y = recovered_cases), colour = "black") + 
  #scale_y_continuous(trans='log10')+
  scale_x_date(date_breaks = "1 day",date_labels = "%d %b")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  labs(y = "Cumulative incidence", title = paste("COVID-19 fitted vs observed cumulative incidence,",province), 
       subtitle = "(red=fitted incidence from SIR model, orange=observed incidence, black=observed recovered)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename=paste("images/sir_",province,".png",sep=""),height = 3, width = 6)

  



