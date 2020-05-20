#remove (list = objects())
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
require(deSolve)
library(ggplot2)
library(scales)

## Get data from John Hopkins dataset
get_jhu_data<-function(region="Argentina",population=40e7){
  
  # Get dataset from Hopkins University 
  jhu_url_confirmed <-
    "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
  jhu_url_recovered <-
    "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
  jhu_url_dead <-
    "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
  
  # Filter 
  confirmed_long_jhu <-
    read_csv(jhu_url_confirmed) %>% rename(province = "Province/State", country_region = "Country/Region") %>%
    filter(country_region == region) %>%
    select(-province, -Lat,-Long) %>%
    pivot_longer(-country_region, names_to = "Date", values_to = "infected_cases")
  
  recovered_long_jhu <-
    read_csv(jhu_url_recovered) %>% rename(province = "Province/State", country_region = "Country/Region") %>%
    filter(country_region == region) %>%
    select(-province, -Lat,-Long) %>%
    pivot_longer(-country_region, names_to = "Date", values_to = "recovered_cases")
  
  dead_long_jhu <-
    read_csv(jhu_url_dead) %>% rename(province = "Province/State", country_region = "Country/Region") %>%
    filter(country_region == region) %>%
    select(-province, -Lat,-Long) %>%
    pivot_longer(-country_region, names_to = "Date", values_to = "dead_cases")

  cumulative_incidence <-
    confirmed_long_jhu %>% select(Date, infected_cases) 
  
  cumulative_recover <-
    recovered_long_jhu %>% select(Date, recovered_cases)
  
  cumulative_dead <- dead_long_jhu %>% select(Date, dead_cases)
  
  cumulative_incidence<-
   cumulative_incidence %>%
   inner_join(cumulative_recover) %>% 
   inner_join(cumulative_dead) %>% 
   select(Date,infected_cases,recovered_cases,dead_cases)
  
  cumulative_incidence <- cumulative_incidence %>% mutate(Date = mdy(Date))
  return(list("region"=region,"cumulative_incidence"=cumulative_incidence,"population"=population))
}

## subset the data for model fitting
create_fit_data <- function(cases_data,start_date="03-26-2020",fitted_date="04-09-2020"){
  # Extract the vectors used for fitting the model
  fit_data<- cases_data$cumulative_incidence %>% filter(
                                               Date >= mdy(start_date),
                                               Date <= mdy(fitted_date)) 
  
  N=cases_data$population
  Quarantine <- (fit_data$infected_cases - fit_data$recovered_cases - fit_data$dead_cases)
  Infected <- fit_data$infected_cases # Not known, but probably not zero
  Exposed <- fit_data$infected_cases  # Not known, but probably not zero
  Recovered <- fit_data$recovered_cases
  Dead <- fit_data$dead_cases
  init <-
    c(
      S = N - Infected[1] - Exposed[1] - Dead[1] - Quarantine[1] - Recovered[1],
      E = Exposed[1],
      I = Infected[1],
      Q = Quarantine[1],
      R = Recovered[1],
      D = Dead[1]
    )
  
  list("fit_data" = fit_data,
        "init" = init,
        "start_date" = start_date,
        "fitted_date" = fitted_date
        )
}

## Peak
calculate_peak<-function(data){
 
   date<-((data$forecast %>% 
               filter(Q > lag(Q,1))  
             %>% count()))[[1]]+mdy(data$start_date)
  value<-(data$forecast %>% 
               filter(Date == date) %>%  select(Q))[[1]]  
  value<-ifelse(length(value)==0,0,value)
  
  return(list("date"=date,"value"=value))
}


## Doubling time
calculate_double_time<-function(data){

  active_initial<-data$forecast$infected_cases[1] -
    data$forecast$dead_cases[1] -
    data$forecast$recovered_cases[1]

  final_time<-as.integer(mdy(data$fitted_date)-mdy(data$start_date))
  active_final<-data$forecast$infected_cases[final_time] -
    data$forecast$dead_cases[final_time] -
    data$forecast$recovered_cases[final_time]
  
  if (active_final > active_initial){
      doubling_time<-round(log(2)*final_time/
                       (log(data$forecast$infected_cases[final_time]) - 
                          log(data$forecast$infected_cases[1])))
  }else
      doubling_time<-0
    
 return(doubling_time)
}

pretty_print_num<-function(num){
  if ( num >= 1000 ) {
    num <- num/1000
    res <- paste(round(num,digits = 1),"k",sep="")
  } else {
    res <- paste(num) 
  }
  return (res)
}

ks <- function (x) { number_format(accuracy = .1,
                                   scale = 1/1000,
                                   suffix = "k",
                                   big.mark = ",")(x) }

# Plot the results
SEIQRDP_plot<-function(data,province="NA Region"){
  forecast<-nrow(data$forecast) - data$forecast %>% 
    filter(Date == mdy(data$fitted_date)) %>% select(time)
  
  double_time <- calculate_double_time(data)
  peak <- calculate_peak(data)
  ymax <- max(data$forecast$Q)
  plot <-
  ggplot(data$forecast %>% filter(Date >= mdy(data$start_date)),
         aes(x = Date)) +
  
  
  geom_line(aes(y = Q, color="active",linetype = "projected"), alpha=0.3) +
  geom_line(data = data$forecast %>% filter(Date<=mdy(data$fitted_date)), aes(y = Q, colour = "active", linetype="fitted"),size=1) +
  
  geom_line(aes(y = D, colour = "dead", linetype = "projected"),alpha=0.3) +
  geom_line(data = data$forecast %>% filter(Date<=mdy(data$fitted_date)), aes(y = D, colour = "dead",linetype="fitted"),size=1) +
  
  geom_line(aes(y = R, colour = "recover", linetype = "projected"),alpha=0.3) +
  geom_line(data = data$forecast %>% filter(Date<=mdy(data$fitted_date)), aes(y = R, colour = "recover",linetype="fitted"),size=1) +
  
  #geom_line(aes(y = I, colour = "infected",linetype = "dashed")) +
  #geom_line(aes(y = E), colour = "green") +
  #annotate(geom="text",x=mdy(data$fitted_date)-0.3,y=ymax,label="Fitted Model",angle=90,size=2.5,alpha=0.3)+
    
  #geom_vline(xintercept = mdy(data$fitted_date), color = "green", size = 0.5) +
  geom_point(aes(y = ( data$forecast$infected_cases - data$forecast$recovered_cases - data$forecast$dead_cases ),
                 shape="observations"), colour = "orange") +
  geom_point(aes(y = data$forecast$recovered_cases), colour = "blue",shape=1) +
  geom_point(aes(y = data$forecast$dead_cases), colour = "black",shape=1) +
  geom_point(aes(y=peak$value,x=peak$date), colour="orange",size=3,shape=6)+
    
  geom_text(data=
              data$forecast[seq(1,length(data$forecast$Q) ,by = 10),] %>% filter(Date>=mdy(data$fitted_date))
            ,
            aes(x=Date, y=Q+(Q*0.22),label=ks(round(Q))),size=3,color='orange'
  )+
    
    geom_point(data=
                data$forecast[seq(1,length(data$forecast$Q) ,by = 5),]
              ,
              aes(x=Date, y=Q),size=1,color='orange'
    )+  
  #geom_text(data=
  #              data$forecast[seq(1,length(data$forecast$R) ,by = 5),]
  #            ,
  #            aes(x=Date, y=R+(R*0.09),label=pretty_print_num(round(R))),size=3,color='blue'
  #  )+
  
  
  scale_x_date(date_breaks = "3 day", date_labels = "%d %b") +
  scale_y_continuous(labels = ks)+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual("", 
                     #breaks = c("Active", "Recover", "Dead"),
                     values = c( 
                                 "active"="orange", 
                                 "dead"="black", 
                                 "recover"="blue"
                                
                                 #"infected"="yellow",
                                 #"exposed"="green"
                                  )) +
    
    scale_linetype_manual("",values=c("projected"="dashed","fitted"="solid"))+
    scale_shape_manual("",values=c("observations"=1))+


  
 labs(
     y = "Number of Cases", x= "Time in Days",
    title = paste(province, "COVID19 forecast. SEIQRDP model (Peng et al. 2020)"), subtitle=
      paste("Fitted with ", mdy(data$fitted_date)- mdy(data$start_date), " days. Forecasted for ",
            forecast," days from ", mdy(data$fitted_date),".\nPeak expected in ",ymd(peak$date) - today()," days",
            ". Doubling every ",double_time," days.",sep=""), caption = "Source: Johns Hopkins CSSE\n@harpolabs"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="top")+
  guides(
    linetype=guide_legend(keywidth = 2, keyheight = 1),
    colour=guide_legend(keywidth = 2, keyheight = 1))
  #ggsave(filename = paste("images/SEIQRDP_", province, ".png", sep = ""),height = 3, width = 6)
  plot
}

SEIQRDP_model_learned <- function (ctime, state_values, parameters)
{
  kappa <- parameters[[2]]
  lambda <- parameters[[3]]
  parameters<-parameters[[1]]
  with (as.list (c(state_values, parameters)),     # variable names within parameters can be used
        {
        
          N=S+I+E+Q+R+D
          # compute derivatives
          dS = -beta * ((S * I) / N) - alpha * S
          dE =  beta *  ((S * I) / N) - delta * E
          dI = delta * E - gamma * I
          dQ = gamma * I - lambda[ctime] * Q - kappa[ctime] * Q
          dR = lambda[ctime] * Q
          dD = kappa[ctime] * Q
          #dP = alpha * S
          results = c (dS, dE, dI, dQ, dR, dD)
          list (results)
        })
}

SEIQRDP_predict  <- function(time,parameters,fit_data,data){

  t <- 1:as.integer(mdy(fit_data$fitted_date) - mdy(fit_data$start_date)+time)
  
  lambda0 <- parameters['lambda0']
  kappa0 <- parameters['kappa0']
  #lambda <-  lambda0 * (1 - exp(-lambda0 * (t*100)))
  #kappa  <- kappa0  *  exp(-kappa0 * (t*10))
  lambda <-  lambda0 * (1 - exp(-lambda0 * (t*100)))
  kappa  <- kappa0  *  exp(-kappa0 * (t))
  
  forecast <- data.frame(ode(
  y = fit_data$init,
  times =  t ,
  func = SEIQRDP_model_learned,
  parms = list(parameters,kappa,lambda)
))
  
  forecast<-forecast %>% mutate(Date=mdy(fit_data$start_date)+ days(t-1)) %>% 
  left_join(data$cumulative_incidence)
  
  return(list("forecast"=forecast,
              "start_date"=fit_data$start_date,
              "fitted_date"=fit_data$fitted_date
              ))
  
}

# SEIQRDP for fitting lambda and kappa
SEIQRDP_model_fit <- function (ctime, state_values, parameters)
{
  with (as.list (c(state_values, parameters)),     # variable names within parameters can be used
        {
          N=S+I+E+Q+R+D
          # compute derivatives
          dS = -beta * ((S * I) / N) - alpha * S
          dE =  beta *  ((S * I) / N) - delta * E
          dI = delta * E - gamma * I
          dQ = gamma * I - lambda0 * Q - kappa0 * Q
          dR = lambda0 * Q
          dD = kappa0 * Q
          results = c (dS, dE, dI, dQ, dR, dD)
          list (results)
        })
}

# define a function to calculate the residual sum of squares
# (RSS), passing in parameters beta and gamma that are to be
# optimised for the best fit to the incidence data
RSS <- function(parameters,init,fit_data) {
  names(parameters) <-
    c("alpha", "beta", "gamma", "delta", "kappa0", "lambda0")
  Day=1:length(fit_data$infected_cases)
  
  out <-
    ode(
      y = init,
      times = Day,
      func = SEIQRDP_model_fit,
      parms = parameters
    )
  #print(out)
  qfit <- out[, 5]
  rfit <- out[, 6]
  dfit <- out[, 7]
  Quarantine <- fit_data$infected_cases -fit_data$recovered_cases - fit_data$dead_cases
  Recovered <- fit_data$recovered_cases
  Dead <- fit_data$dead_cases

  q =  sqrt(mean(sum((Quarantine - qfit) ** 2)))
  r =  sqrt(mean(sum((Recovered - rfit) ** 2)))
  d =  sqrt(mean(sum((Dead - dfit) ** 2)))
  return (0.33333 * d + 0.33333 * q + 0.33333 * r)
  
}

SEIQRDP_fit <- function (parameters_guess,data){

 par<-as.list (c(parameters_guess,data))
 

 names(par)[1:6]<-c("alpha_guess","beta_guess","gamma_guess","delta_guess","kappa_guess","lambda_guess")
 #print(par)
  Opt <-
    optim(
      c(
        par$alpha_guess,
        par$beta_guess,
        par$gamma_guess,
        par$delta_guess,
        par$kappa_guess,
        par$lambda_guess
      ),
      RSS,
      NULL,
      par$init,
      par$fit_data,
      
      method = "L-BFGS-B",
      upper = c(3, 3, 3, 3, 0.01, 0.1),
      lower = c(0.0, 0.0,0, 0, 0, 0.00),
      hessian = T
      
    )
  #print(Opt$message)
  
  Opt_par <-
    setNames(Opt$par,
             c("alpha", "beta", "gamma", "delta", "kappa0", "lambda0"))
  

  return(Opt_par)
  
}


# MAIN

#guess.QT = 0.5; % 5 weeks, quarantine time in days, recovery time, infectious period, delta^(-1)
#guess.LT = 5; % 11, 1-14 days, latent time in days, incubation period, gamma^(-1)   
#alpha represent the protection rate, 
#beta infection rate, 
#gamma average latent time, 
#delta average quarantine time, 
#lambda cure rate, #
#kappa mortality
#guess.lambda = [0.1, 0.05]; % recovery rate
#guess.kappa  = [0.1, 0.05]; % death rate

#region="Brazil"
#population=21e7
#region ="Italy"
#population=61e6

#region = "Spain"
#population = 46e6
#region = "Argentina"
#population =  40e6

#data<-get_jhu_data(region=region,population = population )
#data_fit<-create_fit_data(data,fitted_date = "04-02-2020",start_date = "03-18-2020")

#LT <- 5    #1-14 days, latent time in days, incubation period, gamma^(-1)   
#QT <- 5  #Quarantine time in days, recovery time, infectious period, delta^(-1)
#alpha_guess <- 1.0
#beta_guess <- 1.0
#gamma_guess <- 1 / LT
#delta_guess <- 1 / QT
#kappa_guess <- 0.01 #dead rate
#lambda_guess <- 0.08 # recover rate
#parameters<-c(alpha_guess,beta_guess,gamma_guess,delta_guess,kappa_guess,lambda_guess)
#params<-SEIQRDP_fit(parameters,data_fit)
#forecast <-30
#forecast_data <- SEIQRDP_predict(forecast,params,data_fit,data)
#plot<-SEIQRDP_plot(forecast_data,region)
#print(params)
#plot
