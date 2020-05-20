##
# GIF Converter
# ffmpeg -i output.mp4 -vf "fps=10,scale=1024:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 output.gif
# VIDEO FROM IMAGES
# ffmpeg -i SEIQRDP_ARG%3d.png -vcodec libx264  -vf "setpts=3*PTS" SEIQRDP_ARG.mp4 
#



library(stringr)
LT <- 5    #1-14 days, latent time in days, incubation period, gamma^(-1)   
QT <- 5  #Quarantine time in days, recovery time, infectious period, delta^(-1)
alpha_guess <- 1.0
beta_guess <- 1.0
gamma_guess <- 1 / LT
delta_guess <- 1 / QT
kappa_guess <- 0.01 #dead rate
lambda_guess <- 0.08 # recover rate
region = "Argentina"
population =  46e6

data<-get_jhu_data(region=region,population = population )

for (i in 0:30){

  start_date<-(today()-i-7) 
  fitted_date<-(today()-i) 
  start_date <- start_date %>% format('%m-%d-%Y')
  fitted_date <- fitted_date %>% format('%m-%d-%Y')

  print(start_date)
  print(fitted_date)

  data_fit<-create_fit_data(data,fitted_date = fitted_date,start_date = start_date)
  parameters<-c(alpha_guess,beta_guess,gamma_guess,delta_guess,kappa_guess,lambda_guess)



  params<-SEIQRDP_fit(parameters,data_fit)
  forecast <-90
  forecast_data <- SEIQRDP_predict(forecast,params,data_fit,data)
  plot_arg<-SEIQRDP_plot(forecast_data,region)


  ggsave(plot_arg,filename = paste0("./SEIQRDP_7_ARG",str_pad(30-i, 3, pad = "0"),".png"),height = 4,width = 9) 
}


for (i in 7:14){
  
  start_date<-(today()-i) 
  fitted_date<-(today()) 
  start_date <- start_date %>% format('%m-%d-%Y')
  fitted_date <- fitted_date %>% format('%m-%d-%Y')
  
  print(start_date)
  print(fitted_date)
  
  data_fit<-create_fit_data(data,fitted_date = fitted_date,start_date = start_date)
  parameters<-c(alpha_guess,beta_guess,gamma_guess,delta_guess,kappa_guess,lambda_guess)
  
  
  
  params<-SEIQRDP_fit(parameters,data_fit)
  forecast <-90
  forecast_data <- SEIQRDP_predict(forecast,params,data_fit,data)
  plot_arg<-SEIQRDP_plot(forecast_data,region)
  
  
  ggsave(plot_arg,filename = paste0("./SEIQRDP_fit_ARG",str_pad(i-7, 3, pad = "0"),".png"),height = 4,width = 9) 
}
