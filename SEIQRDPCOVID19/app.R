#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(shinythemes)
library(plotly)
source("SEIQRDP.R")


# Global config
LT <- 5    #1-14 days, latent time in days, incubation period, gamma^(-1)   
QT <- 5  #Quarantine time in days, recovery time, infectious period, delta^(-1)
alpha_guess <- 1.0
beta_guess <- 1.0
gamma_guess <- 1 / LT
delta_guess <- 1 / QT
kappa_guess <- 0.01 #dead rate
lambda_guess <- 0.08 # recover rate

countries_pop <- list(
                  "Austria"=9e6,
                  "Argentina"=46e6,
                  "Bolivia"=11.3e6, 
                  "Brazil"=219e6,
                  "Chile"=18e6,
                  "Czechia"=10e6,
                  "Colombia"=50e6,
                  "Ecuador"=17e6,
                  "Germany"=83e6,
                  "Italy"=61e6,
                  "Korea, South"=60e6,
                  "Paraguay"=7e6,
                  "Peru"=31.9e6,
                  "Russia"=144e6,
                  "Spain"=46e6,
                  "Uruguay"=3.4e6,
                  "US"=328e6,
                  "Venezuela"=28e6
                  
                  
                  )

valueBox <- function(value, subtitle, icon, color) {
  div(class = "col-lg-6 col-md-6",
      div(class = "panel panel-primary",
          div(class = "panel-heading", style = paste0("background-color:", color),
              div(class = "row",
                  div(class = ("col-xs-12 text-left"),
                      div(style = ("font-size: 22px; font-weight: bold; "),
                          textOutput(value)
                      ),
                      div(style = ("font-size: 10px; "),
                        subtitle)
                  )
              )
          ),
          #div(class = "panel-footer",
          #    div(class = "clearfix")
          #)
      )
  )
}



ui <- fluidPage(
  #theme = shinythemes::shinytheme("yeti"), 
  theme = shinythemes::shinytheme("darkly"),
  title = "Yet another SEIQRDP COVID19 Projection",
  
 # plotOutput('modelPlot',height="480px") %>% withSpinner() ,
 plotOutput('modelPlot') %>% withSpinner() ,
 
 fluidRow(
   column(12,
        div(style="text-align:left;
        font-size:9px;",
        HTML("Data Source from Johns Hopkins CSSE. SEIQRDP. Model from (Peng et al. 2020).<BR>
        <a href=http://labsin.org> http://labsin.org</a>. (<a href=https://twitter.com/harpolabs> @harpolabs </a> & 
             <a href=https://twitter.com/RGonzalez_PhD> @RGonzalez_PhD </a>)<br><br><br><br>")
        )
   )
  ),
 
  fluidRow(
    column(3,
           
           
          # fluidRow(
           valueBox(value = "peak",
                    subtitle = "for the reaching the peak",
                    icon = "tachometer",
                    color = "#1979a9"), 
         #  ),
         #  fluidRow(
           valueBox(value = "double",
                    subtitle = "for doubling the number of cases",
                    icon = "tachometer",
                    color = "orange")
       
          # )
       
          ),
    
    column(2,
           
           checkboxGroupInput("compartiment", label = "Compartiments : ", 
                              choices = list("Infected" = 1, "Exposed" = 2 ,"Recovered" =3, "Dead" = 4),
                              selected = 3)
           
    ),
    
    
    column(2,
           #div(style="height: 20px;font-size: 15px;",
           sliderInput("Pdays",
                       "Days considered for projecting:",
                       min = 1,
                       max = 365,
                       value = 30) 
           #)
    ),
    column(2, offset = 0,    
           #div(style="height: 20px;font-size: 15px;",
           sliderInput("Fdays",
                     "Days considered for fitting:",
                     min = 1,
                     max = 60, 
                     value = 10) 
           #)
           
    ),
    column(2,
           #div(style="height: 20px;font-size: 15px;",
           selectInput("region", label = "Select Country", 
                       choices = list(
                                  "Argentina"="Argentina",
                                  "Austria"="Austria",
                                  "Bolivia"="Bolivia",
                                  "Brazil"="Brazil",
                                  "Chile"="Chile",
                                  "Czechia"="Czechia",
                                  "Colombia"="Colombia",
                                  "Ecuador"="Ecuador",
                                  "Germany"="Germany",
                                  "Italy"="Italy",
                                  "Korea, South"="Korea, South",
                                  "Paraguay"="Paraguay",
                                  "Peru"="Peru",
                                  "Russia"="Russia",
                                  "Spain"="Spain",
                                  "Uruguay"="Uruguay",
                                  "US"="US",
                                  "Venezuela"="Venezuela"
                                      
                       ), 
                       selected = "Argentina",selectize=FALSE)
           ),
    column(2, 
           dateInput("in_date", 
                     label = "Current Date",
                                        value = today()-15,
                                         min   = today()-60,
                                         max   = today()-2
                     
                     )
           
           )
      
  )
)


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
#  observe({
#    data<- data<-get_jhu_data(region=input$region,population = countries_pop[[input$region]] )
  
#    #mindate<-ymd(((model$data$cumulative_incidence %>% filter (infected_cases > 50) %>% select (Date) %>% filter(row_number() ==1)))[[1]])
#    date <- today()
#    mindate <- today()
#    #print(mindate)
#    updateDateInput(session, "in_date",
#                    label = "Select Date",
#                    value = date,
#                    min   = today()-30,
#                    max   = today()
#    )
#  })

  
  
  build_model <- reactive ({
    print(paste("Fitted",input$Fdays,"days"))
    start_date<-(input$in_date-(input$Fdays)) 
    #fitted_date<-(start_date+input$Fdays) 
    fitted_date<-(input$in_date) 
    
    start_date <- start_date %>% format('%m-%d-%Y')
    fitted_date <- fitted_date %>% format('%m-%d-%Y')
  
    data<-get_jhu_data(region=input$region,population = countries_pop[[input$region]] )
    
    data_fit<-create_fit_data(data,fitted_date = fitted_date,start_date = start_date)
    parameters<-c(alpha_guess,beta_guess,gamma_guess,delta_guess,kappa_guess,lambda_guess)
    params<-SEIQRDP_fit(parameters,data_fit)
    return (list("data"=data,
                "data_fit"=data_fit,
                "params"=params))
    
  }) 
  
   output$modelPlot <- renderPlot({

     model<-build_model()
     
     #print(model)
     forecast <-input$Pdays+1
     forecast_data <- SEIQRDP_predict(forecast,model$params,model$data_fit,model$data)
     plot_arg<-SEIQRDP_plot(forecast_data,paste(input$region," [Pop. ", ms(countries_pop[[input$region]]) ,"]",sep=""), input$compartiment)
     plot_arg
       
    }) 
   
   output$peak <- renderText({

     model<-build_model()
     forecast <-input$Pdays+1
     forecast_data <- SEIQRDP_predict(forecast,model$params,model$data_fit,model$data)
     double_time <- calculate_double_time(forecast_data)
     peak <- calculate_peak(forecast_data)
     print(paste0(ymd(peak$date) - today(), "d"))
     }) 
   
   output$double <- renderText({
     model<-build_model()
     forecast <-input$Pdays+1
     data <- SEIQRDP_predict(forecast,model$params,model$data_fit,model$data)
     double_time <- calculate_double_time(data)
     print(paste0(double_time,"d"))
     
   }) 
  
}

# Run the application 
shinyApp(ui = ui, server = server)

