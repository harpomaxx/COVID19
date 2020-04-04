COVID19 ARGENTINA Dataset
================

The dataset
===========

The dataset was produced using *only* the official information provided by the Argentine Ministry of Health. [*available here*](https://www.argentina.gob.ar/coronavirus/informe-diario)

The dataset contains daily information from mars 5th, when the first case was reported. A total of 30 variables are included. All the variables contain daily information with the exceptin of the those starting with `total_`, which refers to cumulative values.

``` r
names(covid19_arg)
```

    ##  [1] "date"                  "total_cases"          
    ##  [3] "total_deaths"          "total_recovered"      
    ##  [5] "total_tests_negatives" "total_tests"          
    ##  [7] "total_imported%"       "total_local_contact%" 
    ##  [9] "total_under_analysis%" "new_cases"            
    ## [11] "new_deaths"            "local_contact"        
    ## [13] "imported"              "under_analysis"       
    ## [15] "CABA"                  "Buenos Aires"         
    ## [17] "Chaco"                 "Córdoba"              
    ## [19] "Corrientes"            "Entre Ríos"           
    ## [21] "Jujuy"                 "La Pampa"             
    ## [23] "Mendoza"               "Misiones"             
    ## [25] "Neuquén"               "Río Negro"            
    ## [27] "Salta"                 "San Luis"             
    ## [29] "Santa Cruz"            "Santa Fe"             
    ## [31] "Santiago del Estero"   "Tierra del Fuego"     
    ## [33] "San Juan"              "La Rioja"             
    ## [35] "Tucumán"               "X36"

Example
=======

Just load the dataset first
---------------------------

``` r
covid19_arg<-read_delim("COVID19_ARG.tsv",delim ='\t')
```

Simple plots using `tidyverse`
==============================

Total and daily number of confirmed infections
----------------------------------------------

``` r
covid19_arg<-covid19_arg %>% mutate(date=mdy(date))
covid19_arg %>% ggplot()+
  geom_col(aes(x=date,y=total_cases),fill='orange')+
  geom_line(aes(x=date,y=new_cases),fill='green')+
  geom_point(aes(x=date,y=new_cases),fill='green',size=2)+
    ylab("confirmed infections")+
  theme_bw()+
    scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave("images/total.png",height = 2,width = 6)
```

![](README_files/figure-markdown_github/pressure-1.png)

![](./images/total.png)

Total and daily number of confirmed deaths
------------------------------------------

``` r
covid19_arg %>% ggplot()+
  geom_col(aes(x=date,y=total_deaths),fill='red')+
  geom_line(aes(x=date,y=new_deaths),fill='green')+
  geom_point(aes(x=date,y=new_deaths),fill='green',size=2)+
    ylab("confirmed deaths")+
  scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave("images/totaldeaths.png",height = 2,width =6)
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png) ![](./images/totaldeaths.png)

Total number of deaths per province
-----------------------------------

``` r
covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "total") %>% group_by(Province) %>% summarise(total_cases=sum(total)) %>% arrange((total_cases)) %>%
   ggplot()+
  geom_col(aes(x=Province,y=total_cases),fill='orange')+
    ylab("confirmed infections")+
  #scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave("images/deathsperprovinces.png",height = 2,width =6)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png) ![](./images/deathsperprovinces.png)

Top 6 provinces with confirmed cases
------------------------------------

``` r
top_5<-covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "total") %>% group_by(Province) %>% summarise(total_cases=sum(total)) %>%  top_n(5) %>% select(Province)

covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%group_by(date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "new_cases") %>% group_by(Province) %>% mutate(cumulative_cases=cumsum(new_cases)) %>% arrange(desc(cumulative_cases)) %>% filter( Province %in% unlist(as.list(top_5))) %>% 

  ggplot()+
  geom_line(aes(x=date,y=cumulative_cases,color=Province))+
  geom_point(aes(x=date,y=cumulative_cases,color=Province),size=2)+
  #geom_line(aes(x=date,y=total))+
  ylab("cumulative infections per day [log]")+
   scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme_bw()+
  scale_y_continuous(trans='log10')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave("images/cumulativeperprovincestop5.png",height = 2,width =6)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

![](./images/cumulativeperprovincestop5.png)

Infections per day in Top 6 provinces with confirmed cases
----------------------------------------------------------

``` r
top_5<-covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "total") %>% group_by(Province) %>% summarise(total_cases=sum(total)) %>%  top_n(6) %>% select(Province)

covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "total") %>% filter( Province %in% unlist(as.list(top_5))) %>% 
ggplot()+
  geom_col(aes(x=date,y=total))+
  #geom_point(aes(x=date,y=total,color=Province),size=2)+
  #geom_line(aes(x=date,y=total))+
  ylab("new infections per day")+
   scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Province, ncol = 3)+
  ggsave("images/totalprovinces.png",height = 4,width =6)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

![](./images/totalprovinces.png)

Simple SIR model for ARG Forcasting for next 5 days
---------------------------------------------------

Model adapted from [Tim Churches's](https://timchurches.github.io/blog/posts/2020-02-18-analysing-covid-19-2019-ncov-outbreak-data-with-r-part-1) blog.

![](./images/sir_Argentina.png)

![](./images/geseir_Argentina.png)
