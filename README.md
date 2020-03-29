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
    ##  [3] "total_deaths"          "total_imported%"      
    ##  [5] "total_local_contact%"  "total_under_analysis%"
    ##  [7] "new_cases"             "new_deaths"           
    ##  [9] "local_contact"         "imported"             
    ## [11] "under_analysis"        "CABA"                 
    ## [13] "Buenos Aires"          "Chaco"                
    ## [15] "Córdoba"               "Corrientes"           
    ## [17] "Entre Ríos"            "Jujuy"                
    ## [19] "La Pampa"              "Mendoza"              
    ## [21] "Misiones"              "Neuquén"              
    ## [23] "Río Negro"             "Salta"                
    ## [25] "San Luis"              "Santa Cruz"           
    ## [27] "Santa Fe"              "Santiago del Estero"  
    ## [29] "Tierra del Fuego"      "Tucumán"              
    ## [31] "X31"

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

![](README_files/figure-markdown_github/pressure-1.png)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

Infections per day in some provinces
------------------------------------

``` r
covid19_arg[,c(12:27)] %>% replace(is.na(.), 0) %>% tibble::add_column(date=covid19_arg$date) %>%
  pivot_longer(-date,names_to = "Province", values_to = "total") %>% filter( Province %in% c("CABA","Buenos Aires","Córdoba")) %>% 
ggplot()+
  geom_line(aes(x=date,y=total,color=Province))+
  geom_point(aes(x=date,y=total,color=Province),size=4)+
  #geom_line(aes(x=date,y=total))+
  ylab("new infections per day")+
   scale_x_date(date_breaks = "1 day", date_labels = "%d %b")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)
