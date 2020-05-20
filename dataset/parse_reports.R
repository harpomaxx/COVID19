# 
library(pdftools)
library(stringr)
library(pdftools)
library(stringr)
library(lubridate)
library(dplyr)
library(readr)

arreglar_porcentaje<-function(x){
 x %>% str_replace(",",".") %>% as.numeric() %>% "/"(100)
}


dir_reportes<-"/home/harpo/informes"
files <- list.files(path = dir_reportes, pattern = "pdf$")
report<-data.frame()

for(inf_file in files){
  print(inf_file)
  
  inf<-pdf_text(paste0(dir_reportes,"/",inf_file))
  inf<-str_split(inf,'\n')
  provincias=c()
  
  total_importados<-0
  total_contacto<-0
  total_comunitario<-0
  total_descartados<-0
  
  for (line in unlist(inf)){
    line<-gsub("|•","-",line)
  
      if (grepl("^.* hasta ayer es de (\\d+\\.?\\d+) .*$",line,perl=T))
      total_descartados<-gsub("^.* hasta ayer es de (\\d+\\.?\\d+) .*$","\\1",line,perl=T)
    
    if (grepl("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son importados,.*$",line,perl=T))
      total_importados<-gsub("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son importados,.*$","\\1", line,perl=T)
    
    if (grepl("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son contactos.*$",line,perl=T))
      total_contacto<-gsub("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son contactos.*$","\\1", line,perl=T)
    
    if (grepl("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son casos de circu.*$",line,perl=T))
      total_comunitario<-gsub("^.*[ ]*\\((\\d+\\,?.?\\d*)%\\) son casos de circu.*$","\\1", line,perl=T)
    
    if (grepl("^.* es de (\\d+\\.?\\d+), .*$",line,perl=T))
      total_infectados<-gsub("^.* es de (\\d+\\.?\\d+), .*$","\\1",line,perl=T)
    
    if (grepl("^.* (\\d+\\.?\\d+) fallecieron.*$",line,perl=T))
      total_fallecidos<-gsub("^.* (\\d+\\.?\\d+) fallecieron.*$","\\1",line,perl=T)
    
    if (grepl("^.*altas es de (\\d+\\.?\\d+) .*$",line,perl=T))
      total_recuperados<-gsub("^.*altas es de (\\d+\\.?\\d+) .*$","\\1",line,perl=T)
    
    if (grepl("^.*confirmados (\\d+\\.?\\d+) .*$",line,perl=T))
      total_casos_nuevos<-gsub("^.*confirmados (\\d+\\.?\\d+) .*$","\\1", line,perl=T)
    
    
    #provincias
    if (grepl("^[ ]*(-|•|)?[ ]*([a-zA-Z áéíóú]+) ((-)?\\d+\\.?\\d*)(\\*+)? (\\||\\/) .*$", line,perl=T)){
      provincia<-gsub("^[ ]*(-|•|)?[ ]*([a-zA-Z áéíóú]+) ((-)?\\d+\\.?\\d*)(\\*+)? (\\||\\/) .*$","\\2|\\3" , line,perl=T) %>% str_split("\\|",simplify = T)
      if (grepl(provincia[1],"Aires"))
        provincia[1]<-"Buenos Aires"
      provincias<-rbind(provincias,cbind(provincia=provincia[1],casos_nuevos=provincia[2]))
      
    }
  }
  
  provincias<-t(provincias)
  provincias_col_names<-t(provincias[1,])
  casos_nuevos <-t(as.numeric(provincias[2,]))
  provincias<-casos_nuevos
  colnames(provincias)<-t(provincias_col_names)
  provincias<-provincias %>% as.data.frame()
  
  new_report<-data.frame(file=inf_file,total_infectados=total_infectados %>% str_remove("\\."),
                         total_fallecidos=total_fallecidos %>% str_remove("\\."),
                         total_recuperados=total_recuperados %>% str_remove("\\."),
                         total_casos_nuevos=total_casos_nuevos %>% str_remove("\\."),
                         total_importados=total_importados %>% arreglar_porcentaje(),
                         total_contacto=total_contacto %>% arreglar_porcentaje(),
                         total_comunitario=total_comunitario %>% arreglar_porcentaje(),
                         total_descartados=total_descartados %>% str_remove("\\."))
  new_report<-cbind(new_report,provincias)
  report<-rbind(report,new_report)
}
report<- report %>% mutate(date=dmy(str_sub(file,1,8))-1) %>% arrange(date) %>% readr::write_csv(paste0(dir_reportes,"/covid19_argentina.csv"))
View(report)
