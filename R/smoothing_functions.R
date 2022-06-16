# Moving average 7 days
ma <- function(x, n = 7){
  stats::filter(x, rep(1 / n, n), sides = 2)
  }

# fonction de lissage
fun.smooth<-function(data,cases,time){
  
  smoothed_ma<-ma(data[,cases])
  cor<-0
  # cherche le parametre de lissage qui permet de se rapprocher le plus d'une moyenne mobile de 7 jours
  for (k in seq(0.2,0.7,by=0.1)){
    smoothed_ss<-smooth.spline(x=data[,time],y=data[,cases],spar=k)$y
    cor.data<-cor.test(smoothed_ma,smoothed_ss)$estimate
    if (cor.data>cor){
      cor<-cor.data
      sparam<-k
    }
  }
  smoothed<-smooth.spline(x=data[,time],y=data[,cases],spar=sparam)$y
  
  return(smoothed)
}

# exemple sur un dataframe avec une colonne "hosp" et une colonne "date"
# dat$smooth<-fun.smooth(dat,"hosp","date")