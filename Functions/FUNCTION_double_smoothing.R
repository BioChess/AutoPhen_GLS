double_smoothing <- function(df,sun){
  df$diff.h<-difftime(df$tSecond,df$tFirst,units = "hours")
######First smoothing, essentially compensation for east-west movements [mean lat (lat_smooth1)] and north-south movements [mean lon (lon_smooth1)], noon-midnigths. 
  df<-df[df$diff.h<24,]
  df$date_time<-df$tFirst+(difftime(df$tSecond,df$tFirst, units="hours"))/2                   
  df$date_time <- as.POSIXct(as.character(df$date_time),"UTC")                                         

  latlon <- coord(df$tFirst,df$tSecond,df$type,degElevation=sun)   
  df<-cbind(df,latlon)
  
lat_smooth1<-df$lat
for (i in 2:(length(df$lat))){
  ifelse(df$type[i]!=df$type[i-1] & (difftime(df$tFirst[i],df$tFirst[i-1],units = "hours")<30)
         
         ,lat_smooth1[i]<-(df$lat[i]+df$lat[i-1])/2,lat_smooth1[i]<-df$lat[i]) 
}

lon_smooth1<-df$lon
for (i in 2:(length(df$lon))){
  ifelse(df$type[i]!=df$type[i-1] & (difftime(df$tFirst[i],df$tFirst[i-1],units = "hours")<30)          
         ,lon_smooth1[i]<-(df$lon[i]+df$lon[i-1])/2,lon_smooth1[i]<-df$lon[i]) 
}  

#####################################################################################
####extra loop with regard to positions crossing the Pacific 180 Meridian####
lon_smooth10<-lon_smooth1

for (i in 2:(length(df$lon))){
  ifelse(abs(df$lon[i]-df$lon[i-1])>300 & df$type[i]!=df$type[i-1] & difftime(df$tFirst[i],df$tFirst[i-1],units = "hours")<30
         ,lon_smooth1[i]<-(-1*(180-abs(lon_smooth10[i]))*lon_smooth10[i]/abs(lon_smooth10[i])),lon_smooth1[i]<-lon_smooth10[i])
}
#####################################################################################

df$lon_smooth1<-lon_smooth1                                        
df$lat_smooth1<-lat_smooth1                                        

######Second smoothing, two point moving average of noon-midnigth positions, same as Intiproc and Phillips et al 

meanx<-NA
meany<-NA
meanz<-NA

lat_smooth2<-df$lat_smooth1
lon_smooth2<-df$lon_smooth1

x<-cos((df$lat_smooth1/180)*pi)*cos((df$lon_smooth1/180)*pi)
y<-cos((df$lat_smooth1/180)*pi)*sin((df$lon_smooth1/180)*pi)
z<-sin((df$lat_smooth1/180)*pi)

# 2x Smoothing ala Phillips et al------

NObs<-length(df$date_time)
for(i in 1:(NObs-1)) {
  meanx[i]<-ifelse((df$type[i]!=df$type[i+1] & (difftime(df$tFirst[i+1],df$tFirst[i],units = "hours")<30)),(x[i+1]+x[i])/2, x[i])
  meany[i]<-ifelse((df$type[i]!=df$type[i+1] & (difftime(df$tFirst[i+1],df$tFirst[i],units = "hours")<30)),(y[i+1]+y[i])/2, y[i])
  meanz[i]<-ifelse((df$type[i]!=df$type[i+1] & (difftime(df$tFirst[i+1],df$tFirst[i],units = "hours")<30)),(z[i+1]+z[i])/2, z[i])
  lat_smooth2[i]<-90-(atan2(sqrt(meanx[i]^2+meany[i]^2),meanz[i])*(180/pi))
  lon_smooth2[i]<-1*(90-(atan2(meanx[i],meany[i])*180/pi))
}

##################################################################################################
####smoothing positions between -180 and -90 degrees longitude
lon_smooth3<-lon_smooth2    #copy of lon_smooth2 ,  trigonometry uses -90-270 as scale, but this crashes when lon < -90 .  

#This fixes this problem and positions crossing the Pacific 180 Meridian 
for(i in 1:(NObs-1)) {
  ifelse(lon_smooth3[i]>180,lon_smooth2[i]<-lon_smooth3[i]-360,lon_smooth2[i]<-lon_smooth3[i])}
#################################################################################################

df$lon_smooth2<-lon_smooth2
df$lat_smooth2<-lat_smooth2

output<-df
return(output)
}

