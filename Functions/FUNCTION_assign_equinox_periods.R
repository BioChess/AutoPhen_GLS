assign_equinox_periods <- function(df,breedingloc_lat){
  
median_lat_autumn_start<-median(df$lat_smooth2[month(df$tFirst)%in%8 & day(df$tFirst) %in% c(25:29)])
if(is.na(median_lat_autumn_start))median_lat_autumn_start<-breedingloc_lat
median_lat_autumn_end<-median(df$lat_smooth2[month(df$tFirst)%in%10 & day(df$tFirst) %in% c(18:22)])
if(is.na(median_lat_autumn_end))median_lat_autumn_end<-median_lat_autumn_start
median_lat_spring_start<-median(df$lat_smooth2[month(df$tFirst)%in%2 & day(df$tFirst) %in% c(20:24)])
if(is.na(median_lat_spring_start))median_lat_spring_start<-median_lat_autumn_end
median_lat_spring_end<-median(df$lat_smooth2[month(df$tFirst)%in%4 & day(df$tFirst) %in% c(15:19)])
if(is.na(median_lat_spring_end))median_lat_spring_end<-median_lat_spring_start

lats<-c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80,-90)
start<-c(-7,-10,-14,-18,-21,-21,-21,-21,-22,-23,-22,-21,-21,-21,-21,-18,-14,-10,-7)
end<-abs(start)
adjust_autumn_eq_by<-c(9,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6,-7,-8,-9)
adjust_spring_eq_by<-c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,8)
equinox_table<-data.frame(lats,start,end,adjust_autumn_eq_by,adjust_spring_eq_by)
equinox_table$autumn_equinox<-as.Date(paste(year(df$tFirst[1]),"-09-23",sep=""))+adjust_autumn_eq_by
equinox_table$spring_equinox<-as.Date(paste(year(df$tFirst[1])+1,"-03-21",sep=""))+adjust_spring_eq_by
equinox_table$start_aut_eq<-equinox_table$autumn_equinox+equinox_table$start
equinox_table$end_aut_eq<-equinox_table$autumn_equinox+end
equinox_table$start_spring_eq<-equinox_table$spring_equinox +equinox_table$start
equinox_table$end_spring_eq<-equinox_table$spring_equinox +end

start_aut<-equinox_table$start_aut_eq[which.min(abs(equinox_table$lats - median_lat_autumn_start))]
end_aut<-equinox_table$end_aut_eq[which.min(abs(equinox_table$lats - median_lat_autumn_end))]
start_spring<-equinox_table$start_spring_eq[which.min(abs(equinox_table$lats - median_lat_spring_start))]
end_spring<-equinox_table$end_spring_eq[which.min(abs(equinox_table$lats - median_lat_spring_end))]

df$eqfilter<-1
df$eqfilter[df$tFirst>=(start_aut+1) & df$tFirst<=(end_aut-1)]<-0
df$eqfilter[df$tFirst>=(start_spring+1) & df$tFirst<=(end_spring-1)]<-0

output<-df$eqfilter
return(output)
}