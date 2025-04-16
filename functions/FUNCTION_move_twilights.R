move_twilights <- function(df,minutes_different,sun,show_plot){
  
  require(DescTools)
  require(CircStats)
  conv <- 2*pi/1440 #need to treat times as circular, because of issues around 24:00/00:00 hours 
  
  #make a reference table for selecting thresholds (daily change of light per lat/month in min*2 + 15m variance):
  
  #latitudinal resolution in the reference table
  thresholds<-as.data.frame(c(75,65,55,40,25,0,-25,-40,-55,-65,-75))
  colnames(thresholds)<-c("closest_lat")

  #daylength change in minutes by month and latitude ('closest_lat'):
  thresholds$jan<-c(35,27,20,18,16,15,16,18,20,27,35)
  thresholds$feb<-c(35,27,21,19,16,16,16,19,21,27,35)
  thresholds$mar<-c(30,25,22,20,18,17,18,20,22,25,30)
  thresholds$apr<-c(35,27,22,20,18,17,18,20,22,27,35)
  thresholds$may<-c(35,27,21,19,17,16,17,19,21,27,35)
  thresholds$jun<-c(35,27,20,18,16,15,16,18,20,27,35)
  thresholds$jul<-c(35,27,20,18,16,15,16,18,20,27,35)
  thresholds$aug<-c(35,27,21,19,17,16,17,19,21,27,35)
  thresholds$sep<-c(30,25,22,20,18,17,18,20,22,25,30)
  thresholds$oct<-c(35,27,22,20,18,17,18,20,22,27,35)
  thresholds$nov<-c(35,27,21,19,17,16,17,19,21,27,35)
  thresholds$dec<-c(35,27,20,18,16,15,16,18,20,27,35)
  
  df$time<-strftime(df$tFirst, format="%H:%M:%S")
  df$doy <- as.numeric(strftime(df$tFirst, format = "%j"))
  df$minutes<-minute(df$tFirst)+60*hour(df$tFirst)
  df$diffmin_rowi2<-0
  df$diffmin_rowi3<-0
  df$diff_doy2<-0
  df$diff_doy3<-0
  
  ##define  the threshold in minutes for which A and C is considered the same (twilight A will define the threshold):
  #for this we use the reference table above and select the right month and which reference latitude that is closest to the birds' latitude
  #However, that may not be possible in months affected by equinox (e.g March).
  #during an approx equinox period the function use estimated lats based on 10 lats before and after each equinox period 
  #for a few remaining NA's, the average threshold across lats and months is used (22 min)

  latlon <- coord(df$tFirst,df$tSecond,df$type,degElevation=sun)   
  df<-cbind(df,latlon)

  df$lat<-fillMissing(df$lat, span = 1, max.fill = 5)
  #remove lats during approx. spring equinox period (20 Feb - 5 Apr)
  df$lat[yday(df$tFirst)%in%c(51:95)]<-NA
  #remove lats during approx. autumn equinox period (8 Sep - 20 Oct)
  df$lat[yday(df$tFirst)%in%c(251:293)]<-NA
  #fill missing lats with linear lats based on 10 lats before and after each equinox period
  df$lat<-fillMissing(df$lat, span = 10, max.fill = 130)

  for(i in 1:nrow(df[!is.na(df$lat),])){
    df$thresholds[!is.na(df$lat)][i]<-thresholds[which(thresholds$closest_lat== Closest(thresholds[,1], df$lat[!is.na(df$lat)][i])),(month(df$tFirst[!is.na(df$lat)][i])+1)]
    }
  df$thresholds[is.na(df$thresholds)]<-22 
  
  ##find difference between sunrise A, B, and C of (type = 1) 
  rise<-df[df$type==1,]
  rise2<-NULL
  rise2<-rise
  
  for (i in 1:(length(rise$tFirst))) {
    tryCatch({
      rise2$diffmin_rowi2[i] <- as.numeric(abs(rise$minutes[i + 1]-rise$minutes[i]))
      if(rise$minutes[i + 1]<60 & rise$minutes[i]>1380){rise2$diffmin_rowi2[i] <- as.numeric(abs((rise$minutes[i + 1]+1440)-rise$minutes[i]))}
      if(rise$minutes[i + 1]>1380 & rise$minutes[i]<60){rise2$diffmin_rowi2[i] <- as.numeric(abs(rise$minutes[i + 1]-(rise$minutes[i]+1440)))}
      rise2$diffmin_rowi3[i] <- as.numeric(abs(rise$minutes[i + 2]-rise$minutes[i]))
      if(rise$minutes[i + 2]<60 & rise$minutes[i]>1380){rise2$diffmin_rowi3[i] <- as.numeric(abs((rise$minutes[i + 2]+1440)-rise$minutes[i]))}
      if(rise$minutes[i + 2]>1380 & rise$minutes[i]<60){rise2$diffmin_rowi3[i] <- as.numeric(abs(rise$minutes[i + 2]-(rise$minutes[i]+1440)))}
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  k<-1
  row<-1
  for (k in 1:(nrow(rise) - 1)) {
    rise2[row, 9] <- as.numeric(rise[k + 1, 5]-rise[k, 5])
    rise2[row, 10] <- as.numeric(rise[k + 2, 5]-rise[k, 5])
    row <- row + 1
  }
  
  rise2[is.na(rise2)] <- 0
  rise2$change<-FALSE
  rise2$changeto<-0
  
  for(i in 1:((length(rise2$tFirst))-2)){
    if(rise2$diff_doy2[i]==1 & rise2$diff_doy3[i]==2 & rise2$diffmin_rowi2[i]>(((rise2$threshold[i]-15)/2)+minutes_different) & rise2$diffmin_rowi2[i+1]>(((rise2$threshold[i]-15)/2)+minutes_different) & rise2$diffmin_rowi3[i]<rise2$threshold[i])
    {rise2$change[i+1]<-TRUE}else{rise2$change[i+1]<-FALSE}
    if(rise2$diff_doy2[i]==1 & rise2$diff_doy3[i]==2 & rise2$diffmin_rowi2[i]>(((rise2$threshold[i]-15)/2)+minutes_different) & rise2$diffmin_rowi2[i+1]>(((rise2$threshold[i]-15)/2)+minutes_different) & rise2$diffmin_rowi3[i]<rise2$threshold[i])
    {rise2$changeto[i+1]<-((circ.mean(conv*(c(rise2$minutes[i],rise2$minutes[i+2])-1))/conv) + 1440) %% 1440}
  }
  
  #changes to be made in rise have been defined
  ##########################################
  
  ##find difference between sunsetA, B, and C of (type = 2)
  set<-df[df$type==2,]
  set2<-NULL
  set2<-set
  
  for (i in 1:(length(set$tFirst))) {
    tryCatch({
      set2$diffmin_rowi2[i] <- as.numeric(abs(set$minutes[i + 1]-set$minutes[i]))
      if(set$minutes[i + 1]<60 & set$minutes[i]>1380){set2$diffmin_rowi2[i] <- as.numeric(abs((set$minutes[i + 1]+1440)-set$minutes[i]))}
      if(set$minutes[i + 1]>1380 & set$minutes[i]<60){set2$diffmin_rowi2[i] <- as.numeric(abs(set$minutes[i + 1]-(set$minutes[i]+1440)))}
      set2$diffmin_rowi3[i] <- as.numeric(abs(set$minutes[i + 2]-set$minutes[i]))
      if(set$minutes[i + 2]<60 & set$minutes[i]>1380){set2$diffmin_rowi3[i] <- as.numeric(abs((set$minutes[i + 2]+1440)-set$minutes[i]))}
      if(set$minutes[i + 2]>1380 & set$minutes[i]<60){set2$diffmin_rowi3[i] <- as.numeric(abs(set$minutes[i + 2]-(set$minutes[i]+1440)))}
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  k<-1
  row<-1
  for (k in 1:(nrow(set) - 1)) {
    set2[row, 9] <- as.numeric(set[k + 1, 5]-set[k, 5])
    set2[row, 10] <- as.numeric(set[k + 2, 5]-set[k, 5])
    row <- row + 1
  }
  
  set2[is.na(set2)] <- 0
  set2$change<-FALSE
  set2$changeto<-0
  
  for(i in 1:((length(set2$tFirst))-2)){
    if(set2$diff_doy2[i]==1 & set2$diff_doy3[i]==2 & set2$diffmin_rowi2[i]>(((set2$threshold[i]-15)/2)+minutes_different) & set2$diffmin_rowi2[i+1]>(((set2$threshold[i]-15)/2)+minutes_different) & set2$diffmin_rowi3[i]<set2$threshold[i])
    {set2$change[i+1]<-TRUE}else{set2$change[i+1]<-FALSE}
    if(set2$diff_doy2[i]==1 & set2$diff_doy3[i]==2 & set2$diffmin_rowi2[i]>(((set2$threshold[i]-15)/2)+minutes_different) & set2$diffmin_rowi2[i+1]>(((set2$threshold[i]-15)/2)+minutes_different) & set2$diffmin_rowi3[i]<set2$threshold[i])
    {set2$changeto[i+1]<-((circ.mean(conv*(c(set2$minutes[i],set2$minutes[i+2])-1))/conv) + 1440) %% 1440}
  }
  
  
  #changes to be made in set are now defined!
  ##########################################
  changes<-as.data.frame(rbind(set2[,c(1:10,14:15)],rise2[,c(1:10,14:15)]))
  changes$DTime<- as.POSIXct( changes$tFirst , format = "%d/%m/%y %H:%M" , tz = "GMT")
  changes<-changes[ order(changes$DTime , decreasing = FALSE ),]
  changes<-changes[,1:12]
  ########################### 
  
  before_changes<-changes  #save for plot
  
  ## make changes
  changes$hour<-0
  for(i in 1:length(changes$tFirst)){
    if(changes$changeto[i]!=0 & changes$changeto[i]<60){changes$hour[i]<-0}
    if(changes$changeto[i]!=0 & changes$changeto[i]>60){changes$hour[i]<-1}
    if(changes$changeto[i]!=0 & changes$changeto[i]>120){changes$hour[i]<-2}
    if(changes$changeto[i]!=0 & changes$changeto[i]>180){changes$hour[i]<-3}
    if(changes$changeto[i]!=0 & changes$changeto[i]>240){changes$hour[i]<-4}
    if(changes$changeto[i]!=0 & changes$changeto[i]>300){changes$hour[i]<-5}
    if(changes$changeto[i]!=0 & changes$changeto[i]>360){changes$hour[i]<-6}
    if(changes$changeto[i]!=0 & changes$changeto[i]>420){changes$hour[i]<-7}
    if(changes$changeto[i]!=0 & changes$changeto[i]>480){changes$hour[i]<-8}
    if(changes$changeto[i]!=0 & changes$changeto[i]>540){changes$hour[i]<-9}
    if(changes$changeto[i]!=0 & changes$changeto[i]>600){changes$hour[i]<-10}
    if(changes$changeto[i]!=0 & changes$changeto[i]>660){changes$hour[i]<-11}
    if(changes$changeto[i]!=0 & changes$changeto[i]>720){changes$hour[i]<-12}
    if(changes$changeto[i]!=0 & changes$changeto[i]>780){changes$hour[i]<-13}
    if(changes$changeto[i]!=0 & changes$changeto[i]>840){changes$hour[i]<-14}
    if(changes$changeto[i]!=0 & changes$changeto[i]>900){changes$hour[i]<-15}
    if(changes$changeto[i]!=0 & changes$changeto[i]>960){changes$hour[i]<-16}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1020){changes$hour[i]<-17}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1080){changes$hour[i]<-18}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1140){changes$hour[i]<-19}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1200){changes$hour[i]<-20}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1260){changes$hour[i]<-21}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1320){changes$hour[i]<-22}
    if(changes$changeto[i]!=0 & changes$changeto[i]>1380){changes$hour[i]<-23}
  }
  changes$minute<-changes$changeto-(changes$hour*60)
  changes$minute<-round(changes$minute,digits=0)
  changes$second<-format(changes$tFirst, format = "%S")
  changes$year<-format(changes$tFirst, format = "%Y")
  changes$month<-format(changes$tFirst, format = "%m")
  changes$day<-format(changes$tFirst, format = "%d")
  changes$changetodatetime<-paste(changes$year,"-",changes$month,"-",changes$day," ",changes$hour,":",changes$minute,":",changes$second,sep="")
  changes$changetodatetime <- as.POSIXct(changes$changetodatetime, format= "%Y-%m-%d %H:%M:%S")
  changes$tFirst[changes$change == TRUE] <- changes$changetodatetime[changes$change == TRUE]  
  changes$tSecond[1:(length(changes$tFirst)-1)]<-changes$tFirst[2:length(changes$tFirst)]
  
  
  df<-NULL
  df<-changes[,1:3]
  
  
  
  
  #plot
  if(show_plot==TRUE){
    before_changes<-before_changes[!is.na(before_changes$tFirst),]
    before_changes$hours<-as.numeric(format(before_changes[, 1], "%H")) + as.numeric(format(before_changes[,1], "%M"))/60
    
    after_changes<- changes[changes$change==TRUE,]
    after_changes<-after_changes[!is.na(after_changes$changetodatetime),]
    after_changes$hours<-as.numeric(format(after_changes$changetodatetime, "%H")) + as.numeric(format(after_changes$changetodatetime, "%M"))/60
    
    plot(before_changes$tFirst,before_changes$hours,col="grey",pch=19,cex=0.5, ylim=c(0,24),yaxt="n",xaxt="n",ann = FALSE)
    mtext(side = 1, text = "Month", line = 1,cex=0.7)
    mtext(side = 2, text = "Time of day (GMT)", line = 1,cex=0.7)
    points(before_changes$tFirst[before_changes$change==TRUE],before_changes$hours[before_changes$change==TRUE],cex=0.5,col="firebrick",pch=19)
    points(after_changes$tFirst,after_changes$hours,col="cornflowerblue",pch=19,cex=0.5)
    daterange=c(as.POSIXlt(min(before_changes$tFirst)), as.POSIXlt(max(before_changes$tFirst)))
    axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="month"), format="%b",cex.axis=0.6,tck=-0.02,mgp=c(3, 0, 0))
    axis(side=2,at=c(1:24),labels=c(1:24),tck=-0.02,cex.axis=0.6,las=2,mgp=c(3, 0.3, 0))
  }
  
  
  output<-na.omit(df)
  return(output)
}
