
#df<- twl
#breedingloc_lon <- col_lon
#breedingloc_lat <- col_lat
#days_breeding <- c(105:258)
#months_breeding <- c(5:8)
#months_extensive_legtucking<-NA
#show_plot<-TRUE



twilight_cleanup2 <- function(df,breedingloc_lon,breedingloc_lat,days_breeding,months_extensive_legtucking,show_plot){
  
  #datetime_conversion
  df$time<-strftime(df$tFirst, format="%H:%M:%S")
  df$time_mins<-as.numeric(difftime(as.POSIXct(df$time, format = '%H:%M:%S'), as.POSIXct('00:00:00', format = '%H:%M:%S'), units = 'min'))
  
  # find dates with too many twilights
  ok=NULL
  ok<-df
  ok$remove<-FALSE
  ok$date<-as.Date(df$tFirst, format="%Y-%mm-%dd")
  for(i in 1:length(ok$time[ok$type==1])){
    tryCatch({
      if(ok$date[ok$type==1][i]==ok$date[ok$type==1][i+1]){ok$remove[ok$type==1][i]<-TRUE}
      if(ok$date[ok$type==1][i]==ok$date[ok$type==1][i+1]){ok$remove[ok$type==1][i+1]<-TRUE}
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}
  for(i in 1:length(ok$time[ok$type==2])){
    tryCatch({
      if(ok$date[ok$type==2][i]==ok$date[ok$type==2][i+1]){ok$remove[ok$type==2][i]<-TRUE}
      if(ok$date[ok$type==2][i]==ok$date[ok$type==2][i+1]){ok$remove[ok$type==2][i+1]<-TRUE}
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}
  
  
  #allow dates with two sunrises, if more than 22 hours apart from each other
  rise<-ok[ok$type==1,]
  for(i in 1:length(rise$tFirst)){return_to_false<-rise[rise$date==rise$date[i],]
  if(rise$remove[i]==TRUE & nrow(return_to_false)==2 & (difftime(return_to_false$tFirst[2],return_to_false$tFirst[1], units = "hours"))>22){  rise$remove[rise$date==rise$date[i]]<-FALSE} }
  
  #allow dates with two sunsets, if more than 22 hours apart from each other
  set<-ok[ok$type==2,]
  for(i in 1:length(set$tFirst)){return_to_false<-set[set$date==set$date[i],]
  if(set$remove[i]==TRUE & nrow(return_to_false)==2 & (difftime(return_to_false$tFirst[2],return_to_false$tFirst[1], units = "hours"))>22){ set$remove[set$date==set$date[i]]<-FALSE} }
  
  
  
  ########################### 
  together<-NULL
  together<-as.data.frame(rbind(set,rise))
  together$DTime <- as.POSIXct(together$tFirst , format = "%d/%m/%y %H:%M" , tz = "GMT")
  ok<-NULL
  ok<-together[ order(together$DTime , decreasing = FALSE ),]
  ok2<-ok[,1:6]
  ###########################
  
  ##predict time of sunsets and sunrises:
  #keep unlikely timed twilights out of the predictions with use of loess and standard deviation
  #NB! if the dataset is too short, this can create error when finding outliers with loess and will have no effect!
  #loess-filtering:
  ok2$loess_filter<-TRUE
  tryCatch({
    ok2$loess_filter[ok2$remove==FALSE] <-loessFilter(ok2$tFirst[ok2$remove==FALSE],ok2$tSecond[ok2$remove==FALSE],ok2$type[ok2$remove==FALSE],k=3,plot=F) 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  #Standard deviation: calculate SD every 5th day to keep points with SD> 60mins out when making predictions:
  ok2$sd[ok2$type==1]<-roll_sd((hour(ok2$tFirst[ok2$type==1])*60)+minute(ok2$tFirst[ok2$type==1]), 5, weights = rep(1, 5), center = TRUE,
                               min_obs = 1, complete_obs = FALSE, na_restore = FALSE,
                               online = TRUE)  
  ok2$sd[ok2$type==2]<-roll_sd((hour(ok2$tFirst[ok2$type==2])*60)+minute(ok2$tFirst[ok2$type==2]), 5, weights = rep(1, 5), center = TRUE,
                               min_obs = 1, complete_obs = FALSE, na_restore = FALSE,
                               online = TRUE)  
  ok2$sd[is.na(ok2$sd)]<-10000
  ok2$loess_filter[ok2$remove==FALSE & !is.na(ok2$sd)] <-ok2$sd[ok2$remove==FALSE & !is.na(ok2$sd)]<60 
  
  
  #predict twilights  
  ok3<-ok2
  
  #circular times
  hours<-as.numeric(format(ok3[, 1], "%H")) + as.numeric(format(ok3[,1], "%M"))/60
  for (t in 1:2) {
    cor <- rep(NA, 24)
    for (i in 0:23) {cor[i + 1] <- max(abs((c(hours[ok3$type == t][1], hours[ok3$type == t]) + i)%%24 - (c(hours[ok3$type == t], hours[ok3$type == t][length(hours)]) + i)%%24), na.rm = T)}
    hours[ok3$type == t] <- (hours[ok3$type == t] + (which.min(round(cor,2))) - 1)%%24}
  
  #NB! if the dataset is too short, this can create error in loess prediction and this filter will have no effect!
  ok3$hours<-hours
  ok3$predict<-NA
  for (d in seq(30, 1, length = 5)) {tryCatch({
    ok3$predict[ok3$type == 1 & ok3$loess_filter & ok3$remove == FALSE ] <- predict(loess(ok3$hours[ok3$type == 1 & ok3$loess_filter & ok3$remove == FALSE]~as.numeric(ok3$tFirst[ok3$type == 1 & ok3$loess_filter & ok3$remove == FALSE]), span = 0.05))
    ok3$predict[ok3$type == 2 & ok3$loess_filter & ok3$remove == FALSE ] <- predict(loess(ok3$hours[ok3$type == 2 & ok3$loess_filter & ok3$remove == FALSE]~as.numeric(ok3$tFirst[ok3$type == 2 & ok3$loess_filter & ok3$remove == FALSE]), span = 0.05)) 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}
  
  # instruct max.fill to fill by number of dates, not number of rows:
  fillMisspred1<-ok3[c(1,3,10)]
  fillMisspred1<-fillMisspred1[fillMisspred1$type%in%1,]
  fillMisspred1<-fillMisspred1[!duplicated(as.Date(fillMisspred1$tFirst)),]
  
  fillMisspred2<-ok3[c(1,3,10)]
  fillMisspred2<-fillMisspred2[fillMisspred2$type%in%2,]
  fillMisspred2<-fillMisspred2[!duplicated(as.Date(fillMisspred2$tFirst)),]
  
  fillMisspred1<-as.data.frame(fillMisspred1 %>%
                                 mutate(Date = as.Date(tFirst)) %>%
                                 mice::complete(Date = seq.Date(min(as.Date(fillMisspred1$tFirst)), max(as.Date(fillMisspred1$tFirst)), by="day")))
  fillMisspred2<-as.data.frame(fillMisspred2 %>%
                                 mutate(Date = as.Date(tFirst)) %>%
                                 mice::complete(Date = seq.Date(min(as.Date(fillMisspred2$tFirst)), max(as.Date(fillMisspred2$tFirst)), by="day")))
  
  #fill short holes in the dataset:
  fillMisspred1$predict<-fillMissing(fillMisspred1$predict, span = 1, max.fill = 5)
  fillMisspred2$predict<-fillMissing(fillMisspred2$predict, span = 1, max.fill = 5)
  
  
  #define non-breeding months. For practicalities, add first and last month of breeding to the "months_nonbreeding"
  days_nonbreeding<-1:366
  if(length(days_breeding)>0){
    days_nonbreeding<-days_nonbreeding[!(days_nonbreeding%in%days_breeding)]
    days_nonbreeding=c(days_nonbreeding,min(days_breeding),max(days_breeding))
  }
  
  #fill larger holes in the dataset up to 20 days during non-breeding:
  fillMisspred1$predict[yday(fillMisspred1$Date)%in%days_nonbreeding]<-fillMissing(fillMisspred1$predict[yday(fillMisspred1$Date)%in%days_nonbreeding], span = 10, max.fill = 20)
  fillMisspred2$predict[yday(fillMisspred2$Date)%in%days_nonbreeding]<-fillMissing(fillMisspred2$predict[yday(fillMisspred2$Date)%in%days_nonbreeding], span = 10, max.fill = 20)
  
  fillMisspred1<-fillMisspred1[!is.na(fillMisspred1$type),]
  fillMisspred2<-fillMisspred2[!is.na(fillMisspred2$type),]
  
  #add predictions to working table
  #i=1
  for(i in 1:length(unique(as.Date(fillMisspred1$tFirst))))
  {ok3$predict[ok3$type%in%1 & as.Date(ok3$tFirst)%in%unique(as.Date(fillMisspred1$tFirst))[i]]<-fillMisspred1$predict[as.Date(fillMisspred1$tFirst)%in%unique(as.Date(fillMisspred1$tFirst))[i]]}
  #i=1
  for(i in 1:length(unique(as.Date(fillMisspred2$tFirst))))
  {ok3$predict[ok3$type%in%2 & as.Date(ok3$tFirst)%in%unique(as.Date(fillMisspred2$tFirst))[i]]<-fillMisspred2$predict[as.Date(fillMisspred2$tFirst)%in%unique(as.Date(fillMisspred2$tFirst))[i]]}
  
  
  ok3$hours2<-as.numeric(format(ok3[, 1], "%H")) + as.numeric(format(ok3[,1], "%M"))/60
  ok3$test<-ok3$hours2-ok3$hours
  ok3$date<-as.Date(ok3$tFirst, format="%Y-%mm-%dd")
  
  
  #When breeding months are specified, find time of civil dusk and dawn for the breeding location during the year. 
  
  if(length(days_breeding)>0){
    sun_data<-getSunlightTimes(date = ok3$date, lat = breedingloc_lat, lon = breedingloc_lon,
                               keep = c("dawn", "dusk"), tz = "GMT") #tz changed from "UTC" to "GMT"
    sun_data<-as.data.frame(sun_data)
    colnames(sun_data)<-c("date","lat","lon","dawn","dusk")
    hours_dusk<-as.numeric(format(sun_data$dusk, "%H")) + as.numeric(format(sun_data$dusk, "%M"))/60
    hours_dawn<-as.numeric(format(sun_data$dawn, "%H")) + as.numeric(format(sun_data$dawn, "%M"))/60
    ok3<-cbind(ok3,hours_dusk,hours_dawn)
    ok3$hours_dusk<-as.POSIXct(paste(as.Date(ok3$tFirst)," ",as.character(floor(ok3$hours_dusk)),":",as.character(floor((ok3$hours_dusk - floor(ok3$hours_dusk))*60)),":00",sep=""),format="%Y-%m-%d %H:%M:%S", tz="UTC")
    ok3$hours_dawn<-as.POSIXct(paste(as.Date(ok3$tFirst)," ",as.character(floor(ok3$hours_dawn)),":",as.character(floor((ok3$hours_dawn - floor(ok3$hours_dawn))*60)),":00",sep=""),format="%Y-%m-%d %H:%M:%S", tz="UTC")
  }
  
  
  ok3$predict<-abs(ok3$predict+ok3$test)
  ok3$predict[!is.na(ok3$predict)& ok3$predict>24]<-24-(ok3$predict[!is.na(ok3$predict)& ok3$predict>24]-24)
  
  ok3$predict2<-as.POSIXct(paste(as.Date(ok3$tFirst)," ",as.character(floor(ok3$predict)),":",as.character(floor((ok3$predict - floor(ok3$predict))*60)),":00",sep=""),format="%Y-%m-%d %H:%M:%S", tz="UTC")
  ok3$fillMissing_predict<-as.POSIXct(paste(as.Date(ok3$tFirst)," ",as.character(floor(ok3$predict)),":",as.character(floor((ok3$predict - floor(ok3$predict))*60)),":00",sep=""),format="%Y-%m-%d %H:%M:%S", tz="UTC")
  
  
  #when breeding months are specified use the breeding loactions' civil (sun -6) dawn and dusk as predicted twilights during the summer if predict/predict2 is NA:
  if(length(days_breeding)>0){
    ok3$predict2[yday(ok3$date)%in%days_breeding & is.na(ok3$predict2) & ok3$type==1]<-ok3$hours_dawn[yday(ok3$date)%in%days_breeding & is.na(ok3$predict2) & ok3$type==1]
    ok3$predict2[yday(ok3$date)%in%days_breeding & is.na(ok3$predict2) & ok3$type==2]<-ok3$hours_dusk[yday(ok3$date)%in%days_breeding & is.na(ok3$predict2) & ok3$type==2]
  }
  
  
  
  #################################
  #' This function returns the time in 'timeVector' that is ' closest to 'time'
  closest.time <- function(timeVector, time) {
    times()
    x <- times(timeVector)
    v <- times(time)
    clockwise_distance = abs(x - v) 
    anticlockwise_distance = times("23:59:59") - clockwise_distance + times("00:00:01")
    clockwise_and_anticlockwise <-  matrix(c(anticlockwise_distance,  clockwise_distance), ncol = 2)
    shortest_distance_of_the_two <- apply(clockwise_and_anticlockwise, 1, min)
    indx <- which(shortest_distance_of_the_two == min(shortest_distance_of_the_two))
    x[indx] }
  ################################
  
  ##retain the candidate twilight that is closest to a predicted twilight, remove the others
  
  check<-ok3[!is.na(ok3$predict2) & ok3$type==1,]
  #check<-ok3[ok3$remove==FALSE & ok3$type==1,]
  rise<-ok3[ok3$type==1,]
  
  
  i=1
  for(i in 1:length(rise$tFirst)){
    keep = NULL
    tFirst = NULL
    #if(rise$remove[i]==TRUE){
    x = as.Date(rise$date[i])
    cand_to_keep<-rise[rise$date==rise$date[i],]
    cand_to_keep$time<-strftime(cand_to_keep$tFirst, format="%H:%M:%S")
    
    #in months with extensive leg tucking due to molt, candidates are allowed to be compared to a reference time from other dates up to 10 days 
    #(because predictions might have failed to be calculated for a period due to light disturbance, but the bird keeps to the same area all the time) 
    
    if(nrow(check)>0){
      compare_to_this<-check[which(abs(check$date-x) == min(abs(check$date-x)))[1],]
      compare_to_this$comp_predict_time<-strftime(compare_to_this$predict2, format="%H:%M:%S")
      if(month(x)%in%months_extensive_legtucking & is.na(cand_to_keep$predict2[1]) & as.numeric(compare_to_this$date-x)<11){
        keep<-closest.time(cand_to_keep$time,compare_to_this$comp_predict_time)
        to_keep<-cand_to_keep[cand_to_keep$time==keep[1],]
        if(!is.na(compare_to_this$comp_predict_time)){
          if(rise$remove[i]==TRUE & rise$tFirst[i]==to_keep$tFirst){
            rise$remove[i]<-FALSE}
        }}}
    #for all other months of the year: candidates can only be compared to times within the same date. 
    compare_to_this<-cand_to_keep[1,]
    compare_to_this$comp_predict_time<-strftime(compare_to_this$predict2, format="%H:%M:%S")
    keep<-closest.time(cand_to_keep$time,compare_to_this$comp_predict_time)
    to_keep<-cand_to_keep[cand_to_keep$time==keep[1],]
    if(!is.na(compare_to_this$comp_predict_time)){
      if(rise$remove[i]==TRUE & rise$tFirst[i]==to_keep$tFirst){
        rise$remove[i]<-FALSE}
    }}
  #}
  false_rise_removed<-rise[rise$remove==FALSE,]
  #now only one rise remain for each day
  
  check<-ok3[!is.na(ok3$predict2) & ok3$type==2,]
  #check<-ok3[ok3$remove==FALSE & ok3$type==2,]
  set<-ok3[ok3$type==2,] #rise
  
  
  i=1
  for(i in 1:length(set$tFirst)){
    keep = NULL
    tFirst = NULL
    x = as.Date(set$date[i])
    cand_to_keep<-set[set$date==set$date[i],]
    cand_to_keep$time<-strftime(cand_to_keep$tFirst, format="%H:%M:%S")
    
    #in months with extensive leg tucking due to molt, candidates are allowed to be compared to a reference time from other dates up to 10 days 
    #(because predictions might have failed to be calculated for a period due to light disturbance, but the bird keeps to the same area all the time) 
    
    if(nrow(check)>0){
      compare_to_this<-check[which(abs(check$date-x) == min(abs(check$date-x)))[1],]
      compare_to_this$comp_predict_time<-strftime(compare_to_this$predict2, format="%H:%M:%S")
      if(month(x)%in%months_extensive_legtucking & is.na(cand_to_keep$predict2[1]) & as.numeric(compare_to_this$date-x)<11){
        keep<-closest.time(cand_to_keep$time,compare_to_this$comp_predict_time)
        to_keep<-cand_to_keep[cand_to_keep$time==keep[1],]
        if(!is.na(compare_to_this$comp_predict_time)){
          if(set$remove[i]==TRUE & set$tFirst[i]==to_keep$tFirst){
            set$remove[i]<-FALSE}
        }}}
    #for all other months of the year: candidates can only be compared to times within the same date.       
    compare_to_this<-cand_to_keep[1,]
    compare_to_this$comp_predict_time<-strftime(compare_to_this$predict2, format="%H:%M:%S")
    keep<-closest.time(cand_to_keep$time,compare_to_this$comp_predict_time)
    to_keep<-cand_to_keep[cand_to_keep$time==keep[1],]
    if(!is.na(compare_to_this$comp_predict_time)){
      if(set$remove[i]==TRUE & set$tFirst[i]==to_keep$tFirst){
        set$remove[i]<-FALSE}
    }}
  false_set_removed<-set[set$remove==FALSE,]
  #now only one set remain for each day
  
  
  
  ########################### 
  df2<-NULL
  df2<-as.data.frame(rbind(false_set_removed[,c(1:3)],false_rise_removed[,c(1:3)]))
  df2<-df2[ order(df2$tFirst , decreasing = FALSE ),]
  df2$tSecond[1:(length(df2$tSecond)-1)]<-df2$tFirst[2:length(df2$tSecond)]
  df<-NULL
  df<-df2
  
  ###########################
  
  #plot
  if(show_plot==TRUE){
    base<-ok
    chosen<-ok3
    chosen$remove2<-NA
    chosen$remove2[chosen$type%in%2]<-set$remove
    chosen$remove2[chosen$type%in%1]<-rise$remove
    chosen<-chosen[chosen$remove!=chosen$remove2,]
    chosen$hours<-as.numeric(format(chosen[, 1], "%H")) + as.numeric(format(chosen[,1], "%M"))/60
    base$hours<-as.numeric(format(base[, 1], "%H")) + as.numeric(format(base[,1], "%M"))/60
    result<-df
    result$hours<-as.numeric(format(result[, 1], "%H")) + as.numeric(format(result[,1], "%M"))/60
    set$hours<-as.numeric(format(set[, 1], "%H")) + as.numeric(format(set[,1], "%M"))/60
    rise$hours<-as.numeric(format(rise[, 1], "%H")) + as.numeric(format(rise[,1], "%M"))/60
    predicted_line<-ok3
    predicted_line$fillMissing_predict<-as.numeric(format(predicted_line$fillMissing_predict, "%H")) + as.numeric(format(predicted_line$fillMissing_predict, "%M"))/60
    predicted_line$predict2<-as.numeric(format(predicted_line$predict2, "%H")) + as.numeric(format(predicted_line$predict2, "%M"))/60
    test<-result[as.Date(result$tFirst)%in%as.Date(base$tFirst[!(base$tFirst%in%result$tFirst)]),]
    plot(base$tFirst,base$hours,col="grey",cex=0.5,pch=19, ylim=c(0,24),yaxt="n",xaxt="n",ann = FALSE)
    mtext(side = 1, text = "Month", line = 1,cex=0.7)
    mtext(side = 2, text = "Time of day (GMT)", line = 1,cex=0.7)
    points(rise$tFirst[rise$remove%in%TRUE],rise$hours[rise$remove%in%TRUE],col="firebrick",pch=19,cex=0.5)
    points(set$tFirst[set$remove%in%TRUE],set$hours[set$remove%in%TRUE],col="firebrick",pch=19,cex=0.5)
    points(chosen$tFirst,chosen$hours,col="cornflowerblue",pch=19,cex=0.5)
    lines(predicted_line$tFirst[predicted_line$type==1&predicted_line$remove==FALSE], predicted_line$predict2[predicted_line$type==1&predicted_line$remove==FALSE], type = "l",col="chartreuse3",lwd=2,lty=1,)
    lines(predicted_line$tFirst[predicted_line$type==2&predicted_line$remove==FALSE], predicted_line$predict2[predicted_line$type==2&predicted_line$remove==FALSE], type = "l",col="chartreuse3",lwd=2,lty=1)
    
    lines(predicted_line$tFirst[predicted_line$type==1&predicted_line$remove==FALSE], predicted_line$fillMissing_predict[predicted_line$type==1&predicted_line$remove==FALSE], type = "l",col="orange",lwd=2,)
    lines(predicted_line$tFirst[predicted_line$type==2&predicted_line$remove==FALSE], predicted_line$fillMissing_predict[predicted_line$type==2&predicted_line$remove==FALSE], type = "l",col="orange",lwd=2)
    
    daterange=c(as.POSIXlt(min(base$tFirst)), as.POSIXlt(max(base$tFirst)))
    axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="month"), format="%b",cex.axis=0.6,tck=-0.02,mgp=c(3, 0, 0))
    axis(side=2,at=c(1:24),labels=c(1:24),tck=-0.02,cex.axis=0.6,las=2,mgp=c(3, 0.3, 0))
  }
  
  output<-df
  return(output)
}
