datetime_conversion <- function(datetime,tz="UTC"){
  
  datetime <- as.character(datetime)
  output <- NULL
  # example 04.04.2013 16:04:00
  if(substr(datetime[1],3,3) == "." &
     substr(datetime[1],6,6) == "." &
     substr(datetime[1],11,11) == " " &
     substr(datetime[1],14,14) == ":" &
     substr(datetime[1],17,17) == ":" &
     nchar (datetime[1])     == 19)  output <- as.POSIXct(strptime(datetime,"%d.%m.%Y %H:%M:%S"),tz=tz)
  
  # example 04/04/2013 16:04:00
  if(is.null(output) &
     substr(datetime[1],3,3) == "/" &
     substr(datetime[1],6,6) == "/" &
     substr(datetime[1],11,11) == " " &
     substr(datetime[1],14,14) == ":" &
     substr(datetime[1],17,17) == ":" &
     nchar (datetime[1])     == 19)  output <- as.POSIXct(strptime(datetime,"%d/%m/%Y %H:%M:%S"),tz=tz)
  
  # example 04/04/2013 16:04
  if(is.null(output) &
     substr(datetime[1],3,3) == "/" &
     substr(datetime[1],6,6) == "/" &
     substr(datetime[1],11,11) == " " &
     substr(datetime[1],14,14) == ":" &
     nchar (datetime[1])     == 16)  output <- as.POSIXct(strptime(datetime,"%d/%m/%Y %H:%M"),tz=tz)
  
  # example 04/04/2013 16:04
  if(is.null(output) &
     substr(datetime[1],3,3) == "/" &
     substr(datetime[1],6,6) == "/" &
     substr(datetime[1],11,11) == " " &
     substr(datetime[1],14,14) == ":" &
     nchar (datetime[1])     == 16)  output <- as.POSIXct(strptime(datetime,"%d/%m/%Y %H:%M"),tz=tz)
  
  # example 04.04.13 16:04:00
  if(is.null(output) &
     substr(datetime[1],3,3) == "." &
     substr(datetime[1],6,6) == "." &
     substr(datetime[1],9,9) == " " &
     substr(datetime[1],12,12) == ":" &
     substr(datetime[1],15,15) == ":" &
     nchar (datetime[1])     == 17)  output <- as.POSIXct(strptime(datetime,"%d.%m.%y %H:%M:%S"),tz=tz)
  
  # example 04/04/13 16:04:00
  if(is.null(output) &
     substr(datetime[1],3,3) == "/" &
     substr(datetime[1],6,6) == "/" &
     substr(datetime[1],9,9) == " " &
     substr(datetime[1],12,12) == ":" &
     substr(datetime[1],15,15) == ":" &
     nchar (datetime[1])     == 17)  output <- as.POSIXct(strptime(datetime,"%d/%m/%y %H:%M:%S"),tz=tz)
  
  # example 04.04.13 16:04
  if(is.null(output) &
     substr(datetime[1],3,3) == "." &
     substr(datetime[1],6,6) == "." &
     substr(datetime[1],9,9) == " " &
     substr(datetime[1],12,12) == ":" &
     nchar (datetime[1])     == 14)  output <- as.POSIXct(strptime(datetime,"%d.%m.%y %H:%M"),tz=tz)
  
  # example 04/04/13 16:04
  if(is.null(output) &
     substr(datetime[1],3,3) == "/" &
     substr(datetime[1],6,6) == "/" &
     substr(datetime[1],9,9) == " " &
     substr(datetime[1],12,12) == ":" &
     nchar (datetime[1])     == 17)  output <- as.POSIXct(strptime(datetime,"%d/%m/%y %H:%M"),tz=tz)
  
  
  return(output)
}