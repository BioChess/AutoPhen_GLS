#### FUNCTIONS SGAT MODELING ####

##### EARTH MASK ####

#This function constructs a gridded representation of the world's land masses for the region 
#delimited by xlim and ylim with a resolution of n cells per degree and creates a look-up function
#that returns NA for locations that fall outside the extent of the grid, 
#otherwise it returns TRUE or FALSE depending whether the point corresponds to land or sea.

#Functions from Raul scripts
distribution.mask <- function(xlim, ylim, n = 4, land = TRUE,shape) {
  r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1], 
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(shape))
  r <- cover(rasterize(elide(shape, shift = c(-360, 0)), r, 1, silent = TRUE), 
             rasterize(shape, r, 1, silent = TRUE),
             rasterize(elide(shape, shift = c(360, 0)), r, 1, silent = TRUE))
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  if (land) 
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

#Functions from my script doing SGAT analysis
earthseaMask <- function(xlim, ylim, n = 4, pacific=FALSE) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(proj4string(shape), avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
  mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
               rasterize(wrld_simpl, r, 1, silent = TRUE), 
               rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  
  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
}


log.prior <- function(p) {
  f <- mask2(p)
  ifelse(f | is.na(f), 0, -10)
}



#### SST MASK ####

sstMask <- function(p) {
  f <- apply(cbind(1:nrow(p)), 1, function(x) sstArray[cbind(.bincode(p[x,2], ybin), .bincode(p[x,1], xbin), x)])
  ifelse(f==0 | is.na(f), -1000, log(f))
}


#### SPEED FILTER ####

sp_fun <- function(x){
  require(fossil)
  # x <- ls[[1]]
  x <- x[order(x$Date_Time),]
  x$Dist1b <- NA
  x$Dist2b <- NA
  x$Time1b <- NA
  x$Time2b <- NA
  x$Speed1b <- NA
  x$Speed2b <- NA
  x$Dist1a <- NA
  x$Dist2a <- NA
  x$Time1a <- NA
  x$Time2a <- NA
  x$Speed1a <- NA
  x$Speed2a <- NA
  x$Qspeed <- NA
  for (z in 3:(nrow(x)-2)) {             # exclude first two and last two positions
    # z = 3 
    # creates empty values for later function
    Dist1b<-numeric(nrow(x)); Dist2b<-numeric(nrow(x))
    Time1b<-numeric(nrow(x)); Time2b<-numeric(nrow(x))
    Speed1b<-numeric(nrow(x)); Speed2b<-numeric(nrow(x))
    Dist1a<-numeric(nrow(x)); Dist2a<-numeric(nrow(x))
    Time1a<-numeric(nrow(x)); Time2a<-numeric(nrow(x))
    Speed1a<-numeric(nrow(x)); Speed2a<-numeric(nrow(x))
    QSpeed<-numeric(nrow(x)) 
    lat <- x$lat; long <- x$lon; dtime <- x$Date_Time
    
    # positions before
    Dist1b[z] <- deg.dist(long[z-1], lat[z-1], long[z], lat[z])  # distance between the point and previous point
    x$Dist1b[z]<- Dist1b[z]
    Dist2b[z] <- deg.dist(long[z-2], lat[z-2], long[z], lat[z])  # distance between the point and second previous point
    x$Dist2b[z]<- Dist2b[z]
    Time1b[z] <- -(as.numeric(difftime(dtime[z-1], dtime[z], units = "hours", tz= "GMT")))
    x$Time1b[z]<-Time1b[z]
    Time2b[z] <- -(as.numeric(difftime(dtime[z-2], dtime[z], units = "hours", tz= "GMT")))
    x$Time2b[z]<-Time2b[z]
    Speed1b[z] <- Dist1b[z]/Time1b[z]
    x$Speed1b[z]<-Speed1b[z]
    Speed2b[z] <- Dist2b[z]/Time2b[z]
    x$Speed2b[z]<-Speed2b[z]
    
    # positions after
    Dist1a[z] <- deg.dist(long[z], lat[z], long[z+1], lat[z+1])  # distance between the point and next point
    x$Dist1a[z]<- Dist1a[z]
    Dist2a[z] <- deg.dist(long[z], lat[z], long[z+2], lat[z+2])  # distance between the point and second next point
    x$Dist2a[z]<- Dist2a[z]
    Time1a[z] <- -(as.numeric(difftime(dtime[z], dtime[z+1], units = "hours", tz= "GMT")))
    x$Time1a[z]<-Time1a[z]
    Time2a[z] <- -(as.numeric(difftime(dtime[z], dtime[z+2], units = "hours", tz= "GMT")))
    x$Time2a[z]<-Time2a[z]
    Speed1a[z] <- Dist1a[z]/Time1a[z]
    x$Speed1a[z]<-Speed1a[z]
    Speed2a[z] <- Dist2a[z]/Time2a[z]
    x$Speed2a[z]<-Speed2a[z]
    v2b <- Speed2b[z]
    v1b <- Speed1b[z]
    v1a <- Speed1a[z]
    v2a <- Speed2a[z]
    v<-numeric(nrow(x))
    v[z]<- sqrt(sum(v2b^2, v1b^2, v1a^2, v2a^2)/4)
    x$Qspeed[z] <- v[z]
  }
  return(x)
}
