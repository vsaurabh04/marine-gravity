library(ncdf4)
library(lubridate)
library(dplyr)
library(geosphere)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(gstat)
library(matlib)
library(sp)
library(raster)


read_input <- function(file_name){
  df = NULL
  newlist <- list.files(file_name)
  
  index <- which(newlist=="desktop.ini")
  if (length(index)!=0){
    newlist <- newlist[-index]
  }
  
  for (i in seq_along(newlist)){
    nc <- nc_open(file.path(file_name,newlist[i]))
    lat <- ncvar_get(nc, "glat.00")
    lon <- ncvar_get(nc, "glon.00")
    ssh <- ncvar_get(nc, "ssh.33")
    jdn <- ncvar_get(nc, "jday.00")
    date<-as.Date(jdn+0.5, origin=as.Date("2000-01-01"))
    cycleNo <- ncatt_get(nc, 0, "cycle")
    passNo <- ncatt_get(nc, 0, "pass_number")
    size <- dim(ssh)
    nc_close(nc)
    rbind(df, data.frame(lat, lon, ssh, jdn, date, size, cycleNo, passNo))->df
  }
  names(df)[c(8,10)]<-c("cycleNo", 'passNo')
  
  df[,c(6,7,9)]<-NULL # some unwated inclusions
  
  return(df)
}

# group the readings by pass number and then group them cycle wise for each pass

group_pairs <- function(df){
  byPass <- df %>% group_by(passNo) %>% group_split() 
  byCycle <- lapply(byPass, function(x){a <- x %>% group_by(cycleNo)
  a<- group_split(a)})
  
  # cycle with minimum number of data is chosen as reference cycle so that we find closer point pairs
  ref_cyc <- lapply(byCycle, function(z){
    inds <- which.min(sapply(z, nrow))
    ref <- z[[inds]]; ref}) 
  
  # rest points
  other_cyc <- lapply(byCycle, function(z){
    inds <- which.min(sapply(z, nrow))
    other <- z[-inds];other
  })
  
  
  # Now lets find the pair of points close to each other
  
  pairs <- NULL
  
  for (i in 1:length(ref_cyc)){
    for (j in 1:length(other_cyc[[i]])){
      dist_mat <- distm(matrix(c(ref_cyc[[i]]$lon, ref_cyc[[i]]$lat), ncol=2), matrix(c(other_cyc[[i]][[j]]$lon, other_cyc[[i]][[j]]$lat), ncol=2), fun=distGeo)/1000
      # filter based on the dist_mat
      lat<- other_cyc[[i]][[j]]$lat[max.col(-dist_mat)]
      lon<- other_cyc[[i]][[j]]$lon[max.col(-dist_mat)]
      jdn<- other_cyc[[i]][[j]]$jdn[max.col(-dist_mat)]
      ssh<- other_cyc[[i]][[j]]$ssh[max.col(-dist_mat)]
      cycleNo <- other_cyc[[i]][[j]]$cycleNo[max.col(-dist_mat)]
      passNo <- other_cyc[[i]][[j]]$passNo[max.col(-dist_mat)]
      
      rbind(pairs, data.frame(lon, lat, ssh, jdn, cycleNo, passNo))->pairs
    }
  }
  
  
  pair_points_byPass <- pairs %>% group_by(passNo) %>% group_split()
  pair_points_byCycle <- lapply(pair_points_byPass, function(x){
    a <- x %>% group_by(cycleNo)
    a<- group_split(a)
  })
  
  # Grouping all the common points to find their time series
  final.pairs <- list()
  for (i in 1:length(pair_points_byCycle)){
    my.list <- pair_points_byCycle[[i]]
    
    l1 <- map_df(.x = seq(1:length(my.list)),
                 .f = function(x) bind_rows(my.list[x]) %>%
                   mutate(id = row_number(),
                          df = x)) %>% arrange(id)
    out <- split(l1, f = l1$id)
    final.pairs[[i]] <- out
  } 
  
  return(final.pairs)
}

#calculate dovs
calc_dovs <- function(final.pairs){
  new_df =  NULL
  for (i in 1:length(final.pairs)){
    for (j in 1:length(final.pairs[[i]])){
      curr_grp <- final.pairs[[i]][[j]]
      dist_mat <- distm(matrix(c(curr_grp$lon, curr_grp$lat), ncol=2),matrix(c(curr_grp$lon, curr_grp$lat), ncol=2), fun=distGeo)/1000
      delta_h <- outer(curr_grp$ssh, curr_grp$ssh, '-')
      distance <-  dist_mat[upper.tri(dist_mat)]
      delta_h <- delta_h[upper.tri(delta_h)]
      
      L <- delta_h/distance   #observation matrix 
      
      azimuth_function <- Vectorize(function(a, b){
        phi_mean <- (curr_grp$lat[a] + curr_grp$lat[b])/2
        delta_lambda <- curr_grp$lon[a] - curr_grp$lon[b]
        delta_phi <- curr_grp$lat[a] - curr_grp$lat[b]
        
        atan((cos(phi_mean*delta_lambda))/delta_phi)
      })
      
      azimuth <- outer(seq_len(dim(curr_grp)[1]), seq_len(dim(curr_grp)[1]), FUN = azimuth_function)
      azimuth <- azimuth[upper.tri(azimuth)]
      
      A <- matrix(c(cos(azimuth), sin(azimuth)), ncol=2)   # design matrix 
      
      # Least Square 
      X <- inv(t(A)%*%A) %*% (t(A) %*% L)   # to find out the value of eta and zeta
      
      mean_lon <- mean(curr_grp$lon)
      mean_lat <- mean(curr_grp$lat)
      mean_ssh <- mean(curr_grp$ssh)
      eta <- X[1]
      zeta <- X[2]
      rbind(new_df, data.frame(lon=mean_lon, lat=mean_lat, ssh=mean_ssh, eta, zeta)) -> new_df
    }
  }
  
  return(new_df)
}

# Code for kriging interpolation here 
kriging <- function(df, grid_size)
{
  # construction of spatial grid  of grid_size X grid_size degree
  grid_size <- 0.5
  lon_grids <- 1 + ((max_lon - min_lon)/grid_size)
  lat_grids <- 1 + ((max_lat - min_lat)/grid_size)
  grd <- data.frame(lon = rep(seq(min_lon, max_lon, grid_size), lon_grids), lat = rep(seq(min_lat, max_lat, grid_size), each = lat_grids))
  coordinates(grd) <- ~lon+lat
  grd <- as.data.frame(spsample(grd, "regular", n = (lon_grids-1)*(lat_grids-1), offset = c(grid_size, grid_size)))
  names(grd) <- c("lon", "lat")
  coordinates(grd) <- c("lon", "lat")
  gridded(grd) <- TRUE  # Create SpatialPixel object
  fullgrid(grd) <- TRUE # Create SpatialGrid object
  proj4string(grd) <- "+init=epsg:4326"
  # Add projection information to the grid
  
  # variogram on the de - trended data.
  coordinates(df) <- ~lon + lat # convert data frame to SPDF
  
  # Add projection information to the processed.data
  crs(df) <- "+init=epsg:4326"
  
  var_zeta <- variogram(zeta~1, data = df, cloud = FALSE)
  
  # to fit.variogram() via the vgm() function
  dat_fit_zeta <- fit.variogram(var_zeta, vgm(c("Exp", "Mat", "Sph")))
  
  # The following plot allows us to assess the fit
  plot(var_zeta, dat_fit_zeta)
  
  f1 <- as.formula(zeta~lon + lat)
  zeta_krg <- krige(f1, df, grd, dat_fit_zeta)
  
  # 
  var_eta <- variogram(eta~1, data=df, cloud=FALSE)
  dat_fit_eta <- fit.variogram(var_eta,  vgm(c("Exp", "Mat", "Sph")))
  plot(var_eta, dat_fit_eta) # to visualise the fit 
  
  f2 <- as.formula(eta~lon + lat)
  eta_krg <- krige(f2, df, grd, dat_fit_eta)
  
  interpolated_grd <- data.frame(grd, zeta_krg@data, eta_krg@data)
  
  # interpolated_grd[,c("col num here")]<-NULL # some unwated inclusions
  
  colnames(interpolated_grd) <- c("lon", "lat", "pred_zeta", "zeta.var", "pred_eta", "eta.var")
  return(interpolated_grd)
}



path <- "C:/Users/Vikram/Desktop/ce678/marine-gravity/data"   # add the path to your dataset folder here(SARAL or JASON-2)
# path <- "C:/Users/Vikram/Documents/data_jason2/all"  


df <- read_input(path)
ggplot()+geom_point(data = df, aes(lon, lat), color="blue") #to visualise the data 

# figuring out the the range of latitude and longitude for the data
min_lat <- min(df$lat)
max_lat <- max(df$lat)

min_lon <- min(df$lon)
max_lon <- max(df$lon)

final.pairs <- group_pairs(df) # pairs the close points from all the cycles based on pass number



new_df <- calc_dovs(final.pairs)

# to visualise the data 
ggplot()+geom_point(data = new_df, aes(lon, lat), color="black")+scale_y_continuous(breaks=seq(min_lat, max_lat, 0.5))+
  scale_x_continuous(breaks=seq(min_lon, max_lon, 0.5))



grid_size <- 0.5
interpolated_grd <- kriging(new_df, 0.5)


ggplot()+geom_point(data = interpolated_grd, aes(lon, lat), color="black")+scale_y_continuous(breaks=seq(min_lat, max_lat, grid_size))+
  scale_x_continuous(breaks=seq(min_lon, max_lon, grid_size))

# Now we need do the 2D FFT to find out the gravity anomaly 



# estimate normal gravity wrt GRS80 reference system 
# using somigliana-pizetti equation to find normal gravity modified for GRS80

gamma = 9.780327*(1+0.0053024*(sin(interpolated_grd$lon))^2 - 0.0000058*(sin(2*interpolated_grd$lon))^2)

fft_eta <- fft(interpolated_grd$pred_eta)
fft_zeta <- fft(interpolated_grd$pred_zeta)

# some value for wavenumbers taken 
kx <- 5
ky <- 5

f_delta_G <- -(0+1i)*2*pi*(kx*fft_eta+ky*fft_zeta)/sqrt(kx^2 + ky^2)
delta_g <- Re(fft(f_delta_G, inverse=TRUE) / length(f_delta_G))
 
