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

path <- "C:/Users/Vikram/Desktop/ce678/marine-gravity/data"   # add the path to your dataset folder here 
df = NULL
newlist <- list.files(path)
for (i in seq_along(newlist)){
  nc <- nc_open(file.path(path,newlist[i]))
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

ggplot()+geom_point(data = df, aes(lon, lat), color="blue") #to visualise the data 

# group the readings by pass number and then group them cycle wise for each pass

byPass <- df %>% group_by(passNo) %>% group_split() 
byCycle <- lapply(byPass, function(x){a <- x %>% group_by(cycleNo)
a<- group_split(a)})

# cycle with minimum number of data is chosen as reference cycle so that we find closer point pairs
ref_cyc <- lapply(byCycle, function(z){
  inds <- which.min(sapply(z, nrow))
  ref <- z[[inds]]; ref})

# rest points
other_cyc <- lapply(byCycle, function(z){
  inds <- which.max(sapply(z, nrow))
  other <- z[-inds];other
})

# Now lets find the pair of points close to each other

pairs <- NULL

for (i in 1:length(ref.cyc)){
  for (j in 1:length(other.cyc[[i]])){
    dist_mat <- distm(matrix(c(ref.cyc[[i]]$lon, ref.cyc[[i]]$lat), ncol=2), matrix(c(other.cyc[[i]][[j]]$lon, other.cyc[[i]][[j]]$lat), ncol=2), fun=distGeo)/1000
    # filter based on the dist_mat
    lat<- other.cyc[[i]][[j]]$lat[max.col(-dist_mat)]
    lon<- other.cyc[[i]][[j]]$lon[max.col(-dist_mat)]
    jdn<- other.cyc[[i]][[j]]$jdn[max.col(-dist_mat)]
    ssh<- other.cyc[[i]][[j]]$ssh[max.col(-dist_mat)]
    cycleNo <- other.cyc[[i]][[j]]$cycleNo[max.col(-dist_mat)]
    passNo <- other.cyc[[i]][[j]]$passNo[max.col(-dist_mat)]
    
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


# plotting point pairs for tenth pass 
visualise <- final.pairs[[10]]  # for checking 

ggplot()+geom_point(data = graph_test[["1"]], aes(lon, lat), color="red")+
  geom_point(data = visualise[["2"]], aes(lon, lat), color="blue")+
  geom_point(data = visualise[["3"]], aes(lon, lat), color="black")+
  geom_point(data = visualise[["4"]], aes(lon, lat), color="orange")+
  geom_point(data = visualise[["5"]], aes(lon, lat), color="purple")+
  geom_point(data = visualise[["6"]], aes(lon, lat), color="pink")+
  geom_point(data = visualise[["7"]], aes(lon, lat), color="cyan")+
  geom_point(data = visualise[["8"]], aes(lon, lat), color="yellow")+
  geom_point(data = visualise[["9"]], aes(lon, lat), color="green")



new_df <- list()

# length(final.pairs)

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
    
    # Least Square Collocation
    X <- inv(t(A)%*%A) %*% (t(A) %*% L)   # to find out the value of eta and zeta
    
    mean_lon <- mean(curr_grp$lon)
    mean_lat <- mean(curr_grp$lat)
    mean_ssh <- mean(curr_grp$ssh)
    eta <- X[1]
    zeta <- X[2]
    rbind(new_df, data.frame(lon=mean_lon, lat=mean_lat, ssh=mean_ssh, eta, zeta)) -> new_df
    
  }
}

# to visualise the data 
ggplot()+geom_point(data = new_df, aes(lon, lat), color="black")+scale_y_continuous(breaks=seq(14, 20, 0.5))+
  scale_x_continuous(breaks=seq(64, 70, 0.5))


# Code for kriging interpolation here 

# construction of spatial grid  
grd <- data.frame(lon = rep(seq(64, 70, 0.5), 13), lat = rep(seq(14, 20, 0.5), each = 13))
coordinates(grd) <- ~lon+lat
grd <- as.data.frame(spsample(grd, "regular", n = 144, offset = c(0.5, 0.5)))
names(grd) <- c("lon", "lat")
coordinates(grd) <- c("lon", "lat")
gridded(grd) <- TRUE  # Create SpatialPixel object
fullgrid(grd) <- TRUE # Create SpatialGrid object
proj4string(grd) <- "+init=epsg:4326"
# Add projection information to the grid

# variogram on the de - trended data.
coordinates(new_df) <- ~lon + lat # convert data frame to SPDF

# Add projection information to the processed.data
crs(new_df) <- "+init=epsg:4326"

var_zeta <- variogram(zeta~1, data = new_df, cloud = FALSE)

# to fit.variogram() via the vgm() function
dat_fit_zeta <- fit.variogram(var_zeta, vgm(c("Exp", "Mat", "Sph")))
# The following plot allows us to assess the fit
plot(var_zeta, dat_fit_zeta)

f1 <- as.formula(zeta~lon + lat)
zeta_krg <- krige(f1, new_df, grd, dat_fit_zeta)

# 
var_eta <- variogram(eta~1, data=new_df, cloud=FALSE)
dat_fit_eta <- fit.variogram(var_eta,  vgm(c("Exp", "Mat", "Sph")))
plot(var_eta, dat_fit_eta) # to visualise the fit 
f2 <- as.formula(eta~lon + lat)

eta_krg <- krige(f2, new_df, grd, dat_fit_eta)

interpolated_grd <- data.frame(grd, zeta_krg@data, eta_krg@data)

# interpolated_grd[,c("col num here")]<-NULL # some unwated inclusions

colnames(interpolated_grd) <- c("lon", "lat", "Pred.Zeta", "Zeta.var", "Pred.Eta", "Eta.var")

ggplot()+geom_point(data = interpolated_grd, aes(lon, lat), color="black")+scale_y_continuous(breaks=seq(14, 20, 0.5))+
  scale_x_continuous(breaks=seq(64, 70, 0.5))

# Now we need do the 2d fft to find out the gravity anomaly 


