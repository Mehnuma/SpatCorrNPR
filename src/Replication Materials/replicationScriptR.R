
library(gstat)
library(geostats)
library(npsp)
library(sp)
library(sf)
library(tictoc)

#************************** Dataset: Texas Wildfire********************#
rmse_wildfire_gstat = NULL
runtime2_wildfire_gstat = NULL
ind_wildfire1=1
ind_wildfire2=1
n1_wildfire = 178
n2_wildfire = 76
for(i in 1:10){
  tic()
  train_wildfire = read.csv('Train_wildfire.csv')
  train_wildfire = train_wildfire[ind_wildfire1:(n1_wildfire*i),]
  test_wildfire = read.csv('Test_wildfire.csv')
  test_wildfire = test_wildfire[ind_wildfire2:(n2_wildfire*i),]
  train_wildfire_x = train_wildfire$lat
  train_wildfire_y = train_wildfire$lon
  train_wildfire_z = train_wildfire$exp_area
  train_wildfire_z = (train_wildfire_z-min(train_wildfire_z))/(max(train_wildfire_z)-min(train_wildfire_z))
  train_wildfire_sf <- st_as_sf(train_wildfire, coords = c("lat", "lon"))
  test_wildfire_x = test_wildfire$lat
  test_wildfire_y = test_wildfire$lon
  test_wildfire_z = test_wildfire$exp_area
  test_wildfire_z = (test_wildfire_z-min(test_wildfire_z))/(max(test_wildfire_z)-min(test_wildfire_z))
  test_wildfire_sf <- st_as_sf(test_wildfire, coords = c("lat", "lon"))
  v <- gstat::variogram(train_wildfire_z~1,train_wildfire_sf)
  vm <- gstat::fit.variogram(v, fit.sills = F, fit.ranges = F, model = vgm(model = c("Sph", "Exp", "Gau")))
  gstat_kriging = gstat::krige(train_wildfire_z~1,
                               train_wildfire_sf, newdata = test_wildfire_sf, model = vm)
  yhat_wildfire_gstat <- gstat_kriging$var1.pred
  ind_wildfire1 = n1_wildfire*i+1;
  ind_wildfire2 = n2_wildfire*i+1;
  rmse_wildfire_gstat[i] = sqrt(mean((test_wildfire_z -yhat_wildfire_gstat)^2))
  end_time = toc(quiet=TRUE)
  runtime2_wildfire_gstat[i] = end_time$toc-end_time$tic
}

rmse_wildfire_geostats = NULL
runtime2_wildfire_geostats = NULL
ind_wildfire1=1
ind_wildfire2=1
n1_wildfire = 178
n2_wildfire = 76
for(i in 1:10){
  tic()
  train_wildfire = read.csv('Train_wildfire.csv')
  train_wildfire = train_wildfire[ind_wildfire1:(n1_wildfire*i),]
  test_wildfire = read.csv('Test_wildfire.csv')
  test_wildfire = test_wildfire[ind_wildfire2:(n2_wildfire*i),]
  train_wildfire_x = train_wildfire$lat
  train_wildfire_y = train_wildfire$lon
  train_wildfire_z = train_wildfire$exp_area
  train_wildfire_z = (train_wildfire_z-min(train_wildfire_z))/(max(train_wildfire_z)-min(train_wildfire_z))
  train_wildfire_sf <- st_as_sf(train_wildfire, coords = c("lat", "lon"))
  test_wildfire_x = test_wildfire$lat
  test_wildfire_y = test_wildfire$lon
  test_wildfire_z = test_wildfire$exp_area
  test_wildfire_z = (test_wildfire_z-min(test_wildfire_z))/(max(test_wildfire_z)-min(test_wildfire_z))
  test_wildfire_sf <- st_as_sf(test_wildfire, coords = c("lat", "lon"))
  v <- geostats::semivariogram(train_wildfire_x, train_wildfire_y, train_wildfire_z, model="spherical")
  yhat_wildfire_geostats = geostats::kriging(train_wildfire_x,train_wildfire_y, train_wildfire_z,
                                       test_wildfire_x, test_wildfire_y, svm = v)
  ind_wildfire1 = n1_wildfire*i+1;
  ind_wildfire2 = n2_wildfire*i+1;
  rmse_wildfire_geostats[i] = sqrt(mean((test_wildfire_z -yhat_wildfire_geostats)^2))
  end_time = toc(quiet=TRUE)
  runtime2_wildfire_geostats[i] = end_time$toc-end_time$tic
}

rmse_wildfire_npsp = NULL
runtime2_wildfire_npsp = NULL
ind_wildfire1=1
ind_wildfire2=1
n1_wildfire = 178
n2_wildfire = 76
for(i in 1:10){
  tic()
  train_wildfire = read.csv('Train_wildfire.csv')
  train_wildfire = train_wildfire[ind_wildfire1:(n1_wildfire*i),]
  test_wildfire = read.csv('Test_wildfire.csv')
  test_wildfire = test_wildfire[ind_wildfire2:(n2_wildfire*i),]
  train_wildfire_x = train_wildfire$lat
  train_wildfire_y = train_wildfire$lon
  train_wildfire_z = train_wildfire$exp_area
  train_wildfire_z = (train_wildfire_z-min(train_wildfire_z))/(max(train_wildfire_z)-min(train_wildfire_z))
  test_wildfire_x = test_wildfire$lat
  test_wildfire_y = test_wildfire$lon
  test_wildfire_z = test_wildfire$exp_area
  test_wildfire_z = (test_wildfire_z-min(test_wildfire_z))/(max(test_wildfire_z)-min(test_wildfire_z))
  lp_wildfire = locpol(train_wildfire[,1:2], train_wildfire$exp_area, h = diag(1.8,2), hat.bin = TRUE)
  esvar = np.svariso.corr(lp_wildfire, h = 3, plot = FALSE, max.iter = 100)
  svar_wildfire <- fitsvar.sb.iso(esvar)
  yhat_wildfire_npsp <- kriging.simple(train_wildfire[,1:2], train_wildfire$exp_area, test_wildfire[,1:2], svar_wildfire)
  ind_wildfire1 = n1_wildfire*i+1;
  ind_wildfire2 = n2_wildfire*i+1;
  rmse_wildfire_npsp[i] = sqrt(mean((test_wildfire_z -yhat_wildfire_npsp$kpred)^2))
  end_time = toc(quiet=TRUE)
  runtime2_wildfire_npsp[i] = end_time$toc-end_time$tic
}

#************************** Dataset: US Precipitation ********************#
rmse_precip_gstat = NULL
runtime2_precip_gstat = NULL
ind_precip1=1
ind_precip2=1
n1_precip = 738
n2_precip = 315
for(i in 1:10){
  tic()
  train_precip = read.csv('Train_precipitation.csv')
  train_precip = train_precip[ind_precip1:(n1_precip*i),]
  test_precip = read.csv('Test_precipitation.csv')
  test_precip = test_precip[ind_precip2:(n2_precip*i),]
  train_precip_x = train_precip$lat
  train_precip_y = train_precip$lon
  train_precip_z = train_precip$precipitation
  train_precip_z = (train_precip_z-min(train_precip_z))/(max(train_precip_z)-min(train_precip_z))
  train_precip_sf <- st_as_sf(train_precip, coords = c("lat", "lon"))
  test_precip_x = test_precip$lat
  test_precip_y = test_precip$lon
  test_precip_z = test_precip$precipitation
  test_precip_z = (test_precip_z-min(test_precip_z))/(max(test_precip_z)-min(test_precip_z))
  test_precip_sf <- st_as_sf(test_precip, coords = c("lat", "lon"))
  v <- gstat::variogram(train_precip_z~1,train_precip_sf)
  vm <- gstat::fit.variogram(v, model = vgm(c("Sph", "Exp", "Gau")))
  gstat_kriging = gstat::krige(train_precip_z~1,
                               train_precip_sf, newdata = test_precip_sf, model = vm)
  yhat_precip_gstat <- gstat_kriging$var1.pred
  ind_precip1 = n1_precip*i+1;
  ind_precip2 = n2_precip*i+1;
  rmse_precip_gstat[i] = sqrt(mean((test_precip_z -yhat_precip_gstat)^2))
  end_time = toc(quiet=TRUE)
  runtime2_precip_gstat[i] = end_time$toc-end_time$tic
}

rmse_precip_geostats = NULL
runtime2_precip_geostats = NULL
ind_precip1=1
ind_precip2=1
n1_precip = 738
n2_precip = 315
for(i in 1:10){
  tic()
  train_precip = read.csv('Train_precipitation.csv')
  train_precip = train_precip[ind_precip1:(n1_precip*i),]
  test_precip = read.csv('Test_precipitation.csv')
  test_precip = test_precip[ind_precip2:(n2_precip*i),]
  train_precip_x = train_precip$lat
  train_precip_y = train_precip$lon
  train_precip_z = train_precip$precipitation
  train_precip_z = (train_precip_z-min(train_precip_z))/(max(train_precip_z)-min(train_precip_z))
  test_precip_x = test_precip$lat
  test_precip_y = test_precip$lon
  test_precip_z = test_precip$precipitation
  test_precip_z = (test_precip_z-min(test_precip_z))/(max(test_precip_z)-min(test_precip_z))
  v <- geostats::semivariogram(train_precip_x, train_precip_y, train_precip_z, model="spherical")
  yhat_precip_geostats = geostats::kriging(train_precip_x,train_precip_y, train_precip_z,
                                       test_precip_x, test_precip_y, svm = v)
  ind_precip1 = n1_precip*i+1;
  ind_precip2 = n2_precip*i+1;
  rmse_precip_geostats[i] = sqrt(mean((test_precip_z -yhat_precip_geostats)^2))
  end_time = toc(quiet=TRUE)
  runtime2_precip_geostats[i] = end_time$toc-end_time$tic
}

rmse_precip_npsp = NULL
runtime2_precip_npsp = NULL
ind_precip1=1
ind_precip2=1
n1_precip = 738
n2_precip = 315
for(i in 1:10){
  tic()
  train_precip = read.csv('Train_precipitation.csv')
  train_precip = train_precip[ind_precip1:(n1_precip*i),]
  test_precip = read.csv('Test_precipitation.csv')
  test_precip = test_precip[ind_precip2:(n2_precip*i),]
  train_precip_x = train_precip$lat
  train_precip_y = train_precip$lon
  train_precip_z = train_precip$precipitation
  test_precip_x = test_precip$lat
  test_precip_y = test_precip$lon
  test_precip_z = test_precip$precipitation
  lp_precip = locpol(train_precip[,1:2], train_precip$precipitation, h = diag(60,2), hat.bin = TRUE)
  x = coordinates(train_precip[1:2])
  svar.bin <- svariso(x, residuals(lp_precip))
  svar.h <- h.cv(svar.bin)$h
  esvar = np.svariso.corr(lp_precip, h = svar.h, plot = FALSE, max.iter = 100)
  svar_precip <- fitsvar.sb.iso(esvar)
  yhat_precip_npsp <- kriging.simple(train_precip[,1:2], train_precip$precipitation, test_precip[,1:2],svar_precip)
  ind_precip1 = n1_precip*i+1;
  ind_precip2 = n2_precip*i+1;
  rmse_precip_npsp[i] = sqrt(mean((test_precip_z -yhat_precip_npsp$kpred)^2))
  end_time = toc(quiet=TRUE)
  runtime2_precip_npsp[i] = end_time$toc-end_time$tic
}

#************************** Dataset: Global Earthquake ********************#
rmse_eq_gstat = NULL
runtime2_eq_gstat = NULL
ind_eq1=1
ind_eq2=1
n1_eq = 2309
n2_eq = 989
for(i in 1:10){
  tic()
  train_eq = read.csv('Train_earthquake.csv')
  train_eq = train_eq[ind_eq1:(n1_eq*i),]
  test_eq = read.csv('Test_earthquake.csv')
  test_eq = test_eq[ind_eq2:(n2_eq*i),]
  train_eq_x = train_eq$lat
  train_eq_y = train_eq$lon
  train_eq_z = train_eq$mag
  train_eq_z = (train_eq_z-min(train_eq_z))/(max(train_eq_z)-min(train_eq_z))
  train_eq_sf <- st_as_sf(train_eq, coords = c("lat", "lon"))
  test_eq_x = test_eq$lat
  test_eq_y = test_eq$lon
  test_eq_z = test_eq$mag
  test_eq_z = (test_eq_z-min(test_eq_z))/(max(test_eq_z)-min(test_eq_z))
  test_eq_sf <- st_as_sf(test_eq, coords = c("lat", "lon"))
  v <- gstat::variogram(train_eq_z~1, train_eq_sf)
  vm <- gstat::fit.variogram(v, model = vgm(c("Sph", "Exp", "Gau")))
  gstat_kriging = gstat::krige(train_eq_z~1,
                               train_eq_sf, newdata = test_eq_sf, model = vm)
  yhat_eq_gstat <- gstat_kriging$var1.pred
  ind_eq1 = n1_eq*i+1;
  ind_eq2 = n2_eq*i+1;
  rmse_eq_gstat[i] = sqrt(mean((test_eq_z -yhat_eq_gstat)^2))
  end_time = toc(quiet=TRUE)
  runtime2_eq_gstat[i] = end_time$toc-end_time$tic
}

## The following code takes over four days to run, and therefore we have 
#included the rmse and runtime from the paper. Feel free to run it by uncommenting.

# rmse_eq_geostats = NULL
# runtime2_eq_geostats = NULL
# ind_eq1=1
# ind_eq2=1
# n1_eq = 2309
# n2_eq = 989
# for(i in 1:10){
#   tic()
#   train_eq = read.csv('Train_earthquake.csv')
#   train_eq = train_eq[ind_eq1:(n1_eq*i),]
#   test_eq = read.csv('Test_earthquake.csv')
#   test_eq = test_eq[ind_eq2:(n2_eq*i),]
# 
#   train_eq_x = train_eq$lat
#   train_eq_y = train_eq$lon
#   train_eq_z = train_eq$mag
#   train_eq_z = (train_eq_z-min(train_eq_z))/(max(train_eq_z)-min(train_eq_z))
# 
#   test_eq_x = test_eq$lat
#   test_eq_y = test_eq$lon
#   test_eq_z = test_eq$mag
#   test_eq_z = (test_eq_z-min(test_eq_z))/(max(test_eq_z)-min(test_eq_z))
# 
#   v <- geostats::semivariogram(train_eq_x, train_eq_y, train_eq_z, model="spherical")
#   yhat_eq_geostats = geostats::kriging(train_eq_x,train_eq_y, train_eq_z,
#                                         test_eq_x, test_eq_y, svm = v)
# 
#   ind_eq1 = n1_eq*i+1;
#   ind_eq2 = n2_eq*i+1;
#   rmse_eq_geostats[i] = sqrt(mean((test_eq_z -yhat_eq_geostats)^2))
#   end_time = toc(quiet=TRUE)
#   runtime2_eq_geostats[i] = end_time$toc-end_time$tic
# }
mean_rmse_eq_geostats = 0.135
sd_rmse_eq_geostats = 0.017
mean_runtime2_eq_geostats = 333268.15
sd_runtime2_eq_geostats = 1173.39

rmse_eq_npsp = NULL
runtime2_eq_npsp = NULL
ind_eq1=1
ind_eq2=1
n1_eq = 2309
n2_eq = 989
for(i in 1:10){
  tic()
  train_eq = read.csv('Train_earthquake.csv')
  train_eq = train_eq[ind_eq1:(n1_eq*i),]
  test_eq = read.csv('Test_earthquake.csv')
  test_eq = test_eq[ind_eq2:(n2_eq*i),]
  train_eq_x = train_eq$lat
  train_eq_y = train_eq$lon
  train_eq_z = train_eq$mag
  test_eq_x = test_eq$lat
  test_eq_y = test_eq$lon
  test_eq_z = test_eq$mag
  x = coordinates(train_eq[1:2])
  bin <- binning(x, train_eq_z)
  lp0.h = diag(30,2)
  lp0 <- locpol(bin, h = lp0.h, hat.bin = TRUE)
  svar.bin <- svariso(x, residuals(lp0), nlags = 60, maxlag = 70)
  svar.h <- h.cv(svar.bin)$h
  esvar = np.svariso.corr(lp0, h = svar.h, plot = FALSE, max.iter = 100)
  svar_eq <- fitsvar.sb.iso(esvar)
  yhat_eq_npsp <- kriging.simple(train_eq[,1:2], train_eq$mag, test_eq[,1:2], svar_eq)
  ind_eq1 = n1_eq*i+1;
  ind_eq2 = n2_eq*i+1;
  rmse_eq_npsp[i] = sqrt(mean((test_eq_z -yhat_eq_npsp$kpred)^2))
  end_time = toc(quiet=TRUE)
  runtime2_eq_npsp[i] = end_time$toc-end_time$tic
}

### Print Results
print(paste0("%%%%%%%%%%%%%%%%%%%%%% Table 7 Results %%%%%%%%%%%%%%%%%%%%%%"))
print(paste0("RMSE for Dataset-1 (Column 1 of Table 7): "))
print(paste0("Dataset-1, gstat RMSE: ",mean(rmse_wildfire_gstat),"(",sd(rmse_wildfire_gstat),")\n"))
print(paste0("Dataset-1, geostats RMSE: ",mean(rmse_wildfire_geostats),"(",sd(rmse_wildfire_geostats),")"))
print(paste0("Dataset-1, npsp RMSE: ",mean(rmse_wildfire_npsp),"(",sd(rmse_wildfire_npsp),")"))

print(paste0("RMSE for Dataset-2 (Column 2 of Table 7): "))
print(paste0("Dataset-2, gstat RMSE: ",mean(rmse_precip_gstat),"(",sd(rmse_precip_gstat),")"))
print(paste0("Dataset-2, geostats RMSE: ",mean(rmse_precip_geostats),"(",sd(rmse_precip_geostats),")"))
print(paste0("Dataset-2, npsp RMSE: ",mean(rmse_precip_npsp),"(",sd(rmse_precip_npsp),")"))

print(paste0("RMSE for Dataset-3 (Column 3 of Table 7): "))
print(paste0("Dataset-3, gstat runtime: ",mean(runtime2_eq_gstat),"(",sd(runtime2_eq_gstat),")"))
print(paste0("Dataset-3, geostats RMSE: ",mean_rmse_eq_geostats,"(",sd_rmse_eq_geostats,")"))
print(paste0("Dataset-3, npsp RMSE: ",mean(rmse_eq_npsp),"(",sd(rmse_eq_npsp),")"))
print(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"))

print(paste0("%%%%%%%%%%%%%%%%%%%%%% Table 8 Results %%%%%%%%%%%%%%%%%%%%%%"))
print(paste0("Runtimes for Dataset-1 (Column 1 of Table 8): "))
print(paste0("Dataset-1, gstat runtime: ",mean(runtime2_wildfire_gstat),"(",sd(runtime2_wildfire_gstat),")"))
print(paste0("Dataset-1, geostats runtime: ",mean(runtime2_wildfire_geostats),"(",sd(runtime2_wildfire_geostats),")"))
print(paste0("Dataset-1, npsp runtime: ",mean(runtime2_wildfire_npsp),"(",sd(runtime2_wildfire_npsp),")"))

print(paste0("Runtimes for Dataset-1 (Column 2 of Table 8): "))
print(paste0("Dataset-2, gstat runtime: ",mean(runtime2_precip_gstat),"(",sd(runtime2_precip_gstat),")"))
print(paste0("Dataset-2, geostats runtime: ",mean(runtime2_precip_geostats),"(",sd(runtime2_precip_geostats),")"))
print(paste0("Dataset-2, npsp runtime: ",mean(runtime2_precip_npsp),"(",sd(runtime2_precip_npsp),")"))

print(paste0("Runtimes for Dataset-1 (Column 3 of Table 8): "))
print(paste0("Dataset-3, gstat runtime: ",mean(runtime2_eq_gstat),"(",sd(runtime2_eq_gstat),")"))
print(paste0("Dataset-3, geostats runtime: ",mean_runtime2_eq_geostats,"(",sd_runtime2_eq_geostats,")"))
print(paste0("Dataset-3, npsp runtime: ",mean(runtime2_eq_npsp),"(",sd(runtime2_eq_npsp),")"))
print(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"))
