library(Solar)

calcSun(datetime = Sys.time() + 8*3600)


dailySolar = function(x) {
  time = seq(x, x + 86400, by = 300)
  insolation = calcSun(datetime = time)$insolation
  mean(insolation, na.rm = T)
}

dailySolar(Sys.time())

insolation = data.frame(date = Sys.time() + c(1:(3600*60)), solar = 0)
insolation = calcSun(datetime = as.POSIXct(insolation$date))

plot(insolation$datetime, insolation$insolation,
     type ='l',
     yaxs = 'i',
     lwd = 3,
     ylab = 'Solar Radiation (rel)',
     xlab = '')

plot(insolation$datetime, insolation$solar.azimuth,
     type ='l',
     yaxs = 'i',
     lwd = 3,
     ylab = 'Solar Radiation (rel)',
     xlab = '')

plot(insolation$datetime, insolation$solar.elevation + insolation$refraction,
     type ='l',
     yaxs = 'i',
     lwd = 3,
     ylab = 'Solar Radiation (rel)',
     xlab = '')

a = ecdf(insolation$solar)
plot(insolation$solar/get.solar(Sys.time()), a(insolation$solar),
     type= 'l',
     lwd = 3, yaxs = 'i', xaxs = 'i', xlim = c(0, 3))

grid = expand.grid(lon = -147.8,
                   lat = 64.8,
                   datetime = seq(make.time(2022, 04, 01), make.time(2022,4,2), by = 60))


grid = calc.solar(grid$lon, grid$lat, tz = -10, datetime = grid$datetime)

sum(grid$insolation)

section = build.section(grid$solar.azimuth%%360, grid$solar.elevation, grid$insolation, gridder = gridBin, ylim = c(0, 90))
section.count = build.section(grid$solar.azimuth%%360, grid$solar.elevation, grid$insolation, gridder = gridCount, ylim = c(0, 90))
section.add = build.section(grid$solar.azimuth%%360, grid$solar.elevation, grid$insolation, gridder = gridAdd, ylim = c(0, 90))
section.time = build.section(grid$solar.azimuth%%360, as.numeric(grid$datetime), grid$insolation, gridder = gridAdd)


layer = expand.grid(lon = c(0:359), lat = c(-70:70))
layer = calc.solar(layer$lon, layer$lat, tz = -1, datetime = make.time(2022, 02, 28, 12,1,1, tz = 'US/Alaska'))

summary(layer)

map = make.map('coastlineWorldMedium', -170, -70, 30, 65)
map = make.map()
add.map.layer(map,
              lon = layer$lon,
              lat = layer$lat,
              z = layer$insolation,
              pal = 'inferno',
              trim = F)


add.radial.axis = function(p) {
  
  for (i in c(1e-3, 30, 60)) {
    add.map.line(list(p = p), lon = c(-180, 180), lat = rep(i, 2))
  }
  for (i in c(0, 90, 180, 270)) {
    add.map.line(list(p = p), lon = rep(i, 2), lat = c(0, 90))
  }
  mtext('0', side = 3, cex = 1)
  mtext('90', side = 4, cex = 1, las = 1)
  mtext('180', side = 1, cex = 1, las = 1)
  mtext('270', side = 2, cex = 1, las = 1)
}

p = '+proj=nsper +R=1 +lat_0=90 +lat_1=90 +h=1e6 + lon_0=180'
plot.default(NULL, NULL, xlim = c(-1,1), ylim = c(-1,1), axes = F, xlab = '', ylab = '')
add.map.layer(list(p = p), section$x, section$y, section$grid$z1, pal = 'parula', trim = F)
add.radial.axis(p)

plot.default(NULL, NULL, xlim = c(-1,1), ylim = c(-1,1), axes = F, xlab = '', ylab = '')
add.map.layer(list(p = p), section.count$x, section.count$y, section.count$grid$z1/365, pal = 'parula', trim = F)
add.radial.axis(p)

plot.default(NULL, NULL, xlim = c(-1,1), ylim = c(-1,1), axes = F, xlab = '', ylab = '')
add.map.layer(list(p = p), section.add$x, section.add$y, section.add$grid$z1/section.count$grid$z1, pal = 'cubicl', trim = F)
add.radial.axis(p)


plot.default(NULL, NULL, xlim = c(-1,1), ylim = c(-1,1), axes = F, xlab = '', ylab = '')
add.map.layer(list(p = p), section.add$x, section.add$y, section.add$grid$z1 * section.count$grid$z1, pal = 'inferno', trim = F)
add.radial.axis(p)




gridCount = function(gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {
  gz = rep(NA, length(gx))
  
  if (xscale == 0) { xscale = Inf }
  if (yscale == 0) { yscale = Inf }
  
  for (i in 1:length(gz)) {
    gz[i] = sum(abs(gx[i] - x) < xscale/2 & abs(gy[i] - y) < yscale/2)
  }
  
  gz ## Return
}

gridAdd = function(gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {
  gz = rep(NA, length(gx))
  
  if (xscale == 0) { xscale = Inf }
  if (yscale == 0) { yscale = Inf }
  
  for (i in 1:length(gz)) {
    gz[i] = sum(z[abs(gx[i] - x) < xscale/2 & abs(gy[i] - y) < yscale/2], na.rm = T)
  }
  
  gz ## Return
}

plot(grid$datetime, grid$insolation, type = 'l', ylim = c(0, 1))



day1 = calc.solar(datetime = seq(make.time(2022,02,27), make.time(2022,02,28), by = 60))
day2 = calc.solar(datetime = seq(make.time(2022,03,06), make.time(2022,03,07), by = 60))

plot(get.hour(day1$datetime) + get.minutes(day1$datetime)/60, day1$insolation, type = 'l')
lines(x = get.hour(day2$datetime) + get.minutes(day2$datetime)/60, y = day2$insolation, col = 'blue')
