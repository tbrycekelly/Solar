library(Solar)

calcSun(datetime = Sys.time() + 8*3600)


dailySolar = function(x) {
  time = seq(x, x + 86400, by = 300)
  insolation = calcSun(datetime = time, lon = -72.3, lat = 42.7)$insolation
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



grid = grid[grid$insolation > 0.1, ]

sum(grid$insolation)

library(SimpleGridder)


gridCount = function(tree, z, gx, gy, neighborhood = 25, weight.func = function(x) {x < 1}) {
  
  grid = data.frame(x = gx, y = gy)
  tmp = tree$query(grid, neighborhood)
  
  grid$z = NA
  
  for (i in 1:nrow(grid)) {
    w = weight.func(tmp$nn.dist[i,])
    grid$z[i] = sum(w)
  }
  
  grid
}



grid = expand.grid(lon = -72.3,
                   lat = 42.7,
                   datetime = seq(as.POSIXct('2022-01-01'), as.POSIXct('2023-12-31'), by = 60*10))


grid = calcSun(lon = grid$lon, lat = grid$lat, datetime = grid$datetime)

grid$insolation = grid$insolation * cos(pi / 180 * (30 - grid$solar.elevation)) * cos(pi / 180 * (180 - grid$solar.azimuth))
grid = grid[grid$insolation > 0.1, ]

#### Calculate location in the sky (i.e. count)
product = SimpleGridder::buildGrid(xlim = c(0, 360), ylim = c(0, 90), nx = 60, ny = 20, x.factor = 1, y.factor = 1)
product = SimpleGridder::setGridder(product, gridder = gridCount, neighborhood = 400)
product = SimpleGridder::appendData(product,
                                    x = grid$solar.azimuth,
                                    y = grid$solar.elevation,
                                    z = grid$insolation,
                                    label = 'insolation')
product = SimpleGridder::interpData(product)


zlim = c(0, max(pretty(product$interp$insolation)))
zlim = c(0, 5000)
SimpleGridder::plotGrid(product, 'insolation', pal = pals::inferno(8), zlim = zlim)
mtext(paste(zlim, collapse = ' -> '))


#### Calculate insolation in the sky (i.e. count)
light = SimpleGridder::buildGrid(xlim = c(0, 360), ylim = c(0, 90), nx = 60, ny = 20, x.factor = 1, y.factor = 1)
light = SimpleGridder::setGridder(light, gridder = gridWeighted, neighborhood = 100)
light = SimpleGridder::appendData(light,
                                    x = grid$solar.azimuth,
                                    y = grid$solar.elevation,
                                    z = grid$insolation,
                                    label = 'insolation')
light = SimpleGridder::interpData(light)


zlim = c(0, max(pretty(light$interp$insolation)))
SimpleGridder::plotGrid(light, 'insolation', pal = pals::inferno(8), zlim = zlim)
mtext(paste(zlim, collapse = ' -> '))

light$interp$light = light$interp$insolation * product$interp$insolation
zlim = c(0, max(pretty(light$interp$light)))
zlim = c(0, 250 * 60/2)
SimpleGridder::plotGrid(light, 'light', pal = pals::inferno(8), zlim = zlim)
mtext(paste(zlim, collapse = ' -> '))







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
