## solar calculator


calc.solar = function(lon = -147.8,
                      lat = 64.8,
                      tz = -10,
                      datetime = Sys.time()) {

  pi = 3.1415926535897

  julian.day = conv.time.excel(datetime, rev = T) - tz / 24 + 2415018.5
  julian.century = (julian.day - 2451545) / 36525

  geom.mean.long.sun = (280.46646 + julian.century * (36000.76983 + julian.century * 0.0003032)) %% 360
  geom.mean.anom.sun = 357.52911 + julian.century * (35999.05029 - 0.0001537 * julian.century)

  eccentric.earth = 0.016708634 - julian.century * (0.000042037 + 0.0000001267 * julian.century)
  sun.center = sin(geom.mean.anom.sun * pi / 180) * (1.914602 - julian.century * (0.004817 + 0.000014 * julian.century)) +
    sin(2 * geom.mean.anom.sun * pi / 180) * (0.019993 - 0.000101 * julian.century) +
    sin(3 * geom.mean.anom.sun * pi / 180) * (0.000289)

  sun.true.lon = geom.mean.long.sun + sun.center
  sun.true.anom = geom.mean.anom.sun + sun.center
  sun.rad.vector = 1.000001018 * (1 - eccentric.earth^2) / (1 + eccentric.earth * cos(sun.true.anom * pi /180))
  sun.apparent.lon = sun.true.lon - 0.00569 - 0.00478 * sin((125.04 - 1934.136 * julian.century) * pi / 180)

  mean.obl.ecliptic = 23 + (26 + ((21.448 - julian.century * (46.815 + julian.century * (0.00059 - julian.century * 0.001813))))/60)/60
  obl.correction = mean.obl.ecliptic + 0.00256 * cos((125.04 - 1934.136 * julian.century) * pi / 180)

  sun.right.ass = 180 / pi * atan2(cos(obl.correction * pi / 180) * sin(sun.apparent.lon * pi / 180), cos(sun.apparent.lon * pi / 180))
  sun.dec = 180 / pi * asin(sin(pi / 180 * obl.correction) * sin(pi / 180 * sun.apparent.lon))

  y = tan(obl.correction / 2 * pi / 180)^2

  eq.time.min = 4 * 180 / pi * ( y * sin(2 * geom.mean.long.sun * pi / 180) -
                                   2 * eccentric.earth * sin(geom.mean.anom.sun * pi / 180) +
                                   4 * eccentric.earth * y * sin(geom.mean.anom.sun * pi / 180) * cos(2 * geom.mean.long.sun * pi / 180) -
                                   0.5 * y^2 * sin(4 * geom.mean.long.sun * pi / 180) -
                                   1.25 * eccentric.earth^2 * sin(geom.mean.anom.sun * pi / 180))
  HA.sunrise = 180 / pi * acos(cos(90.833 * pi / 180) / cos(lat * pi / 180) / cos(pi / 180 * sun.dec) - tan(lat * pi / 180) * tan(sun.dec * pi / 180))

  solar.noon = (720 - 4 * lon - eq.time.min + tz * 60) / 1440
  sunrise = solar.noon - HA.sunrise * 4 / 1440
  sunset = solar.noon + HA.sunrise * 4 / 1440

  sunlight.durration = 8 * HA.sunrise

  true.solar.time = ((get.julian(datetime) - floor(get.julian(datetime))) * 1440 + eq.time.min + 4 * lon - 60 * tz) %% 1440

  hour.angle = rep(0, length(datetime))

  if (true.solar.time < 0) {
    hour.angle = 180 + true.solar.time/4
  } else {
    hour.angle = true.solar.time/4 - 180
  }

  solar.zenith = 180 / pi * acos(sin(pi / 180 * lat) * sin(pi / 180 * sun.dec) + cos(lat * pi / 180) * cos(sun.dec * pi / 180) * cos(hour.angle * pi / 180))
  solar.elevation = 90 - solar.zenith

  ## Calculate atmospheric refraction and attenuation
  refraction = 0
  if (solar.elevation <=85 & solar.elevation > 5) {
    refraction = 58.1 / tan(solar.elevation * pi / 180) - 0.07 / tan(solar.elevation * pi / 180)^3 + 0.000086 / tan(solar.elevation * pi / 180)^5
  } else if (solar.elevation <= 5 & solar.elevation > -0.575) {
    refraction = 1735 + solar.elevation * (-518.2 + solar.elevation * (103.4 + solar.elevation * (-12.79 + solar.elevation * 0.711)))
  } else if (solar.elevation <= -.0575) {
    refraction = -20.772 / tan(solar.elevation * pi / 180)
  }
  refraction = refraction / 3600

  solar.elevation.true = solar.elevation + refraction

  solar.azimuth = rep(0, length(datetime))
  l = hour.angle > 0
  if (hour.angle > 0) {
    solar.azimuth = (180 / pi * acos((sin(lat * pi / 180) * cos(solar.zenith * pi / 180) - sin(sun.dec * pi / 180)) / (cos(lat * pi / 180) * sin(solar.zenith * pi / 180)))) + 180 %% 360
  } else {
    solar.azimuth = (540 - 180 / pi * acos((sin(lat * pi / 180) * cos(solar.zenith * pi / 180) - sin(sun.dec * pi / 180)) / (cos(lat * pi / 180) * sin(solar.zenith * pi / 180)))) %% 360
  }

  AM = 1 / cos(pi / 180 * (90 - solar.elevation.true))
  AM[is.na(AM)] = 0
  ID = rep(0, length(datetime))
  l = AM > 0
  ID = 1.353 * 0.7 ^ (AM^0.678)
  ID[is.na(ID)] = 0

  data.frame(lon = lon,
             lat = lat,
             datetime = datetime,
             tz = tz,
             solar.elevation = solar.elevation,
             solar.dec = sun.dec,
             solar.hour = hour.angle,
             solar.zenith = solar.zenith,
             solar.azimuth = solar.azimuth,
             insolation = ID)
}


get.solar = function(x) {
  time = seq(x, x + 86400, by = 300)
  insolation = rep(0, length(time))
  for (i in 1:length(time)) {
    insolation[i] = calc.solar(datetime = time[i])$insolation
  }
  mean(insolation, na.rm = T)
}

get.solar(Sys.time())

insolation = data.frame(date = make.time(2022) + c(1:365) * 86400, solar = 0)
for (i in 1:nrow(insolation)) {
  insolation$solar[i] = get.solar(insolation$date[i])
}

plot(insolation$date, insolation$solar/ get.solar(Sys.time()),
     type ='l',
     ylim = c(0,3),
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
