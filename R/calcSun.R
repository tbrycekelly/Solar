
#' @export
calcSun = function(lon = -147.8,
                   lat = 64.8,
                   datetime = Sys.time()
                   ) {
  
  if (length(lon) == 1) {
    lon = rep(lon, length(datetime))
  }
  if (length(lat) == 1) {
    lat = rep(lat, length(datetime))
  }
  
  jules = julian(datetime)
  
  geom.mean.long.sun = (280.46646 + jules$century * (36000.76983 + jules$century * 0.0003032)) %% 360
  geom.mean.anom.sun = 357.52911 + jules$century * (35999.05029 - 0.0001537 * jules$century)
  
  eccentric.earth = 0.016708634 - jules$century * (0.000042037 + 0.0000001267 * jules$century)
  sun.center = sin(geom.mean.anom.sun * pi / 180) * (1.914602 - jules$century * (0.004817 + 0.000014 * jules$century)) +
    sin(2 * geom.mean.anom.sun * pi / 180) * (0.019993 - 0.000101 * jules$century) +
    sin(3 * geom.mean.anom.sun * pi / 180) * (0.000289)
  
  sun.true.lon = geom.mean.long.sun + sun.center
  sun.true.anom = geom.mean.anom.sun + sun.center
  sun.rad.vector = 1.000001018 * (1 - eccentric.earth^2) / (1 + eccentric.earth * cos(sun.true.anom * pi /180))
  sun.apparent.lon = sun.true.lon - 0.00569 - 0.00478 * sin((125.04 - 1934.136 * jules$century) * pi / 180)
  
  mean.obl.ecliptic = 23 + (26 + ((21.448 - jules$century * (46.815 + jules$century * (0.00059 - jules$century * 0.001813))))/60)/60
  obl.correction = mean.obl.ecliptic + 0.00256 * cos((125.04 - 1934.136 * jules$century) * pi / 180)
  
  sun.right.ass = 180 / pi * atan2(cos(obl.correction * pi / 180) * sin(sun.apparent.lon * pi / 180), cos(sun.apparent.lon * pi / 180))
  sun.dec = 180 / pi * asin(sin(pi / 180 * obl.correction) * sin(pi / 180 * sun.apparent.lon))
  
  y = tan(obl.correction / 2 * pi / 180)^2
  
  eq.time.min = y * sin(2 * geom.mean.long.sun * pi / 180) -
                                   2 * eccentric.earth * sin(geom.mean.anom.sun * pi / 180) +
                                   4 * eccentric.earth * y * sin(geom.mean.anom.sun * pi / 180) * cos(2 * geom.mean.long.sun * pi / 180) -
                                   0.5 * y^2 * sin(4 * geom.mean.long.sun * pi / 180) -
                                   1.25 * eccentric.earth^2 * sin(2 * geom.mean.anom.sun * pi / 180)
  eq.time.min = 4 * 180 / pi * eq.time.min
  
  tmp = cos(90.833 * pi / 180) / cos(lat * pi / 180) / cos(pi / 180 * sun.dec) - tan(lat * pi / 180) * tan(sun.dec * pi / 180)
  HA.sunrise = 180 / pi * acos(pmin(1, pmax(-1, tmp)))
  
  #solar.noon = (720 - 4 * lon - eq.time.min + tz * 60) / 1440 # local
  solar.noon = (720 - 4 * lon - eq.time.min) / 1440 # UTC
  sunrise = solar.noon - HA.sunrise * 4 / 1440
  sunset = solar.noon + HA.sunrise * 4 / 1440
  
  sunlight.durration = 8 * HA.sunrise
  
  true.solar.time = ((as.numeric(format(datetime, '%H')) + as.numeric(format(datetime, '%M'))/60)/24 * 1440 + eq.time.min + 4 * lon - tz * 60) %% 1440
  
  hour.angle = rep(0, length(datetime))
  
  k = true.solar.time < 0
  hour.angle[k] = 180 + true.solar.time[k]/4
  hour.angle[!k] = true.solar.time[!k]/4 - 180
  
  tmp = sin(pi / 180 * lat) * sin(pi / 180 * sun.dec) + cos(lat * pi / 180) * cos(sun.dec * pi / 180) * cos(hour.angle * pi / 180)
  solar.zenith = 180 / pi * acos(pmin(1,pmax(-1,tmp)))
  solar.elevation = 90 - solar.zenith
  
  ## Calculate atmospheric refraction and attenuation
  refraction = rep(0, length(datetime))
  
  k = (solar.elevation <=85 & solar.elevation > 5)
  refraction[k] = 58.1 / tan(solar.elevation[k] * pi / 180) - 0.07 / tan(solar.elevation[k] * pi / 180)^3 + 0.000086 / tan(solar.elevation[k] * pi / 180)^5
  
  k = (solar.elevation <= 5 & solar.elevation > -0.575)
  refraction[k] = 1735 + solar.elevation[k] * (-518.2 + solar.elevation[k] * (103.4 + solar.elevation[k] * (-12.79 + solar.elevation[k] * 0.711)))
  
  k = (solar.elevation <= -.0575)
  refraction[k] = -20.772 / tan(solar.elevation[k] * pi / 180)
  
  refraction = refraction / 3600
  
  solar.elevation.true = solar.elevation + refraction
  
  solar.azimuth = rep(0, length(datetime))
  k = hour.angle > 0
  tmp = pmax(-1,pmin(1,(sin(lat[k] * pi / 180) * cos(solar.zenith[k] * pi / 180) - sin(sun.dec[k] * pi / 180)) / (cos(lat[k] * pi / 180) * sin(solar.zenith[k] * pi / 180))))
  solar.azimuth[k] = (180 / pi * acos(tmp)) + 180 %% 360
  tmp = pmin(1, pmax(-1,(sin(lat[!k] * pi / 180) * cos(solar.zenith[!k] * pi / 180) - sin(sun.dec[!k] * pi / 180)) / (cos(lat[!k] * pi / 180) * sin(solar.zenith[!k] * pi / 180))))
  solar.azimuth[!k] = (540 - 180 / pi * acos(tmp)) %% 360
  
  AM = 1 / cos((90 - solar.elevation.true) * pi / 180)
  AM[is.na(AM)] = 0
  ID = 1.353 * 0.7 ^ (AM^0.678)
  ID[is.na(ID)] = 0
  
  data.frame(lonitude = lon,
             latitude = lat,
             datetime = format(datetime, usetz = T),
             julesDay = jules$day,
             julestCentury = jules$century,
             solarNoon = solar.noon,
             sunrise = sunrise,
             sunset = sunset,
             solarHour = hour.angle,
             solarZenithDeg = solar.zenith,
             solarAzimuthDeg = solar.azimuth,
             solarElevationDeg = solar.elevation,
             solarElevationApparentDeg = solar.elevation.true,
             refractionDeg = refraction,
             insolation = ID)
}