
sun.mean = function(datetime) {
  century = julian.century(datetime)
  mean.lon = 280.46646 + century * (36000.76983 + century * 0.0003032) %% 360
  mean.anom = 357.52911 + century * (35999.05029 - 0.0001537 * century) # Radians
  eccentric = 0.016708634 - century * (0.000042037 + 0.0000001267 * century)
  sun.center = sin(mean.anom) * (1.914602 - century * (0.004817 + 0.000014 * century)) + sin(2 * rad(mean.anom)) * (0.019993 - 0.000101 * century) + sin(3 * mean.anom) * 0.000289
  
  ## Calc lon
  lon = mean.lon + sun.center
  lon.app = lon - 0.00569 - 0.00478 * sin(rad(125.04 - 1934.136 * century))
  
  anom = mean.anom + sun.center
  
  sun.dist = (1.000001018 * (1 - eccentric^2)) / (1 + eccentric * cos(rad(anom)))
  mean.obl.ecl = 23 + (26 + ((21.448 - century * (46.815 + century * (0.00059 - century * 0.001813)))) / 60) / 60
  obl.corr = mean.obl.ecl + 0.00256 * cos(rad(125.04 - 1934.136 * century))
  
  sun.right.asc = deg(atan2(cos(rad(lon.app)), cos(rad(obl.corr)) * sin(rad(lon.app))))
  declin = deg(asin(sin(rad(obl.corr)) * sin(rad(lon.app))))
  
  y = tan(rad(obl.corr / 2))^2
  
  
  time = 4 * deg(y * sin(2 * rad(mean.lon)) - 2 * eccentric * sin(rad(mean.anom)) + 4 * eccentric * y * sin(rad(mean.anom)) * cos(2 * rad(mean.lon)) - 0.5* y^2 * sin(4*rad(mean.lon)) - 1.25 * eccentric^2 * sin(2*rad(mean.anom)))
  HA = 
  
  list(time = time, HA = HA)
}

times = function(sun, lon, lat, tz = 0) {
  noon = (720 - 4 * lon - sun$time + tz * 60) / 1440
  rise = noon - sun$HA * 4 / 1440
  
} 

