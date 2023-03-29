

julian.day = function(datetime) {
  as.numeric(difftime(datetime, as.POSIXct("1899-12-30 00:00:00", tz = 'GMT'), unit = 'days')) + 2415018.5
}

julian.century = function(datetime) {
  (julian.day(datetime) - 2451545) / 36525
}

rad = function(deg) {
  deg * 3.1415926535 / 180
}

deg = function(rad) {
  rad * 180 / 3.1415926535
}
