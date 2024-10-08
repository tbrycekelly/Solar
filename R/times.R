#' Julian Day definition
julian = function(datetime) {
  day = as.numeric(difftime(datetime, as.POSIXct("1899-12-30 00:00:00", tz = 'GMT'), unit = 'days')) + 2415018.5
  century = (day - 2451545) / 36525
  
  list(day = day, century = century)
}
