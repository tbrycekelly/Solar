library(Solar)
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
                   datetime = seq(as.POSIXct('2022-01-01'), as.POSIXct('2023-12-31'), by = 30))


grid = calcSun(lon = grid$lon, lat = grid$lat, datetime = grid$datetime)

grid$insolation = grid$insolation * cos(pi / 180 * (30 - grid$solar.elevation)) * cos(pi / 180 * (180 - grid$solar.azimuth))
#grid = grid[grid$insolation > 0.1, ]

#### Calculate location in the sky (i.e. count)
product = SimpleGridder::buildGrid(xlim = c(0, 360), ylim = c(0, 90), nx = 160, ny = 100, x.factor = 1, y.factor = 1)
product = SimpleGridder::setGridder(product, gridder = gridCount, neighborhood = 400, weight.func = function(x) {x < 1})
product = SimpleGridder::appendData(product,
                                    x = grid$solar.azimuth,
                                    y = grid$solar.elevation,
                                    z = grid$insolation,
                                    label = 'insolation')
product = SimpleGridder::interpData(product)


zlim = c(0, max(pretty(product$interp$insolation)))
SimpleGridder::plotGrid(product, 'insolation', pal = pals::inferno(8), zlim = zlim)
mtext(paste(zlim, collapse = ' -> '))


#### Calculate insolation in the sky (i.e. count)
light = SimpleGridder::buildGrid(xlim = c(0, 360), ylim = c(0, 90), nx = 160, ny = 100, x.factor = 1, y.factor = 1)
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

zlim = c(1, max(pretty(light$interp$light)))
SimpleGridder::plotGrid(light, 'light', pal = pals::inferno(8), zlim = zlim, ztrim = c(NA, NA), xlim = c(75, 285), ylim = c(0, 75))
mtext(paste(zlim, collapse = ' -> '))
SimpleGridder::addContour(light, 'light', lwd = 3, col = 'white', levels = c(100, 200, 300))

light$interp$light[light$interp$light < 0] = 0

sum(light$interp$light[light$interp$light > 100]) / sum(light$interp$light)
sum(light$interp$light[light$interp$light > 200]) / sum(light$interp$light)
sum(light$interp$light[light$interp$light > 300]) / sum(light$interp$light)

sum(light$interp$light[light$grid$y > 10]) / sum(light$interp$light)
sum(light$interp$light[light$grid$y > 20]) / sum(light$interp$light)
sum(light$interp$light[light$grid$y > 30]) / sum(light$interp$light)


