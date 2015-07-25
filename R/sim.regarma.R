generateAR <- function ( n = 1, l = -1, u = 1, min.distance = .Machine$double.eps, sort.coeff = FALSE ) { 
  if ( n > 0 ) { 
    if ( ( ( u-l ) /n ) <=  min.distance ) { stop ( 'Data are not possible to be generated. please decrease min.distance!' ) } 
    repeat
    { 
      minroots  = 0
      if ( sort.coeff == FALSE ) { 
        ar = runif ( n, l, u ) 
      } else { 
        ar = sort ( runif ( n, l, u ), decreasing  =  TRUE ) 
      }       
      minroots  <-  min ( abs ( polyroot ( z = c ( 1, -ar ) ) ) ) 
      if ( n > 1 ) { 
        if ( minroots > 1 & min ( diff ( abs ( sort ( ar ) ) ) ) > min.distance ) { break  } 
      } else { 
        if ( minroots > 1 ) { break  } 
      } 
    } 
    return ( ar ) 
  } else { 
    return ( 0 ) 
  } 
} 


sim.regARMA <- function ( n = 5, beta = c ( .62 ), x.independent  = TRUE, phi = c ( .3 ), theta = c ( .5 ), var.error = 1, n.z.coeffs = 0, shuffle = FALSE, draw.plot = FALSE ) {    
  if ( n.z.coeffs > 0 ) { 
    beta = c ( beta, rep ( 0, n.z.coeffs ) ) 
    if ( shuffle == TRUE ) { 
      beta = sample ( beta, length ( beta ), replace  =  FALSE ) 
    } 
  } 
  if ( length ( theta ) == 0 ) { t2etha = .1  } else { t2etha = theta  } 
  if ( length ( phi ) == 0 ) { p2hi = .1  } else { p2hi = phi  } ;m = c ( Mod ( polyroot ( z = c ( 1, -p2hi ) ) ), Mod ( polyroot ( z = c ( 1, -t2etha ) ) ) ) 
  if ( length ( m ) > 0 && min ( m ) <= 1 ) stop ( 'Model is not stationary!' ) 
  
  p = length ( phi ) 
  q = length ( theta ) 
  lbeta = length ( beta ) 
  
  if ( x.independent  == TRUE ) { 
    x = rnorm ( n*length ( beta ),0, sd = sqrt( abs ( beta ) ) ) 
    x = matrix ( x, ncol = lbeta, nrow = n ) 
    x = rbind ( matrix ( 0, ncol = lbeta, nrow = p + q ), x ) 
    colnames ( x ) = paste ( 'X', 1:length ( beta ), sep = '.' ) 
  } else { 
    #if ( max ( abs ( beta ) ) > 1 ) stop ( 'Beta is not stationary ( <1 ) ' ) 
    x = matrix ( 0, ncol = lbeta, nrow = n ) 
    for ( i in 1:lbeta ) { 
      if ( beta [ i ] == 0 ) { 
        x [, i ] = rnorm ( n, 0, 1 ) 
      } else { 
        #x [, i ] = arima.sim ( n = n, list ( ar = beta [ i ] ) ) 
        x [, i ] = arima.sim ( n = n, list ( ar = runif(1) ) , sd = sqrt( abs ( beta ) ) ) 
      } 
    } 
    x = rbind ( matrix ( 0, ncol = lbeta, nrow = p + q ), x ) 
    colnames ( x ) = paste ( 'X', 1:length ( beta ), sep = '.' ) 
  } 
  y = rep ( 0, p + q ) 
  t0 = p + q
  distrurbance = c () 
  for ( t in 1:n ) { 
    #print ( t ) 
    rs = 0
    rs = x [ t0 + t, ]  %*% as.matrix ( beta ) 
    
    js = 0
    js = phi %*% y [ ( t + t0 - 1 ) : ( t + t0 - p ) ] 
    
    
    ks = 0
    ks = theta  %*% ( as.matrix ( y [ ( t + t0 - 1 ) : ( t + t0 - q ) ] )  - x [ ( t + t0 - 1 ) : ( t + t0 - q ), ]  %*% as.matrix ( beta ) ) 
    
    s2k = 0
    for ( l in 1:q ) { 
      for ( i in 1:p ) { 
        s2k = s2k + theta [ l ] *phi [ i ] *y [ t + t0 - l - i ] 
      } 
    } 
    ks = ks - s2k
    distrurbance [ t ] = rnorm ( 1, 0, sqrt ( var.error ) ) 
    y [ t + t0 ] = rs + js + ks + distrurbance [ t ] 
  }  
  
  y = y [  -  ( 1:t0 ) ] 
  x = x [  -  ( 1:t0 ), ] 
  #ts.plot ( y, main = paste ( 'AR ( ', toString ( round ( phi, 2 ) ), ' ) \nMA ( ', toString ( round ( theta, 2 ) ), sep = ' ', ' ) ' ) ) 
  if ( draw.plot == TRUE ) { 
    plot ( spline ( y ), type = 'l', ylab = 'y', xlab = 'Time', main = paste ( 'AR ( ', toString ( round ( phi, 2 ) ), ' ) \nMA ( ', toString ( round ( theta, 2 ) ), sep = ' ', ' ) ' ) ) 
  } 
  otpMatrix = cbind ( 1:length ( y ), y, x ) 
  colnames ( otpMatrix ) = c ( 'T', 'Y', paste ( 'X.', 1:lbeta, sep = '' ) ) 
  
  return ( list ( rawdata = otpMatrix, y = y, x = x, beta = beta, phi = phi, theta = theta, error = distrurbance ) ) 
} 
