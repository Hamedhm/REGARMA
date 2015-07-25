Ndiff<-function ( x , df1 = 0 , df2 = 0 , ses.l = 1) {
  if ( df1 > 0 ) {
    if ( df2 > 0 ) {
      return ( diff ( diff ( x , lag = ses.l , differences = df2 ), lag = 1, differences = df1 )  )
    }else{
      return ( diff ( x , lag = 1, differences =  df1 )  ) 
    }
  }else{
    if ( df2 > 0 ) {
      return ( diff ( x , lag = ses.l , differences = df2 )  ) 
    }else{
      return ( x ) 
    }
  }
}


apply.by.factor<-function ( x , fun = scale ) {
  x = as.matrix ( x ) 
  factors = names ( table ( x[ , 1] )  ) 
  ind = as.vector ( which ( table ( x[ , 1] )  == 1 )  ) 
  ind = which ( x[ , 1] == ind ) 
  
  z = by ( x[ , -1] , x[ , 1] , fun ) 
  lf = length ( names ( z )  ) 
  x2 =  ( z[[1]] ) 
  if ( lf > 1 ) {
    for  (  i in 2:lf ) {    
      x2 = rbind ( x2 , as.matrix ( z[[i]] )  )     
    }  
  }
  x2 = cbind ( x[ , 1] , x2 ) 
  if ( length ( ind )  > 0 ) {
    x2[ind , ] = x[ind , ]
  }
  x2 = as.matrix ( x2 ) 
  return ( x2 ) 
}

apply.diff.by.factor<-function ( x , d1 = 0 , d2 = 0, ses.l = 1) {    
  x = as.matrix ( x ) 
  factors = names ( table ( x[ , 1] )  )   
  Hz = matrix ( 0 , ncol =  ( ncol ( x ) -1 ) , nrow = 1 ) 
  #colnames ( Hz )  = colnames ( x ) 
  for  ( i in factors ) {
    z = subset ( x[ , -1] , x[ , 1] == i ) 
    z = Ndiff ( z , d1 , d2 ,ses.l ) 
    Hz = rbind ( Hz , z ) 
  }
  Hz = Hz[-1 , ]
  return ( Hz )   
}


iregarma<-function ( rawdata , ar = 2 , ma = 2 , method = 'enet' , normalize = TRUE , mselection = 'BIC' , alpha = 0 , df1 = 0 , df2 = 0 , ses.l = 1, pacfs = FALSE , BLACKLIST.ID = c (  ) , debug = FALSE ) {
  if ( debug == TRUE ) {Rprof ( tf <- "liregarma.log" , memory.profiling = TRUE ) } 
  FH = matrix ( 0 , ncol =  ( ncol ( rawdata ) +ar+ma-1 )  ) 
  if ( ar > 0 ) {
    if ( ma > 0 ) {
      colnames ( FH )  = c ( colnames ( rawdata ) [-1] , paste ( 'Hp' , 1:ar , sep = '.' ) , paste ( 'Hq' , 1:ma , sep = '.' )  )     
    }else{
      colnames ( FH )  = c ( colnames ( rawdata ) [-1] , paste ( 'Hp' , 1:ar , sep = '.' )  ) 
    }
  }else{
    if ( ma > 0 ) {
      colnames ( FH )  = c ( colnames ( rawdata ) [-1] , paste ( 'Hq' , 1:ma , sep = '.' )  ) 
    }else{
      colnames ( FH )  = c ( colnames ( rawdata ) [-1] ) 
    }
  }
  # Start spiting data into different groups   #
  fac = as.vector ( names ( table ( rawdata[ , 1] )  )  ) 
  if ( is.null ( BLACKLIST.ID )  == FALSE && length ( BLACKLIST.ID )  > 0 ) {
    rElements = which ( fac %in% BLACKLIST.ID ) 
    if ( length ( rElements )  > 0 ) {
      fac = fac[- ( rElements ) ]
      Ncat ( c ( 'REMOVED GROUPS: ' , rElements , '\n' ) , TRUE ) 
    }
  }
  R = length ( fac ) 
  pb <- txtProgressBar ( min  =  1 , max  =  R+2 , style  =  3 ) 
  for  ( nf in 1:R ) {
    xDframe = subset ( rawdata , rawdata[ , 1] == fac[nf] ) 
    tdata = as.matrix ( xDframe ) 
    
    tdata = cbind ( tdata[ , 1] , tdata[ , -1] ) 
    tdata = Ndiff ( tdata , df1 , df2 , ses.l)      
    
    if ( normalize == TRUE ) {
      x = scale ( as.matrix ( tdata[ , - ( 1:2 ) ] )  ) ;y = scale ( as.matrix ( tdata[ , 2] )  ) 
    }else{
      x = as.matrix ( tdata[ , - ( 1:2 ) ] ) ;y = as.matrix ( tdata[ , 2] ) 
    }
    
    # JUST SOME PACFs ARE SHOWN
    if  ( pacfs == TRUE ) {
      par ( mfrow = c ( round ( log ( length ( fac )  ) +1 ) , round ( log ( length ( fac )  ) +1 )  )  ) 
      pacf ( y ) 
    }
    
    Ncat ( c ( nf , ' for the FACTOR:' , fac[nf] , '  ' , dim ( tdata ) , '\n' ) , display  =  FALSE ) 
    if ( ar > 0 ) {
      xar1 = as.matrix ( x[- ( 1:ar ) , ] ) 
      yar1 = y #y[- ( 1:ar ) ]
      ar1matrix = as.matrix ( makematrix ( yar1 , p = ar , beforestart.mean = MASS::huber ( as.vector ( yar1 )  )$mu , rname = 'HpAR1' ) [- ( 1:ar ) , ] ) 
      mar1 = msgps::msgps ( cbind ( xar1 , ar1matrix ) , as.vector ( yar1[- ( 1:ar ) ] ) , penalty = method , alpha = alpha , intercept = FALSE ) 
      hyar1 =  cbind ( xar1 , ar1matrix )   %*%  mar1$dfbic_result$coef #  
      rar1 = yar1[- ( 1:ar ) ]-hyar1
      if  ( pacfs == TRUE ) {
        pacf ( rar1 ) 
      }
    }else{
      xar1 = as.matrix ( x ) 
      yar1 = y
      mar1 = msgps::msgps ( cbind ( xar1 ) , as.vector ( yar1 ) , penalty = method , alpha = alpha , intercept = FALSE ) 
      hyar1 =  cbind ( xar1 )   %*%  mar1$dfbic_result$coef  
      rar1 = yar1-hyar1
      if  ( pacfs == TRUE ) {
        pacf ( rar1 ) 
      }
    }
    
    
    if ( ma > 0 ) {
      if ( ar > 0 ) {
        xma1 = cbind ( xar1[- ( 1:ma ) , ] , ar1matrix[- ( 1:ma ) , ] )   
      }else{
        xma1 = cbind ( xar1[- ( 1:ma ) , ] )   
      }
      
      yma1 = y[- ( 1: ( ar+ma )  ) ]
      ma1matrix = makematrix ( rar1 , p = ma , beforestart.mean = MASS::huber ( as.vector ( rar1 )  )$mu , rname = 'HqMA1' ) [- ( 1:ma ) , ]
      H = cbind ( yma1 , xma1 , ma1matrix ) 
    }else{
      if ( ar > 0 ) {
        xma1 = cbind ( xar1 , ar1matrix )   
        yma1 = y[- ( 1:ar ) ]
      }else{
        xma1 = cbind ( xar1 ) 
        yma1 = y
      }      
      H = cbind ( yma1 , xma1 ) 
    }
    FH = rbind ( FH , H ) 
    #print ( dim ( FH )  ) 
    setTxtProgressBar ( pb , nf ) 
  }
  FH = FH[-1 , ] #<-- remove the first row!
  par ( mfrow = c ( 1 , 1 )  ) 
  r1 = regarma ( x = FH[ , -1] , y = FH[ , 1] , ar = 0 , ma = 0 , method = method , mselection = mselection , alpha = alpha , auto.order = FALSE ) 
  nb = length ( r1$Betas ) 
  if ( nb<1 ) {stop ( 'No coefficient or raw data are empty!' ) }
  if ( ma > 0 ) {
    if ( ar > 0 ) {
      vMA = r1$Betas[ ( nb-ma+1 ) :nb];
      vAR = r1$Betas[ ( nb-ar-ma+1 ) : ( nb-ma ) ]
      btas = r1$Betas[- (  ( nb-ar-ma+1 ) : ( nb )  ) , ]
    }else{
      vMA = r1$Betas[ ( nb-ma+1 ) :nb];vAR = 0  
      btas = r1$Betas[- (  ( nb-ma+1 ) :nb ) , ]
    }
  }else{
    if ( ar > 0 ) {
      vMA = 0; vAR = r1$Betas[ ( nb-ar+1 ) :nb]
      btas = r1$Betas[- (  ( nb-ar+1 ) :nb ) , ]
    }else{
      vAR = vMA = 0
      btas = r1$Betas
    }
  }
  setTxtProgressBar ( pb , R+2 ) ;close ( pb ) 
  if ( debug == TRUE ) {Rprof ( NULL ) ;cat ( '\n------- DEBUG DATA -------\n' ) ;print ( summaryRprof ( tf )  ) }
  return ( list
           ( 
           #--------- Data ---------#
           y = r1$y , 
           x = r1$x , 
           Factors = fac , 
           Rawdata = rawdata , 
           nAR = ar , nMA = ma , 
           AR.Value = vAR , MA.Value = vMA , 
           #--------- Model ---------#
           Selection.method = mselection , 
           Penalty.type = method , 
           Alpha = alpha , 
           #--------- Residuals ---------#
           Residuals = r1$Residuals , 
           #--------- Estimations ---------#
           Estimated.y = r1$Estimated.y , 
           #--------- Coefficients ---------#
           Betas = btas , 
           #--------- BIC ---------#
           Sel.Criteria.Steps = r1$Sel.Criteria.Steps , 
           #--------- MSE's ---------#
           MSE = r1$MSE , 
           MAE = r1$MAE
           #--------- Finish ---------#
           ) 
  ) 
}



