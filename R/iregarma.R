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


iregarma<-function ( rawdata , ar = 0 , ma = 0 , method = 'alasso' , normalize = FALSE , 
                     mselection = 'BIC' , alpha = 0 , df1 = 0 , df2 = 0 , ses.l = 1, 
                     auto.prune = TRUE , zero.find.force = FALSE , rep = 7 ,
                     pacfs = TRUE , BLACKLIST.ID = c (  ) , debug = FALSE
) 
{
  mpreg = FALSE
  internalRep = FALSE # this is just an indicator
  if ( debug == TRUE ) {Rprof ( tf <- "liregarma.log" , memory.profiling = TRUE ) } 
  ini.ar = ar; ini.ma = ma ; r = ncol(rawdata) - 2 
  x=list() ; y= list()
  fdt = table ( rawdata[ , 1])
  rfac= which(fdt < 3 + df1 + df2 + ses.l);  fdt = fdt[-rfac] ; cat(rfac, ' are removed due to lack of observations!' )
  fac = as.vector ( names ( fdt )  ) ; if (length(fac) == 0) stop('Each group must have more than 3 observations!',call. = FALSE)   
  
    #This is pruning loop and it is controled by the end of secuence
  for ( prune.iter in 1: ( rep*50 ) ) { 
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
    
    #print(fac)
    #if((is.numeric(fac) == FALSE) | (min(fac) <= 0) | (is.integer(fac) == FALSE) ) stop ('Factors must be positive numeric integers!',call. = FALSE)
    if ( is.null ( BLACKLIST.ID )  == FALSE && length ( BLACKLIST.ID )  > 0 ) {
      rElements = which ( fac %in% BLACKLIST.ID ) 
      if ( length ( rElements )  > 0 ) {
        fac = fac[- ( rElements ) ]
        Ncat ( c ( 'REMOVED GROUPS: ' , rElements , '\n' ) , FALSE ) 
      }
    }
    R = length ( fac ) 
    pb <- txtProgressBar ( min  =  1 , max  =  R+2 , style  =  3 ) 
    for  ( nf in 1:R ) {
      if (internalRep == FALSE){
        #######--------------------------
        xDframe = subset ( rawdata , rawdata[ , 1] == fac[nf] ) 
        tdata = as.matrix ( xDframe ) 
        
        tdata = cbind ( tdata[ , 1] , tdata[ , -1] ) 
        tdata = Ndiff ( tdata , df1 , df2 , ses.l)      
        
        if ( normalize == TRUE ) {
          x[[nf]] = scale ( as.matrix ( tdata[ , - ( 1:2 ) ] )  ) 
          y[[nf]] = scale ( as.matrix ( tdata[ , 2] )  ) 
        }else{
          x[[nf]] = as.matrix ( tdata[ , - ( 1:2 ) ] ) 
          y[[nf]] = as.matrix ( tdata[ , 2] ) 
        }
        
        # JUST SOME PACFs ARE SHOWN
        if  ( pacfs == TRUE ) {
          par ( mfrow = c ( round ( log ( length ( fac )  ) + 1 ) , round ( log ( length ( fac )  ) + 1 )  )  ) 
          acf  ( y[[nf]] , main=' ACF of y' )
          pacf ( y[[nf]] , main='PACF of y' ) 
        }
        
        Ncat ( c ( '\r ',nf , ' for the FACTOR:' , fac[nf] , '  ' , dim ( tdata ) , '' ) , display  =  1 ) 
        #######--------------------------
      }
      if ( ar > 0 ) {
        xar1 = as.matrix ( x[[nf]][- ( 1:ar ) , ] ) 
        yar1 = y[[nf]] #y[- ( 1:ar ) ]
        ar1matrix = as.matrix ( makematrix ( yar1 , p = ar , beforestart.mean = 0 , rname = 'HpAR1' ) [- ( 1:ar ) , ] ) 
        if( internalRep == FALSE){
          if(mpreg == FALSE ) {
            mar1 = msgps::msgps ( cbind ( xar1 , ar1matrix ) , as.vector ( yar1[- ( 1:ar ) ] ) , penalty = method , alpha = alpha , intercept = FALSE ) 
            coef = mar1$dfbic_result$coef 
          }else{
            mar1 = REGARMA:::mplm (x= cbind ( xar1 , ar1matrix ) , y = as.vector ( yar1[- ( 1:ar ) ] ) )
            coef = mar1$coef
          }
          hyar1 =  cbind ( xar1 , ar1matrix )   %*%  coef#  
        }else{
          new.y = yar1[- ( 1:ar ) ] - ar1matrix %*% vAR
          if(mpreg == FALSE){
            mar1 = msgps::msgps ( cbind ( xar1 ) , as.vector ( new.y ) , penalty = method , alpha = alpha , intercept = FALSE ) 
            coef = mar1$dfbic_result$coef 
          }else{
            mar1 = REGARMA:::mplm (x= cbind ( xar1 ) ,y = as.vector ( new.y ) ) 
            coef = mar1$coef 
          }
          
          hyar1 =  cbind ( xar1 )   %*%  coef# 
        }
        rar1  =  yar1[- ( 1:ar ) ] - hyar1
        if  ( pacfs == TRUE ) {
          pacf ( rar1 , main = paste('residual ACF in iteration /#', prune.iter %% rep) , ylim='AR residuals ACF') 
        }
      }else{
        xar1 = as.matrix ( x[[nf]] ) 
        yar1 = y[[nf]]
        if (mpreg == FALSE){
          mar1 = msgps::msgps ( cbind ( xar1 ) , as.vector ( yar1 ) , penalty = method , alpha = alpha , intercept = FALSE ) 
          coef = mar1$dfbic_result$coef 
        }else{
          mar1 = REGARMA:::mplm ( x= cbind ( xar1 ) , y= as.vector ( yar1 )  ) 
          coef = mar1$coef 
        }
        hyar1 =  cbind ( xar1 )   %*%  coef  
        rar1 = yar1 - hyar1
        if  ( pacfs == TRUE ) {
          pacf ( rar1 ) 
        }
      }
      ################ STEP 2
      if ( ma > 0 ) {
        if ( ar > 0 ) {
          xma1 = cbind ( xar1[- ( 1:ma ) , ] , ar1matrix[- ( 1:ma ) , ] )   
        }else{
          xma1 = cbind ( xar1[- ( 1:ma ) , ] )   
        }
        
        yma1 = y[[nf]][- ( 1: ( ar+ma )  ) ]
        ma1matrix = makematrix ( rar1 , p = ma , beforestart.mean = 0  , rname = 'HqMA1' ) [- ( 1:ma ) , ]
        H = cbind ( yma1 , xma1 , ma1matrix ) 
      }else{
        if ( ar > 0 ) {
          xma1 = cbind ( xar1 , ar1matrix )   
          yma1 = y[[nf]][- ( 1:ar ) ]
        }else{
          xma1 = cbind ( xar1 ) 
          yma1 = y[[nf]]
        }      
        H = cbind ( yma1 , xma1 ) 
      }
      FH = rbind ( FH , H ) 
      #print ( dim ( FH )  ) 
      setTxtProgressBar ( pb , nf ) 
    }
    FH = FH[-1 , ] #<-- remove the first row!
    par ( mfrow = c ( 1 , 1 )  ) 
    r1 = regarma ( x = FH[ , -1] , y = FH[ , 1] , ar = 0 , ma = 0 , method = method , mselection = mselection , alpha = alpha , auto.order = FALSE , auto.prune = FALSE) 
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
    ########
    if( prune.iter %% rep != 0 && ar>0 && ma>0) {
      print(vAR)
      print(vMA)
      internalRep = TRUE
      cat('\r', prune.iter %% rep, ' Recalculating in progress ...')
    }else{
      internalRep = FALSE
      if ( auto.prune  ==  TRUE)
      { 
        ar.temp = findlastZero ( x  =  vAR , force = zero.find.force)
        ma.temp = findlastZero ( x  =  vMA , force = zero.find.force)
        if (  ( ar.temp != ar ) |  ( ma.temp != ma ) )
        { 
          cat  ( '\n Pruning in progress... , AR: ', ar, '-to->', ar.temp, ' , MA: ', ma, '-to->', ma.temp , '\n' )
          ar = ar.temp
          ma = ma.temp
        }else{ 
          break
        } 
      }else{ 
        break
      }
    }
  }
  
  #if ( rep > 1 && ar > 0 && ma>0) { setTxtProgressBar ( pb3, ar+ma+1 ); close ( pb3 ) } 
  setTxtProgressBar ( pb , R+2 ) ;close ( pb ) 
  if ( debug == TRUE ) {Rprof ( NULL ) ;cat ( '\n------- DEBUG DATA -------\n' ) ;print ( summaryRprof ( tf )  ) }
  return ( list
           ( 
           #--------- Data ---------#
           y = r1$y , 
           x = r1$x , 
           Factors = fac , 
           Rawdata = rawdata , 
           nAR = ini.ar , nMA = ini.ma , 
           AR.Value = fillWithZero ( vAR, ini.ar ), 
           MA.Value = fillWithZero ( vMA, ini.ma ), 
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




