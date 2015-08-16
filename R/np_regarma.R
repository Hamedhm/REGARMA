`npregarma`<-function ( x, y, ar = 0, ma = 0, normalize = FALSE, auto.order = FALSE, debug = FALSE, rep = 1000 , stop.criteria = 10^ -10) {
  mpreg = FALSE  
  if ( debug == TRUE ) {Rprof ( tf <- "nplregarma.log",  memory.profiling = TRUE ) } 
  pb2 <- txtProgressBar ( min = 1, max = 3, style = 3 ) 
  if ( rep>1 && ar>0 && ma>0) {pb3 <- txtProgressBar ( min = 1, max = rep + 1, style = 3, label = 'iterations' ) }
  
  
  #------ Initializing ------  #  
  x = as.matrix ( x ) ; temX = x; temY = y
  if ( normalize == TRUE ) {x = scale ( x ) ; y = scale ( y ) }
  print.output = FALSE
  r = ncol ( x ) ; ly = length ( y ) 
  ar.value = c () ; ma.value = c ()  #initial values
  if ( ar + ma >= ly ) {stop ( 'AR-MA orders are higher than the number of datapoints!' ) }
  if  ( auto.order == TRUE ) {     # Auto order 1
    ic = c ( 'aic', 'bic', 'bic' )  [ 2 ] 
    ar.model = auto.arima ( x = y, d = 0, D = 0, max.p = 5, max.q = 0, max.P = 0, max.Q = 0, seasonal = TRUE, ic = ic  ) 
    ar = length ( ar.model$model$theta ) 
    ar.value = ar.model$model$theta
    Ncat ( display = print.output, c ( 'AR AIC,  BIC:', ar.model$aic, '|', ar.model$bic, '\n' )  ) 
    if  ( ar <= 0 || is.null ( ar )  == TRUE )  {
      ar = 0; ar.value = 0
    }
  }  
  
  ini.ar = ar; ini.ma = ma
  coef.matrix=matrix(c(1, 1, countTonumber(ini.ar),  countTonumber(ini.ma), countTonumber(r) ) ,ncol=1)
  
  if ( ar>0 ) {
    s1.ar.seq =  ( 1:ar )   
  }else {
    s1.ar.seq = 0
  }
  
  #--------- REGARMA STEP 1 ----------
  setTxtProgressBar ( pb2, 1 ) 
  if ( ar>0 ) {
    Hp = makematrix ( y, ar, 0, 'Hp' ) 
    if ( mpreg == FALSE ) {
      s1.regar = lm ( as.vector ( y [ -s1.ar.seq ]  ) ~cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  )  + 0 ) 
      s1.haty = cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  )  %*% s1.regar$coefficients # + s1.regar$intercept
    }else{
      #s1.regar = REGARMA:::mplm2 ( x = cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  ) , y = as.vector ( y [ -s1.ar.seq ]  )  ) 
      s1.regar = lm ( as.vector ( y [ -s1.ar.seq ]  ) ~cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  )  + 0 ) 
      s1.haty = cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  )  %*% s1.regar$coefficients 
    }
    s1.error = y [ -s1.ar.seq ] -s1.haty  # residuals
    #     plot ( s1.regar$dfbic_result$result, type = 'l' ) 
  }else{   
    if ( mpreg == FALSE ) {
      s1.lm.result = lm ( as.vector ( y ) ~cbind ( x )  + 0 ) 
      s1.haty = cbind ( x )  %*% s1.lm.result$coefficients   # +  s1.lm.result$intercept
    }else{
      #s1.lm.result = REGARMA:::mplm2 ( x = cbind ( x ) , y = as.vector ( y )  ) 
      s1.lm.result = lm ( as.vector ( y ) ~cbind ( x )  + 0 ) 
      s1.haty = cbind ( x )  %*% s1.lm.result$coefficients   # +  s1.lm.result$intercept
    }
    s1.error = y-s1.haty # residuals
    #     plot ( s1.lm.result$dfbic_result$result, type = 'l' ) 
  }
  # auto order for Y  
  if  ( auto.order == TRUE ) {
    s2.ma.model = auto.arima ( x = s1.error, d = 0, D = 0, max.p = 5, max.q = 0, max.P = 0, max.Q = 0, seasonal = TRUE, ic = ic  ) 
    ma = length ( s2.ma.model$model$theta ) 
    ma.value = s2.ma.model$model$theta
    Ncat ( display = print.output, c ( 'MA AIC BIC: ', s2.ma.model$aic, '|', s2.ma.model$bic, '\n' )  ) 
    if  ( ma <= 0 || is.null ( ma )  == TRUE )  {ma = 0; ma.value = 0}
  }
  if ( ma>0 ) {
    s2.ma.seq =  ( 1:ma ) 
  }else{
    s2.ma.seq = 0
  }
  if  ( ar>0 ) {
    if ( ma>0 ) {
      narPma =  ( 1: ( ar + ma )  )   
    }else{
      narPma =  ( 1:ar ) 
    }
  }else{
    if ( ma == 0 ) {
      narPma = 0
    }else{
      narPma =  ( 1:ma ) 
    }
  }
  
  setTxtProgressBar ( pb2, 2 ) 
  #--------- REGARMA STEP 2 -----#
  for ( j in 1:rep ) {    
    if ( ma>0 ) {
      if ( ar>0 ) {
        #Hq = makematrix ( s1.error, ma, MASS::huber ( as.vector ( s1.error )  )$mu*0, 'Hq' ) 
        Hq = makematrix ( s1.error, ma, 0, 'Hq' ) 
        s2.regarma = lm ( as.vector ( y [ -narPma ]  ) ~cbind ( Hp [ -narPma,  ] , Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  )  + 0 ) 
        s2.regarma.coeffs = s2.regarma$coefficients  
        s2.y.hat = cbind ( Hp [ -narPma,  ] , Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  )  %*% s2.regarma.coeffs # + s2.regarma$intercept
        s2.error = y [ -narPma ] -s2.y.hat
        HREGARMA = cbind ( Hp [ -narPma,  ] , Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  ) 
        #       plot ( s2.regarma$dfbic_result$result, type = 'l' ) 
      }else{
        #Hq = makematrix ( s1.error, ma, MASS::huber ( as.vector ( s1.error )  )$mu*0, 'Hq' ) 
        Hq = makematrix ( s1.error, ma, 0 , 'Hq' ) 
        s2.regarma = lm ( as.vector ( y [ -narPma ]  ) ~cbind ( Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  )  + 0 ) 
        s2.regarma.coeffs = s2.regarma$coefficients 
        s2.y.hat = cbind ( Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  )  %*% s2.regarma.coeffs # + s2.regarma$intercept
        s2.error = y [ -narPma ] -s2.y.hat      
        HREGARMA = cbind ( Hq [ -s2.ma.seq,  ] , x [ -narPma,  ]  ) 
        #       plot ( s2.regarma$dfbic_result$result, type = 'l' ) 
      }
    }else{
      if  ( ar>0 ) {   
        s2.regarma = lm ( as.vector ( y [ -narPma ]  ) ~cbind ( Hp [ -narPma,  ] , x [ -narPma,  ]  )  + 0 ) 
        s2.regarma.coeffs = s2.regarma$coefficients 
        s2.y.hat = cbind ( Hp [ -narPma,  ] , x [ -narPma,  ]  )  %*% s2.regarma.coeffs # + s2.regarma$intercept
        s2.error = y [ -narPma ] -s2.y.hat
        HREGARMA = cbind ( Hp [ -narPma,  ] , x [ -narPma,  ]  ) 
        #       plot ( s2.regarma$dfbic_result$result, type = 'l' ) 
      }else{
        s2.regarma = lm ( as.vector ( y ) ~cbind ( x )  + 0 ) 
        s2.regarma.coeffs = s2.regarma$coefficients
        s2.y.hat = cbind ( x )  %*% s2.regarma.coeffs # + s2.regarma$intercept
        s2.error = y-s2.y.hat
        HREGARMA = cbind ( x ) 
        #       plot ( s2.regarma$dfbic_result$result, type = 'l' ) 
      }
    }
    #setTxtProgressBar ( pb2, 3 ) 
    REGARMA.BIC = '' 
    #   plot ( REGARMA.BIC, type = 'l' ) 
    if ( ar>0 ) {
      if ( ma>0 ) {
        ar.value = as.vector ( s2.regarma.coeffs )  [ 1:ar ] 
        ma.value = as.vector ( s2.regarma.coeffs )  [  ( max ( ar )  + 1 ) : ( max ( ar )  + ma )  ] 
      }else{
        ar.value = as.vector ( s2.regarma.coeffs )  [ 1:ar ] 
        ma.value = 0
      }
    }else{
      if  ( ma>0 ) {
        ar.value = 0
        ma.value = as.vector ( s2.regarma.coeffs )  [ 1:ma ]     
      }else{
        ar.value = 0; ma.value = 0
        # print ( 'No ARMA coefficients' ) 
      }
    }
    Ncat ( display = print.output, c ( '\n AR-MA ORDERS: ', ar, '-', ma, '\n' )  ) 
    regarma.MSE = sum ( abs ( s2.error ) ^2 ) /length ( s2.error ) 
    regarma.MAE = sum ( abs ( s2.error )  )   /length ( s2.error ) 
    if ( ar>0 ) {
      Ncat ( display = print.output, c ( 'ROOT OF AR, MA COEFFICIENTS', '\n' )  ) 
      if ( mpreg == FALSE ) {
        ar.coeffs = s1.regar$coefficients [ 1:ar ] 
      }else{
        ar.coeffs = s1.regar$coefficients [ 1:ar ] 
      }
      ar.root =  ( abs ( polyroot ( c ( 1, ar.coeffs )  )  )  ) 
      if ( length ( ar.root ) >0 ) {
        ar.root = min ( ar.root ) 
      }else{
        ar.root = 'inf'
      }
      Ncat ( display = print.output, c ( 'Stationary of AR: ', ar.root, '\n' )  ) 
    }
    #setTxtProgressBar ( pb2, 4 ) 
    if ( ma>0 ) {    
      if ( ar>0 ) {
        mmodel = s2.regarma
        armmodel = mmodel$coefficients   [ 1:ar ] 
        mammodel = mmodel$coefficients  [  ( ar + 1 ) : ( ar + ma )  ]         
        ma.root =  ( abs ( polyroot ( c ( 1,  ( mammodel )  )  )  )  ) 
        if ( length ( ma.root ) >0 ) {
          ma.root = min ( ma.root ) 
        }else{
          ma.root = 'inf'
        }
        Ncat ( display = print.output, c ( 'Stationary of MA: ', ma.root, '\n' )  ) 
      }else{
        m2model = s2.regarma$coefficients  [ 1:ma ] 
        ma.root =  ( abs ( polyroot ( c ( 1, m2model )  )  )  ) 
        Ncat ( display = print.output, c ( ma.root, ' : POLYNOMIAL \n' )  ) 
        if ( length ( ma.root ) >0 ) {
          ma.root = min ( ma.root ) 
        }else{
          ma.root = 'inf'
        }      
        Ncat ( display = print.output, c ( 'Stationary of MA: ', ma.root, '\n' )  ) 
      }
    }
    if ( ar>0 ) {
      if  ( ma>0 ) {
        btas = as.matrix ( s2.regarma$coefficients   ) 
        if  ( length ( colnames ( x )  )  != 0 ) {
          rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , paste ( 'Hq.', 1:ma, sep = '' ) , colnames ( x )  ) 
        }else{
          rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , paste ( 'Hq.', 1:ma, sep = '' ) , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
        }
        btas = as.matrix ( btas ) ; btas = btas [ - ( 1: ( ar + ma )  ) ,  ]  #REMOVE AR/MA Coefficients
      }else{
        btas = as.matrix ( s2.regarma$coefficients   ) 
        if  ( length ( colnames ( x )  )  != 0 ) {
          rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , colnames ( x )  ) 
        }else{
          rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
        }       
        btas = as.matrix ( btas ) 
        btas = btas [ - ( 1: ( ar + ma )  ) ,  ]  #REMOVE AR/MA Coefficients
      }
    }else{
      if  ( ma>0 ) {
        btas = as.matrix ( s2.regarma$coefficients  ) 
        if  ( length ( colnames ( x )  )  != 0 ) {
          rownames ( btas ) <-c ( paste ( 'Hq.', 1:ma, sep = '' ) , colnames ( x )  ) 
        }else{
          rownames ( btas ) <-c ( paste ( 'Hq.', 1:ma, sep = '' ) , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
        }
        btas = as.matrix ( btas ) 
        btas = btas [ - ( 1: ( ar + ma )  ) ,  ]  #REMOVE AR/MA Coefficients
      }else{
        btas = as.matrix ( s2.regarma$coefficients  ) 
        if  ( length ( colnames ( x )  )  != 0 ) {
          rownames ( btas ) <-c ( colnames ( x )  ) 
        }else{
          rownames ( btas ) <-c ( paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
        }
        btas = as.matrix ( btas ) 
      }
    }
    # This repetition part
    new.coefs=c( 0, j,  fillWithZero ( ar.value, ini.ar ) , fillWithZero ( ma.value, ini.ma ), btas )
    coef.matrix=cbind ( coef.matrix , new.coefs )
    
    if ( ar>0 && rep>1 && ma>0) {
      s1.haty = cbind ( Hp [ -s1.ar.seq,  ] , x [ -s1.ar.seq,  ]  )  %*% c ( ar.value, btas ) ; 
      s1.error = y [ -s1.ar.seq ] -s1.haty ; setTxtProgressBar ( pb3, j ) 
      if(ncol(coef.matrix) > 2){
        stop.criteria.value = sum(abs( abs( new.coefs[-c(1,2)] ) - abs(coef.matrix[ -c(1,2), ncol(coef.matrix) - 1 ] ) )^2)
        cat('\r Stop criteria : ',stop.criteria.value,'')
        if (stop.criteria.value < stop.criteria) {
          setTxtProgressBar ( pb3, rep+1 )
          cat('\n Stop criteria (',stop.criteria,') reached in ',j,' iterations with squared error = ', stop.criteria.value ,'')
          break
        }
      }
    }else{
      break
    }
  }
  if ( rep>1 && ar>0 && ma>0) {  setTxtProgressBar ( pb3, rep + 1 ) ; close ( pb3 ) }
  setTxtProgressBar ( pb2, 3 ) ; close ( pb2 ) 
  
  if ( debug == TRUE ) {Rprof ( NULL ) ; cat ( '\n------- DEBUG DATA -------\n' ) ; print ( summaryRprof ( tf )  ) }
  return ( list
           ( 
           #--------- Data ---------#
           y = temY,  x = temX, 
           nAR = ar,  nMA = ma, 
           AR.Value = ar.value,  MA.Value = ma.value, 
           #--------- Model ---------#
           Selection.method = '- No Selection method -', 
           Penalty.type = '-No penalty-', 
           Alpha = '- No alpha -', 
           #--------- Residuals ---------#
           Residuals = s2.error, 
           #--------- Estimations ---------#
           Estimated.y = s2.y.hat, 
           #--------- Coefficients ---------#
           Betas = btas, 
           coeff.trend=as.matrix(coef.matrix[,-1]),
           #--------- BIC AIC and GCV ---------#
           Sel.Criteria.Steps = '- No results -', 
           #--------- MSE and MAE ---------#
           MSE = regarma.MSE, 
           MAE = regarma.MAE
           #--------- Finish ---------#
           ) 
  ) 
}
