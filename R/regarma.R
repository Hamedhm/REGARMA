fillWithZero <- function ( x, l ) { 
  nx = length ( x )
  if ( l > nx ) { 
    if ( nx != l ) { 
      otp = c ( x, rep ( 0, l-nx ) )
    }else{ 
      otp = x
    } 
  }else{ 
    otp = x    
  } 
  return ( otp )
} 

countTonumber<-function(x){
  if(x <= 0){
    return(0)
  }else{
    return(1:x)
  }
}

`regarma` <- function ( x, y, ar = 0, ma = 0, method = c ( 'enet', 'alasso' ), mselection = c ( 'CP', 'AIC', 'GCV', 'BIC' ), alpha = 0, normalize = FALSE, auto.order = FALSE, auto.prune = TRUE, rep = 15, stop.criteria = 10^ -6, zero.find.force=FALSE, debug = FALSE ) { 
  mpreg = FALSE  
  if ( debug  ==  TRUE ) { Rprof ( tf <- "lregarma.log", memory.profiling = TRUE ) }  
  pb2 <- txtProgressBar ( min  =  1, max  =  3, style  =  3 )
  if ( rep > 1 && ar > 0 && ma>0) { pb3 <- txtProgressBar ( min  =  1, max  =  rep+1, style  =  3, label  =  'iterations' ) } 
  
  
  #------ Initializing ------  #
  if ( length ( method ) > 1 || is.null ( method ) ==  TRUE ) { method = 'alasso'; print ( 'Penalty type is set to Adaptive-lasso.' ) } 
  if ( length ( mselection ) > 1 || is.null ( mselection ) ==  TRUE ) { mselection = 'CP'; print ( 'Selection method is set to Mallows CP.' ) } ;  
  
  selection = which ( mselection  ==  c ( 'CP', 'AIC', 'GCV', 'BIC' ) )
  if ( selection %in%  ( 1:4 ) ) { selection = selection+3 }else{ selection = 7 }  #set default on BIC
  x = as.matrix ( x ); temX = x; temY = y
  if ( normalize  ==  TRUE ) { x = scale ( x ); y = scale ( y ) } 
  print.output = FALSE
  r = ncol ( x ); ly = length ( y )
  ar.values = c (); ma.values = c () #initial values
  if ( ar+ma >=  ly ) { stop ( 'AR-MA orders are higher than the number of datapoints!' ) } 
  if  ( auto.order  ==  TRUE ) {      # Auto order 1
    ic = c ( 'aic', 'bic', 'bic' ) [ selection-3 ] 
    ar.auto.model = auto.arima ( x = y, d  =  0, D  =  0, max.p  =  5, max.q  =  0, max.P  =  0, max.Q  =  0, seasonal  =  TRUE, ic  = ic  )
    ar = length ( ar.auto.model$model$theta )
    ar.values = ar.auto.model$model$theta
    Ncat ( display = print.output, c ( 'AR AIC, BIC:', ar.auto.model$aic, '|', ar.auto.model$bic, '\n' ) )
    if  ( ar <= 0 || is.null ( ar ) ==  TRUE )  { 
      ar = 0; ar.values = 0
    } 
  } 
  ini.ar = ar; ini.ma = ma
  #This is pruning loop and it is controled by the end of secuence
  coef.matrix=matrix(c(1, 1, countTonumber(ini.ar),  countTonumber(ini.ma), countTonumber(r) ) ,ncol=1)
  for ( prune.iter in 1: ( ar+ma+1 ) ) { 
    if ( ar > 0 ) { 
      s1.ar.seq =  ( 1:ar ) 
    } else  { 
      s1.ar.seq = 0
    } 
    #--------- REGARMA STEP 1 ----------
    setTxtProgressBar ( pb2, 1 )
    if ( ar > 0 ) { 
      Hp = makematrix ( y, ar, MASS::huber ( as.vector ( y ) )$mu, 'Hp' )
      if ( mpreg  ==  FALSE ) { 
        s1.regar = msgps::msgps ( cbind ( Hp [ -s1.ar.seq, ], x [ -s1.ar.seq, ]), as.vector ( y [ -s1.ar.seq ] ), penalty = method, intercept = FALSE, alpha = alpha, STEP.max  =  10^6 )
        s1.y.hat = cbind ( Hp [ -s1.ar.seq, ], x [ -s1.ar.seq, ]) %*% s1.regar [[ selection ]]$coef #+s1.regar$intercept
      }else{ 
        s1.regar = mplm ( x = cbind ( Hp [ -s1.ar.seq, ], x [ -s1.ar.seq, ]), y = as.vector ( y [ -s1.ar.seq ] ) )
        s1.y.hat = cbind ( Hp [ -s1.ar.seq, ], x [ -s1.ar.seq, ]) %*% s1.regar$coef 
      } 
      s1.error = y [ -s1.ar.seq ] -s1.y.hat  # residuals
      #     plot ( s1.regar$dfbic_result$result, type = 'l' )
    }else{    
      if ( mpreg  ==  FALSE ) { 
        lm.result = msgps::msgps ( cbind ( x ), as.vector ( y ), penalty = method, intercept = FALSE, alpha = alpha, STEP.max  =  10^6 )
        s1.y.hat = cbind ( x ) %*% lm.result [[ selection ]]$coef   #+ lm.result$intercept
      }else{ 
        lm.result = mplm ( x = cbind ( x ), y = as.vector ( y ) )
        s1.y.hat = cbind ( x ) %*% lm.result$coef   #+ lm.result$intercept
      } 
      s1.error = y-s1.y.hat # residuals
      #     plot ( lm.result$dfbic_result$result, type = 'l' )
    } 
    # auto order for Y  
    if  ( auto.order   ==   TRUE ) { 
      auto.ma = auto.arima ( x = s1.error, d  =  0, D  =  0, max.p  =  5, max.q  =  0, max.P  =  0, max.Q  =  0, seasonal  =  TRUE, ic  = ic  )
      ma = length ( auto.ma$model$theta )
      ma.values = auto.ma$model$theta
      Ncat ( display = print.output, c ( 'MA AIC BIC: ', auto.ma$aic, '|', auto.ma$bic, '\n' ) )
      if  ( ma <= 0 || is.null ( ma ) ==  TRUE )  { ma = 0; ma.values = 0 } 
    } 
    
    
    
    if ( ma > 0 ) { 
      ma.seq =  ( 1:ma )
    }else{ 
      ma.seq = 0
    } 
    if  ( ar > 0 ) { 
      if ( ma > 0 ) { 
        ar.ma.seq =  ( 1: ( ar+ma ) ) 
      }else{ 
        ar.ma.seq =  ( 1:ar )
      } 
    }else{ 
      if ( ma  ==  0 ) { 
        ar.ma.seq = 0
      }else{ 
        ar.ma.seq =  ( 1:ma )
      } 
    } 
    
    setTxtProgressBar ( pb2, 2 )
    #--------- REGARMA STEP 2 -----#
    for ( j in 1:rep ) { 
      if ( ma > 0 ) { 
        if ( ar > 0 ) { 
          Hq = makematrix ( s1.error, ma, MASS::huber ( as.vector ( s1.error ) )$mu, 'Hq' )
          s2.regarma.model = msgps::msgps ( cbind ( Hp [ -ar.ma.seq, ], Hq [ -ma.seq, ], x [ -ar.ma.seq, ]), as.vector ( y [ -ar.ma.seq ] ), penalty = method, intercept = FALSE, alpha = alpha, STEP.max  =  10^6 )
          s2.regarma.coeffs = s2.regarma.model [[ selection ]]$coef  
          s2.regarma.y.hat =  cbind ( Hp [ -ar.ma.seq, ], Hq [ -ma.seq, ], x [ -ar.ma.seq, ]) %*% s2.regarma.coeffs #+s2.regarma.model$intercept
          regarma.error = y [ -ar.ma.seq ] -s2.regarma.y.hat
          HREGARMA  =  cbind ( Hp [ -ar.ma.seq, ], Hq [ -ma.seq, ], x [ -ar.ma.seq, ])
          #       plot ( s2.regarma.model$dfbic_result$result, type = 'l' )
        }else{ 
          Hq = makematrix ( s1.error, ma, MASS::huber ( as.vector ( s1.error ) )$mu, 'Hq' )
          s2.regarma.model = msgps::msgps ( cbind ( Hq [ -ma.seq, ], x [ -ar.ma.seq, ]), as.vector ( y [ -ar.ma.seq ] ), penalty = method, intercept = FALSE, alpha = alpha, STEP.max  =  10^6 )
          s2.regarma.coeffs = s2.regarma.model [[ selection ]]$coef  
          s2.regarma.y.hat =  cbind ( Hq [ -ma.seq, ], x [ -ar.ma.seq, ]) %*% s2.regarma.coeffs #+s2.regarma.model$intercept
          regarma.error = y [ -ar.ma.seq ] -s2.regarma.y.hat      
          HREGARMA  =  cbind ( Hq [ -ma.seq, ], x [ -ar.ma.seq, ])
          #       plot ( s2.regarma.model$dfbic_result$result, type = 'l' )
        } 
      }else{ 
        if  ( ar > 0 ) {    
          s2.regarma.model = msgps::msgps ( cbind ( Hp [ -ar.ma.seq, ], x [ -ar.ma.seq, ]), as.vector ( y [ -ar.ma.seq ] ), penalty = method, intercept = FALSE, alpha = 0, STEP.max  =  10^6 )
          s2.regarma.coeffs = s2.regarma.model [[ selection ]]$coef  
          s2.regarma.y.hat =  cbind ( Hp [ -ar.ma.seq, ], x [ -ar.ma.seq, ]) %*% s2.regarma.coeffs #+s2.regarma.model$intercept
          regarma.error = y [ -ar.ma.seq ] -s2.regarma.y.hat
          HREGARMA  =  cbind ( Hp [ -ar.ma.seq, ], x [ -ar.ma.seq, ])
          #       plot ( s2.regarma.model$dfbic_result$result, type = 'l' )
        }else{ 
          s2.regarma.model = msgps::msgps ( cbind ( x ), as.vector ( y ), penalty = method, intercept = FALSE, alpha = alpha, STEP.max  =  10^6 )
          s2.regarma.coeffs = s2.regarma.model [[ selection ]]$coef  
          s2.regarma.y.hat =  cbind ( x ) %*% s2.regarma.coeffs #+s2.regarma.model$intercept
          regarma.error = y-s2.regarma.y.hat
          HREGARMA  =  cbind ( x )
          #       plot ( s2.regarma.model$dfbic_result$result, type = 'l' )
        } 
      }   
      #setTxtProgressBar ( pb2, 3 )
      REGARMA.BIC = s2.regarma.model [[ selection ]]$result  
      #   plot ( REGARMA.BIC, type = 'l' )
      if ( ar > 0 ) { 
        if ( ma > 0 ) { 
          ar.values = as.vector ( s2.regarma.coeffs ) [ 1:ar ] 
          ma.values = as.vector ( s2.regarma.coeffs ) [ ( max ( ar )+1 ): ( max ( ar )+ma ) ] 
        }else{ 
          ar.values = as.vector ( s2.regarma.coeffs ) [ 1:ar ] 
          ma.values = 0
        } 
      }else{ 
        if  ( ma > 0 ) { 
          ar.values = 0
          ma.values = as.vector ( s2.regarma.coeffs ) [ 1:ma ]   
        }else{ 
          ar.values = 0; ma.values = 0
          # print ( 'No ARMA coefficients' )
        } 
      } 
      Ncat ( display = print.output, c ( '\n AR-MA ORDERS: ', ar, '-', ma, '\n' ) )
      regarma.MSE = sum ( abs ( regarma.error )^2 )/length ( regarma.error )
      regarma.MAE = sum ( abs ( regarma.error ) ) /length ( regarma.error )
      if ( ar > 0 ) { 
        Ncat ( display = print.output, c ( 'ROOT OF AR, MA COEFFICIENTS', '\n' ) )
        if ( mpreg  ==  FALSE ) { 
          ar.values.root = s1.regar [[ selection ]]$coef [ 1:ar ] 
        }else{ 
          ar.values.root = s1.regar$coef [ 1:ar ] 
        } 
        ar.value.for.root =  ( abs ( polyroot ( c ( 1, ar.values.root ) ) ) )
        if ( length ( ar.value.for.root ) > 0 ) { 
          ar.value.for.root = min ( ar.value.for.root )
        }else{ 
          ar.value.for.root = 'inf'
        } 
        Ncat ( display = print.output, c ( 'Stationary of AR: ', ar.value.for.root, '\n' ) )
      } 
      #setTxtProgressBar ( pb2, 4 )
      if ( ma > 0 ) {     
        if ( ar > 0 ) { 
          ma.values.root = s2.regarma.model
          #armmodel = ma.values.root [[ selection ]]$coef   [ 1:ar ] 
          ma.values.in.root = ma.values.root [[ selection ]]$coef   [ ( ar+1 ): ( ar+ma ) ]       
          min.root =  ( abs ( polyroot ( c ( 1,  ( ma.values.in.root ) ) ) ) )
          if ( length ( min.root ) > 0 ) { 
            min.root = min ( min.root )
          }else{ 
            min.root = 'inf'
          } 
          Ncat ( display = print.output, c ( 'Stationary of MA: ', min.root, '\n' ) )
        }else{ 
          m2model = s2.regarma.model [[ selection ]]$coef   [ 1:ma ] 
          min.root =  ( abs ( polyroot ( c ( 1, m2model ) ) ) )
          Ncat ( display = print.output, c ( min.root, ' : POLYNOMIAL \n' ) )
          if ( length ( min.root ) > 0 ) { 
            min.root = min ( min.root )
          }else{ 
            min.root = 'inf'
          }       
          Ncat ( display = print.output, c ( 'Stationary of MA: ', min.root, '\n' ) )
        } 
      } 
      if ( ar > 0 ) { 
        if  ( ma > 0 ) { 
          btas = as.matrix  (  s2.regarma.model [[ selection ]]$coef   )
          if   (  length  (  colnames  (  x ) ) !=  0 ) { 
            rownames  (  btas ) <- c ( paste ( 'Hp.', 1:ar, sep = '' ), paste ( 'Hq.', 1:ma, sep = '' ), colnames ( x ) )
          }else{ 
            rownames ( btas ) <- c ( paste ( 'Hp.', 1:ar, sep = '' ), paste ( 'Hq.', 1:ma, sep = '' ), paste ( 'Var.', 1:ncol ( x ), sep = '' ) )
          } 
          btas = as.matrix ( btas ); btas = btas [ - ( 1: ( ar+ma ) ), ]#REMOVE AR/MA Coefficients
        }else{ 
          btas = as.matrix ( s2.regarma.model [[ selection ]]$coef   )
          if  ( length ( colnames ( x ) ) !=  0 ) { 
            rownames ( btas ) <- c ( paste ( 'Hp.', 1:ar, sep = '' ), colnames ( x ) )
          }else{ 
            rownames ( btas ) <- c ( paste ( 'Hp.', 1:ar, sep = '' ), paste ( 'Var.', 1:ncol ( x ), sep = '' ) )
          }        
          btas = as.matrix ( btas )
          btas = btas [ - ( 1: ( ar+ma ) ), ]#REMOVE AR/MA Coefficients
        } 
      }else{ 
        if  ( ma > 0 ) { 
          btas = as.matrix ( s2.regarma.model [[ selection ]]$coef   )
          if  ( length ( colnames ( x ) ) !=  0 ) { 
            rownames ( btas ) <- c ( paste ( 'Hq.', 1:ma, sep = '' ), colnames ( x ) )
          }else{ 
            rownames ( btas ) <- c ( paste ( 'Hq.', 1:ma, sep = '' ), paste ( 'Var.', 1:ncol ( x ), sep = '' ) )
          } 
          btas = as.matrix ( btas )
          btas = btas [ - ( 1: ( ar+ma ) ), ]#REMOVE AR/MA Coefficients
        }else{ 
          btas = as.matrix ( s2.regarma.model [[ selection ]]$coef   )
          if  ( length ( colnames ( x ) ) !=  0 ) { 
            rownames ( btas ) <- c ( colnames ( x ) )
          }else{ 
            rownames ( btas ) <- c ( paste ( 'Var.', 1:ncol ( x ), sep = '' ) )
          } 
          btas = as.matrix ( btas )
        } 
      }
      new.coefs=c( prune.iter-1, j,  fillWithZero ( ar.values, ini.ar ) , fillWithZero ( ma.values, ini.ma ), btas )
      coef.matrix=cbind ( coef.matrix , new.coefs )
      # This repetition part
      if ( ar > 0 && rep > 1 && ma>0) { 
        s1.y.hat = cbind ( Hp [ -s1.ar.seq, ], x [ -s1.ar.seq, ]) %*% c ( ar.values, btas )
        s1.error = y [ -s1.ar.seq ] -s1.y.hat 
        setTxtProgressBar ( pb3, j ); 
        if(ncol(coef.matrix) > 2){
          stop.criteria.value = sum(abs( abs( new.coefs[-c(1,2)] ) - abs(coef.matrix[ -c(1,2), ncol(coef.matrix) - 1 ] ) )^2)
          cat('\r Stop criteria : ',stop.criteria.value,'')
          #print(new.coefs[-c(1,2)])
          #print(coef.matrix[ -c(1,2), ncol(coef.matrix) - 1 ])
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
    if ( auto.prune  ==  TRUE & auto.order  ==  FALSE )
    { 
      ar.temp = findlastZero ( x  =  ar.values , force = zero.find.force)
      ma.temp = findlastZero ( x  =  ma.values , force = zero.find.force)
      if (  ( ar.temp != ar ) |  ( ma.temp != ma ) )
      { 
        if ( ar > 0 && rep > 1 && ma>0) { setTxtProgressBar ( pb3, rep+1 ) }         
        setTxtProgressBar ( pb2, 3 )
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
  
  if ( rep > 1 && ar > 0 && ma>0) { setTxtProgressBar ( pb3, rep+1 ); close ( pb3 ) } 
  setTxtProgressBar ( pb2, 3 ); close ( pb2 )
  
  
  if ( debug  ==  TRUE ) { Rprof ( NULL ); cat ( '\n------- DEBUG DATA -------\n' ); print ( summaryRprof ( tf ) ) } 
  return ( list
           ( 
           #--------- Data ---------#
           y = temY, x = temX, 
           nAR = ini.ar, nMA = ini.ma, 
           AR.Value = fillWithZero ( ar.values, ini.ar ), MA.Value = fillWithZero ( ma.values, ini.ma ), 
           #--------- Model ---------#
           Selection.method = mselection, 
           Penalty.type = method, 
           Alpha = alpha, 
           #--------- Residuals ---------#
           Residuals = regarma.error, 
           #--------- Estimations ---------#
           Estimated.y = s2.regarma.y.hat, 
           #--------- Coefficients ---------#
           Betas = btas, 
           coeff.trend=as.matrix(coef.matrix[,-1]),
           #--------- BIC AIC and GCV ---------#
           Sel.Criteria.Steps = REGARMA.BIC, 
           #--------- MSE and MAE ---------#
           MSE = regarma.MSE, 
           MAE = regarma.MAE
           #--------- Finish ---------#
           )
  )
} 


