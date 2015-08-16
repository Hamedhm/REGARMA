findlastZero <- function  ( x ,force = FALSE) { 
  lx = length  ( x ) 
  out.p = lx
  for  ( i in 1:lx ) { 
    if  ( round( x [ lx - i + 1], 6) == 0 ) { 
      out.p = out.p - 1
    } else { 
      if(force == FALSE){
        break
      }
    } 
  } 
  return  ( out.p ) 
} 
regARMA <- function  ( data,  ar = 0,  ma = 0,  method = c  ( "alasso" ) ,  
                       mselection  =  "BIC",  alpha = 0,  Ndf = 0,  Sdf = 0,  
                       ses.l = 1 , normalize = FALSE,  
                       with.rep = NA,  non.penalized = FALSE,  iteration = 1, 
                       auto.prune = TRUE, zero.find.force=FALSE , debug = FALSE ) { 
  auto.order  =  FALSE
  BLACKLIST.ID  =  c() 
  if  ( is.matrix  ( data ) == FALSE ) { stop  ( 'Input data must be a matrix!' ) } 
  
  if  ( non.penalized == FALSE ) { 
    repTable = table  ( data [ , 1] ) 
    MinReplication = min  ( repTable ) 
    MaxReplication = max  ( repTable ) 
    
    if  ( is.na  ( with.rep ) == TRUE ) { 
      if      ( ( MinReplication > 1 ) &&   ( MaxReplication > 1 ) ) { 
        with.rep = 1
        print  ( 'Replications are detected in dataset,  regARMA with replications is applied.' ) 
      } else { 
        with.rep = 0
        print ( 'No replication is detected in dataset,  regARMA without replications is applied.' ) 
      } 
    } 
    if ( with.rep == 1 ) { 
      cat('\n penalized REGARMA(',ar,',',ma,') with replications\n' )
      if ( auto.order == TRUE ) { stop ( 'Auto order is not active when there are replications! please set auto.order to FALSE.' ) } 
      result = iregarma ( rawdata  =  data,  ar  =  ar,  ma  =  ma,  
                          method  =  method,  normalize  =  normalize,  mselection  =  mselection,  
                          alpha  =  alpha,  df1  =  Ndf,  df2  =  Sdf,  ses.l=ses.l, pacfs  =  FALSE,
                          BLACKLIST.ID  =  BLACKLIST.ID,  debug = debug , auto.prune = auto.prune ,
                          zero.find.force = zero.find.force , rep = iteration) 
      return ( result )  
    } else if ( with.rep == 0 ) { 
      cat('\n penalized REGARMA(',ar,',',ma,')\n' )
      data = Ndiff ( data,  Ndf,  Sdf ,ses.l)  
      y = data [ ,  2]; x = data [ ,   - c ( 2, 1 ) ]  
      
      result = regarma ( x  =  x,  y  =  y,  ar  =  ar,  ma  =  ma,  method  =  method,  zero.find.force = zero.find.force,
                         mselection  = mselection,  alpha  =  alpha,  normalize = normalize,  
                         auto.order  =  auto.order,  auto.prune = auto.prune, rep  =  iteration, debug  =  debug ) 
      result  [[  'Rawdata' ]]  = data
      result  [[  'Factors' ]]  = data [ , 1]
      
      return ( result ) 
    } else { 
      stop ( 'with.rep must be 0, 1 or NA' ) 
    } 
    
  } else if ( non.penalized == TRUE ) { 
    cat('\n non-penalized REGARMA(',ar,',',ma,')\n' )
    data = Ndiff ( data,  Ndf,  Sdf , ses.l) 
    if ( nrow ( data ) <ncol ( data ) + ar + ma + 1 ) { stop ( 'Then number of observation must be greater than the number of variables' ) } 
    if ( auto.prune == TRUE ) { print ( 'Auto prune is not activated in this routine!' ) } 
    y = data [ ,  2]; x = data [ ,   - c ( 2, 1 ) ] 
    result = npregarma ( x  =  x, y  =  y, ar  =  ar, ma  =  ma,  normalize = normalize,  auto.order  =  auto.order,  debug = debug,  rep = iteration ) 
    result  [[  'Rawdata' ]]  = data
    result  [[  'Factors' ]]  = data [ , 1]
    return ( result ) 
    
  } else { 
    stop ( 'Are you kidding me!' ) 
  } 
} 
