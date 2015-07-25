\name{naics5809}
\alias{naics5809}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
    NBER-CES Manufacturing Industry Database(1958-2009)
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
    This database is a joint effort between the National Bureau of Economic Research (NBER) and U.S. Census Bureau's Center for Economic Studies (CES), containing annual industry-level data from 1958-2009 on output, employment, payroll and other input costs, investment, capital stocks, TFP, and various industry-specific price indexes. 
}
  \usage{data(naics5809)}
  \format{
  A data frame with 23868 observations on the following 24 variables.
  \describe{
  \item{\code{sic}}{: NAICS 4-digit Codes}
  \item{\code{year}}{: Year ranges from 1958 to 2009 }
  \item{\code{emp}}{: Total employment in 1000s }
  \item{\code{pay}}{: Total payroll in $1m}
  \item{\code{prode}}{: Production workers in 1000s }
  \item{\code{prodh}}{: Production worker hours in 1m }
  \item{\code{prodw}}{: Production worker wages in $1m}
  \item{\code{vship}}{: Total value of shipments in $1m }
  \item{\code{matcost}}{: Total cost of materials in $1m}
  \item{\code{vadd}}{: Total value added in $1m }
  \item{\code{invest}}{: Total capital expenditure in $1m }
  \item{\code{invent}}{: End-of-year inventories in $1m}
  \item{\code{energy}}{: Cost of electricity & fuels in $1m}
  \item{\code{cap}}{: Total real capital stock in $1m }
  \item{\code{equip}}{: Real capital: equipment in $1m}
  \item{\code{plant}}{: Real capital: structures in $1m}
  \item{\code{piship}}{: Deflator for VSHIP 1997=1.000}
  \item{\code{pimat}}{: Deflator for MATCOST 1997=1.000}
  \item{\code{piinv}}{: Deflator for INVEST 1997=1.000}
  \item{\code{pien}}{: Deflator for ENERGY 1997=1.000 }
  \item{\code{dtfp5}}{: 5-factor TFP annual growth rate }
  \item{\code{tfp5}}{: 5-factor TFP index 1997=1.000 }
  \item{\code{dtfp4}}{: 4-factor TFP annual growth rate }
  \item{\code{tfp4}}{: 4-factor TFP index 1997=1.000}
  }
  }
  \details{
  %%  ~~ If necessary, more details than the __description__ above ~~
  }
  \source{
  %%  ~~ reference to a publication or URL from which the data were obtained ~~
  \url{http://www.nber.org/nberces/}
  }
  \references{
  %%  ~~ possibly secondary sources and usages ~~
  Official website: \url{http://www.nber.org/data/nberces5809.html} \cr
  Variable description and summary statistics: \url{http://www.nber.org/nberces/nberces5809/nberces_5809_summary_stats.pdf} \cr
  Updated documentation: \url{http://www.nber.org/nberces/nberces5809/nberces_5809_technical_notes.pdf} \cr
  }
  \examples{
  data(naics5809)
  str(naics5809)
  }
  \keyword{datasets}
  
