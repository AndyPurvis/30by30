# Harmonize IAM data with LUH2 using 2010 as the harmonization years
# Andy Purvis
#
# This code reads in historical LUH2 data (850-2015).
# Secondary is converted to an effective area of primary in each year, using bii_sec.
# They are harminized at the earliest possible year (2010)
# An annual time-series is produced and converted to an array.
# The array is written to a file. 
# Warning: this code block takes over 1 hour to run on my laptop

harmonize_with_luh2_2010 <- function(iam_file, which_iam){
  stopifnot(
    which_iam %in% c("AIM", "GLOBIOM", "IMAGE", "MAGPIE"),
    file.exists(iam_file)
  )
  
  # Read in historical LUH2 data
  luh2_ena <- stack("output/ENA_LUH2.nc")
  luh2_years <- seq(850, 2015, by=5)
  names(luh2_ena) <- luh2_years
  
  # Preparing to meld the land-use time series
  meld_year <- 2010 
  iam_ena <- stack(iam_file)
  iam_target <- subset(iam_ena, 1) #It is the first layer we harmonize on
  luh2_0 <- subset(luh2_ena, 1)
  luh2_target <- subset(luh2_ena, which(luh2_years == 2010))
  
  diff <- iam_target-luh2_target
  plot(diff, main="ENA: IAM - LUH2 for 2010", sub=iam_file)
  
  # Adjust luh2_fraction so that its meld_year values match those in v_base_fraction
  luh.harmonized <- clamp(overlay(luh2_ena, iam_target, luh2_target, cell_area_raster,
                                  fun=function(x, y, z, c){c-(c-x)*(c-y)/(c-z)}), lower=0, useValues=TRUE)
  # Note that the calculation does lead to some slightly negative areas, when LUH2 data rose before 2018 but the Vivid 
  # data for 1985 has natural habitat at or close to zero. The cells affected 
  
  historic_area <- interpolateTemporal(s=luh.harmonized, xin=luh2_years, xout=c(min(luh2_years):max(luh2_years)), 
                                       outdir="figs", prefix="luh_annual_area", progress=TRUE, 
                                       writechange=FALSE, returnstack=TRUE, overwrite=TRUE) 
  historic_annual_area <- as.array(historic_area)
  saveRDS(historic_annual_area, file=paste("output/historic_annual_area_", which_iam, ".Rds", sep=""))
  
  haa <- apply(historic_annual_area, MARGIN=3, sum, na.rm=TRUE)
  plot(haa ~ c(min(luh2_years):max(luh2_years)), type="l", lwd=2, ylab="Natural habitat area", xlab="Year", cex.lab=1.3)
  abline(v=2010, col="red", lty=2)

  return(historic_annual_area)
}
