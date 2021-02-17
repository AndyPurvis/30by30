# Produce annual historic data
# Andy Purvis
#
# This code reads in historical LUH2 data (850-2015).
# Secondary is converted to an effective area of primary in each year, using bii_sec.
# They are harmonized at the year where they most closely agree with the IMAGE data (1970).
# An annual time-series is produced and converted to an array.
# The array is written to a file. 
# Warning: this code block takes over 1 hour to run on my laptop

# Read in historical LUH2 data
pri <- stack("data/historical-primary.tif")
sec <- stack("data/historical-secondary.tif")
luh2_fraction <- clamp(pri + bii_sec * sec, lower=0, upper=1, useValues=TRUE)
luh2_fraction <- subset(luh2_fraction, c(1:225)) # Drop everything after 1970
luh2_years <- seq(850, 1970, by=5)

# Preparing to meld the land-use time series
meld_year <- 1970 # Ascertained in previous analyses
v_target <- subset(stack("output/ena/ENA-IMAGE-1970-2000.nc"), 1)
luh2_target <- subset(luh2_fraction, 225) * cell_areas
luh2_area <- luh2_fraction * cell_areas
luh2_0 <- subset(luh2_area, 1)

diff <- v_target-luh2_target
plot(diff, main="Vivid estimate - LUH2 estimate for 1970")

# Adjust luh2_fraction so that its meld_year values match those in v_base_fraction
luh.harmonized <- clamp(overlay(luh2_area, v_target, luh2_target, cell_areas,
                                fun=function(x, y, z, c){c-(c-x)*(c-y)/(c-z)}), lower=0, useValues=TRUE)
# Note that the calculation may lead to some slightly negative areas, when LUH2 data rose before the meld year but
# the corresponding IMAGE layer has natural habitat at or close to zero. 

plot(v_target, subset(luh.harmonized, 225))
abline(a=0,b=1,col="red")
plot(v_target - subset(luh.harmonized, 225), main="IMAGE - harmonized LUH2: 1970")
map(add=TRUE)

# concatenate the IMAGE historical data for 1980, 1990 and 2000, and the data for 2010, 2015 and 2020 from one IMAGE scenario, on to this
lu.harmonized <- stack(luh.harmonized, subset(stack("output/ena/ENA-IMAGE-1970-2000.nc"), c(2:4)), 
                        subset(stack("output/ena/ENA-IMAGE-RCPref_SSP2_BID-06Feb2020_FinalMask.nc"), c(1:3)))
lu.years <- c(luh2_years, 1980, 1990, 2000, 2010, 2015, 2020)
historic_area <- interpolateTemporal(s=lu.harmonized, xin=lu.years, xout=c(850:2020), 
                                   outdir="figs", prefix="lu_annual_area", progress=TRUE, 
                                   writechange=FALSE, returnstack=TRUE, overwrite=TRUE) 
if (require(ncdf4)) {	
  rnc <- writeRaster(historic_area, filename="output/historic_area.nc", format="CDF", overwrite=TRUE) 
}

historic_annual_area <- as.array(historic_area)
saveRDS(historic_annual_area, file="output/historic_annual_area.Rds")

haa <- apply(historic_annual_area, MARGIN=3, sum, na.rm=TRUE)
plot(haa ~ c(850:2020), type="l", lwd=2, ylab="Effective natural area", xlab="Year", cex.lab=1.3)
abline(v=1970, col="red", lty=2)
