# Produce annual historic data
# Andy Purvis
#
# This code reads in historical LUH2 data (850-2015).
# Secondary is converted to an effective area of primary in each year, using bii_sec.
# They are harminized at the year where they most closely agree (1985).
# An annual time-series is produced and converted to an array.
# The array is written to a file. 
# Warning: this code block takes over 1 hour to run on my laptop

# Read in historical LUH2 data
pri <- stack("data/historical-primary.tif")
sec <- stack("data/historical-secondary.tif")
luh2_fraction <- clamp(pri + bii_sec * sec, lower=0, upper=1, useValues=TRUE)
luh2_years <- seq(historical_span[1], historical_span[2], by=historical_step)

# Preparing to meld the land-use time series
meld_year <- 1985 # Ascertained in previous analyses
v_target <- clamp(raster("data/vivid-base-primary.tif", layer=which(v_years==meld_year)) + 
                    bii_sec * raster("data/vivid-base-secondary.tif", layer=which(v_years==meld_year)), 
                    lower=0, upper=1, useValues=TRUE) * cell_area_raster
luh2_target <- subset(luh2_fraction, which(luh2_years==meld_year)) * cell_area_raster
luh2_area <- luh2_fraction * cell_area_raster
luh2_0 <- subset(luh2_area, 1)

diff <- v_target-luh2_target
plot(diff, main=paste("Vivid estimate - LUH2 estimate for", meld_year))

# Adjust luh2_fraction so that its meld_year values match those in v_base_fraction
luh.harmonized <- clamp(overlay(luh2_area, v_target, luh2_target, cell_area_raster,
                                fun=function(x, y, z, c){c-(c-x)*(c-y)/(c-z)}), lower=0, useValues=TRUE)
# Note that the calculation does lead to some slightly negative areas, when LUH2 data rose before 2018 but the Vivid 
# data for 1985 has natural habitat at or close to zero. The cells affected 

historic_area <- interpolateTemporal(s=luh.harmonized, xin=luh2_years, xout=c(min(luh2_years):max(luh2_years)), 
                                   outdir="figs", prefix="luh_annual_area", progress=TRUE, 
                                   writechange=FALSE, returnstack=TRUE, overwrite=TRUE) 
historic_annual_area <- as.array(historic_area)
saveRDS(historic_annual_area, file="output/historic_annual_area.Rds")

haa <- apply(historic_annual_area, MARGIN=3, sum, na.rm=TRUE)
plot(haa ~ c(min(luh2_years):max(luh2_years)), type="l", lwd=2, ylab="Natural habitat area", xlab="Year", cex.lab=1.3)
abline(v=1985, col="red", lty=2)
