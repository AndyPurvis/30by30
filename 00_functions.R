# Functions needed for analysis for HM Treasury
# Andy Purvis

interpolated.species.at.t <- function(t, S0, k, c, z, Area)		
{
  ## This function finds the number of species at time t given initial species number S0, the history of Area, 
  ## the values of c and z in the species-area relationship and the k-parameter that governs relaxation rate. 
  ## The function originally came from the R script provided to AP by Isabel Rosa on 30th April 2020. Note that 
  ## the function does not contain anything within it that prevents S rising over time, though that is not 
  ## biologically possible for the endemic species that are the sole focus of this analysis. Therefore, when this 
  ##  function is called, it is called in a way that prevents a year's number of species being higher than the 
  ##previous year's number.
  
  if (is.numeric(Area) & !is.na(sum(Area)) & (max(Area, na.rm=TRUE)>0)) #handling of Area data that are absent, include NA or are all zero
  {  
    integrand <- function(tao, Area, k, c, z)
    { 
      fraction.tao <- tao - floor(tao)
      interpolated.area <- (1 - fraction.tao) * Area[floor(tao + 1)] + fraction.tao * Area[floor(tao + 2)]    
      return(k * c * (interpolated.area ^ z) * exp(k * tao)) 
    }
    integral <- integrate(integrand, lower=0, upper=t, Area, k, c, z, stop.on.error=FALSE)  	#Integrand must be a function
    return(S0 * exp(-k * t) + exp(-k * t) * integral$value)
  }
  else (NA)
}


add_res_fraction <-function(natural, restored, biis, plots=TRUE){
  ## Combines BII-discounted restored land with natural land
  ## natural = rasterStack of natural habitat per grid cell for the years in v_years
  ## restored = rasterStack of the fraction restored per grid cell within each 5-year period from 2020 to 2060
  ## biis = set of values of how BII increases with age of restored habitat
  
  #TODO - Add feature-rich stopifnot
  
  #Force predictable and meaningful names on the rasterStacks, as the names provided are not consistent
  names(natural) <- paste("N",v_years, sep="_")
  names(restored) <- paste("R", v_res_years, sep="_")
  
  # Add restored land to natural fraction, discounted by its age-specific BII
  natural$N_2020 <- natural$N_2020 + biis[1] * restored$R_2020
  natural$N_2025 <- natural$N_2025 + biis[2] * restored$R_2020 +
    biis[1] * restored$R_2025
  natural$N_2030 <- natural$N_2030 + biis[3] * restored$R_2020 +
    biis[2] * restored$R_2025 + biis[1] * restored$R_2030
  natural$N_2035 <- natural$N_2035 + biis[4] * restored$R_2020 +
    biis[3] * restored$R_2025 + biis[2] * restored$R_2030 + biis[1] * restored$R_2035
  natural$N_2040 <- natural$N_2040 + biis[5] * restored$R_2020 +
    biis[4] * restored$R_2025 + biis[3] * restored$R_2030 + biis[2] * restored$R_2035 +
    biis[1] * restored$R_2040
  natural$N_2045 <- natural$N_2045 + biis[6] * restored$R_2020 +
    biis[5] * restored$R_2025 + biis[4] * restored$R_2030 + biis[3] * restored$R_2035 +
    biis[2] * restored$R_2040 + biis[1] * restored$R_2045
  natural$N_2050 <- natural$N_2050 + biis[7] * restored$R_2020 +
    biis[6] * restored$R_2025 + biis[5] * restored$R_2030 + biis[4] * restored$R_2035 +
    biis[3] * restored$R_2040 + biis[2] * restored$R_2045 + biis[1] * restored$R_2050
  natural$N_2055 <- natural$N_2055 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[6] * restored$R_2030 + biis[5] * restored$R_2035 +
    biis[4] * restored$R_2040 + biis[3] * restored$R_2045 + biis[2] * restored$R_2050 +
    biis[1] * restored$R_2055
  natural$N_2060 <- natural$N_2060 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[7] * restored$R_2030 + biis[6] * restored$R_2035 +
    biis[5] * restored$R_2040 + biis[4] * restored$R_2045 + biis[3] * restored$R_2050 +
    biis[2] * restored$R_2055 + biis[1] * restored$R_2060
  natural$N_2070 <- natural$N_2070 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[7] * restored$R_2030 + biis[7] * restored$R_2035 +
    biis[7] * restored$R_2040 + biis[6] * restored$R_2045 + biis[5] * restored$R_2050 +
    biis[4] * restored$R_2055 + biis[3] * restored$R_2060
  natural$N_2080 <- natural$N_2080 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[7] * restored$R_2030 + biis[7] * restored$R_2035 +
    biis[7] * restored$R_2040 + biis[7] * restored$R_2045 + biis[7] * restored$R_2050 +
    biis[6] * restored$R_2055 + biis[5] * restored$R_2060
  natural$N_2090 <- natural$N_2090 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[7] * restored$R_2030 + biis[7] * restored$R_2035 +
    biis[7] * restored$R_2040 + biis[7] * restored$R_2045 + biis[7] * restored$R_2050 +
    biis[7] * restored$R_2055 + biis[7] * restored$R_2060
  natural$N_2100 <- natural$N_2100 + biis[7] * restored$R_2020 +
    biis[7] * restored$R_2025 + biis[7] * restored$R_2030 + biis[7] * restored$R_2035 +
    biis[7] * restored$R_2040 + biis[7] * restored$R_2045 + biis[7] * restored$R_2050 +
    biis[7] * restored$R_2055 + biis[7] * restored$R_2060
  
  if (max(maxValue(natural) > 1)) warning(paste("Returned raster includes values > 1. Layer maxima:", 
                                                paste(maxValue(natural), collapse=", "),"\n"))
  if (max(minValue(natural) < 0)) warning(paste("Returned raster includes values < 0. Layer minima:", 
                                                paste(minValue(natural), collapse=", "),"\n"))

  if (plots==TRUE){
    par(mfrow=c(1,2))
    totals <- as.numeric(cellStats(natural, sum))
    plot(totals ~ v_years, ylab="Effective natural fraction", xlab="Year", cex.lab=1.3)
    bluered.plot(raster1 = subset(natural, which(v_years==2050)),
                 raster2 = subset(natural, which(v_years==2015)), 
                 main="Change in effective habitat fraction, 2015-2050")
    par(mfrow=c(1,1))
  }
  return(natural)
}


get.annual.area <- function(primary, secondary, years, cell.area, bii_sec){
  ## Function to interpolate annual effective area of natural habitat from time series of primary and secondary
  # stopifnot specifies required input formats; returns a 3-D array with same x and y resolution as inputs but one layer per year
  stopifnot(class(primary)=="RasterStack",
            class(secondary)=="RasterStack",
            class(years)=="numeric",
            class(cell.area)=="RasterLayer",
            class(bii_sec)=="numeric",
            dim(primary)==dim(secondary),
            length(years)== nlayers(primary),
            all.equal(years, sort(years)),
            length(unique(years))==length(years),
            dim(cell.land.area)[1:2]==dim(primary)[1:2],
            length(bii_sec)==1)
  area <- overlay(primary, secondary, cell.area, fun=function(x,y,z){return((x+bii_sec*y)*z)}) #SLOW if there are many years
  area_array <- as.array(area)
  
  annual_a <- array(NA, dim=c(dim(area)[1], dim(area)[2], max(years)-min(years)+1))
  for (i in 1:nrow(annual_a)){
    for (j in 1:ncol(annual_a)){
      if (is.na(sum(area_array[i, j, ]))){
        annual_a[i, j, ] <- NA
      }else{
        annual_a[i, j, ] <- approx(y=area_array[i, j, ], x=years, xout = c(min(years):max(years)))$y
      }
    }
  }
  return(annual_a)
}


annual_scenario_lu <- function(hist=historic_annual_area, pri, sec, mf, sf, historical_span=historical_span,
                               scenario_years=v_years, restored_years=v_res_years, bii_sec=bii_sec, 
                               bii_res=bii_res, cell_area_raster=cell_area_raster, overlap.plots=TRUE){
  
  # Function to produce annual effective area of primary habitat for scenarios
  
  natural <- stack(pri) + bii_sec * stack(sec)
  cat(paste("Maximum fraction of cell that is natural =",max(cellStats(natural, "max")),"\n"))
  restored <- stack(mf) + stack(sf)
  cat(paste("Maximum fraction of cell that is restored =",max(cellStats(restored, "max")),"\n"))
  fraction <- clamp(add_res_fraction(natural=natural, restored=restored, biis=bii_res), 
                    lower=0, upper=1, useValues=TRUE)
  area <- fraction * cell_area_raster
  annual_area <- interpolateTemporal(s=area, xin=scenario_years,
                                     xout=c(min(scenario_years):max(scenario_years)), 
                                     outdir="from_interpolateTemporal", prefix="annual_area", progress=TRUE, writechange=FALSE, 
                                     returnstack=TRUE, overwrite=TRUE) 
  annual_area_array <- as.array(annual_area)
  how_many_luh2_years <- 1985-historical_span[1]
  pre_scenario <- historic_annual_area[,,c(1:how_many_luh2_years)]
  annual_lu <- abind(pre_scenario, annual_area_array, along=3)
  
  if (overlap.plots==TRUE){
    alu <- apply(annual_lu, MARGIN=3, sum, na.rm=TRUE)
    plot(alu[1131:1151] ~ each_year[1131:1151] , type="l", lwd=3, col="lightgrey", xlab="Year", ylab="Area", main="LUH2: blue, Vivid: red, Returned: grey")
    vlu <- apply(annual_area_array, MARGIN=3, sum, na.rm=TRUE)
    points(vlu[1:16] ~ c(1985:2000), type="l", col="red")
    points(historic_annual_area[1131:1151] ~ each_year[1131:1151], type="l", col="blue" )
  }
  
  #Set attributes to help traceability
  attr(annual_lu, which="primary") <- pri
  attr(annual_lu, which="secondary") <- sec
  attr(annual_lu, which="mf") <- mf
  attr(annual_lu, which="sf") <- sf
  attr(annual_lu, which="historical_span") <- historical_span
  attr(annual_lu, which="scenario_years") <- v_years
  attr(annual_lu, which="restored_years") <- v_res_years
  attr(annual_lu, which="bii_sec") <- bii_sec
  attr(annual_lu, which="bii_res") <- bii_res
  attr(annual_lu, which="cell_area_raster") <- attr(cell_area_raster, which="file")@name
  attr(annual_lu, which="meld_year") <- 1985
  attr(annual_lu, which="box") <- NULL
  attr(annual_lu, which="when_lu_made") <- Sys.time()
  return(annual_lu)
}



## Function for linear interpolation into a raster stack
#### Temporal Interpolation ####################################################
# This function comes from https://gist.github.com/johnbaums/10465462
# Perform cell-wise linear interpolation between multiple raster layers, and 
# extrapolation beyond the upper limit of input data. Output is saved in .tif 
# format.
#
# Arguments
# s: a rasterStack containing the time slices to be interpolated 
#
# xin: a numeric vector that indicates the times associated with layers in s (in
# the same order as the layers of s - see names(s))
#
# xout: a numeric vector that indicates the times to be interpolated to (NB: if
# xout extends beyond the latest time slice in xin, it will be extrapolated
# using the rate of change from the period between the last and second to last
# time in xin.) 
#
# outdir: the directory to which files will be written (recursively created if 
# not already in existence) (character string)
#
# prefix: the output files will have pattern prefix_x.tif, where x is the
# timestep (potentially multiple digits), and prefix is a string that you 
# specify here (character string) 
#
# progress: show progress bar (TRUE/FALSE) 
#
# writechange: write the change grids that define the change in cell value per
# timestep between each pair of time slices, (TRUE/FALSE). If TRUE, these are
# written to outdir.
#
# returnstack: if TRUE, returns the interpolated grids (at timesteps xout) as a 
# rasterStack (TRUE/FALSE)
#
# ...: additional arguments passed to writeRaster

interpolateTemporal <- function(s, xin, xout, outdir, prefix, progress=FALSE, 
                                writechange=TRUE, returnstack=FALSE, ...) {
  require(raster)
  require(rgdal)
  if(missing(outdir)) stop('Please specify outdir')
  if(missing(prefix)) stop('Please specify prefix')
  if(nlayers(s) != length(xin)) stop('Length of xin must equal the number of layers in s.')
  if(nlayers(s) < 2) stop('stack s must have at least 2 layers.')
  if(!all(findInterval(xout, range(xin), rightmost.closed=TRUE) == 1)) {
    if(any(xout < min(xin))) {
      stop('This function does not extrapolate backwards (i.e. below the earliest element in xin). All elements of xout must be greater that min(xin).')
    } else {      
      warning('Some values of xout require extrapolation beyond the range of xin.\nThe rate of change for extrapolation is assumed to be equal to that for the period between the last and second-to-last elements of xin (after temporal ordering).')
    }
  }
  outdir <- normalizePath(sub('/$|\\\\$', '', outdir), winslash='/', 
                          mustWork=FALSE)
  if(!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
  xout <- unique(xout)
  if(is.unsorted(xin)) {
    s <- s[[order(xin)]]
    xin <- sort(xin)
  }
  len <- diff(xin)
  base <- findInterval(xout, xin)
  lower <- unique(base[base < nlayers(s)])
  s.change <- stack(sapply(if(length(lower) > 0) lower else nlayers(s) - 1, 
                           function(x) {
                             # message(sprintf('Calculating change grid for %s to %s.', xin[x], xin[x+1]))
                             overlay(s[[x]], s[[x+1]], fun=function(x1, x2) (x2-x1)/len[x],
                                     filename=ifelse(writechange, 
                                                     file.path(outdir, sprintf('changegrid_%s_%s', xin[x], xin[x+1])), 
                                                     ''), recycle=FALSE, format='GTiff', ...)
                           }))
  
  multi <- xout - xin[base]
  chg.ind <- ifelse(base > nlayers(s.change), nlayers(s.change), base)
  #message('Calculating grids for years specified in xout...')
  #if(progress) pb <- txtProgressBar(0, length(xout), style=3)
  invisible(sapply(seq_along(xout), function(x) {
    out.rast <- if(xout[x] %in% xin) {
      s[[base[x]]]
    } else {
      overlay(s[[base[x]]], s.change[[chg.ind[x]]],
              fun=function(x1, x2) x1 + (x2*multi[x]))
    }
    writeRaster(out.rast, 
                filename=file.path(outdir, sprintf('%s_%s', prefix, xout[x])), 
                format='GTiff', ...)
    #if(progress) setTxtProgressBar(pb, x)
  }))
  if(isTRUE(returnstack)) {
    f <- file.path(outdir, paste0(prefix, '_', xout, '.tif'))
    return(stack(f[order(as.numeric(sub('.*_(\\d+)\\.tif$', '\\1', f)))]))
  }
}

crop <- function(array, box=bounding.box){
  # Simple function for extracting a bounding box from array of land use or anything else
  if (length(dim(array))==2){
    cropped <- array[box$x, box$y]
  }else{
    cropped <- array[box$x, box$y,]
  }
  dims <- dim(array)
  
  mostattributes(cropped)<-attributes(array)
  attr(cropped, which="box") <- bounding.box
  to.return <- array(data=cropped, dim=c(length(box$x), length(box$y), dim(array)[3]))
  return(to.return)
}

relax <- function(natarea, years, wanted, c_array=array(1, dim=c(nrow(natarea), ncol(natarea))), rsr_2015, k, z){
  ## Performs the relaxation of endemic species richness based on the trajectory of natural area within each grid cell
  ## natarea = array of natural area each year
  ## years = years (consecutive) for which natural area is provided
  ## wanted = years for which S_t must be estimated
  ## c_array = array holding estimates of c for each grid cell; defaul is all 1
  ## z = exponent of species-area relationship
  ## k = rate constant of relaxation process
  ## rsr_2015 = array of range-size rarity in 2015
  ## Returns an array of relative species richness in each grid cell for each year in wanted
  
  stopifnot(is.array(natarea),
            length(dim(natarea))==3,
            is.numeric(years),
            years==sort(years), 
            is.numeric(wanted),
            wanted==sort(wanted),
            dim(c_array)==c(nrow(natarea), ncol(natarea)),
            dim(rsr_2015)[1:2]==c(nrow(natarea), ncol(natarea)))
  if (length(dim(rsr_2015))==2) dim(rsr_2015)<-c(dim(rsr_2015),1) # Varies depending on where it comes from
  
  year_0 <- years[1]
  
  S.0 <- array(data = c_array[,] * natarea[,,1] ^ z, dim=c(nrow(natarea), ncol(natarea))) #SAR
  # print(mean(S.0, na.rm=TRUE)) #Just a useful check that nothing's gone wrong
  
  S_t <- array(NA, dim=c(nrow(natarea), ncol(natarea), length(wanted))) #To hold species numbers remaining at time t
  for (i in 1:nrow(natarea)){
    for (j in 1:ncol(natarea)){
       if (is.na(sum(natarea[i, j, ]))){
        S_t[i, j, ] <- NA
      }else{
        S_t[i, j, 1] <- min(S.0[i,j],
                            interpolated.species.at.t(t=wanted[1]-year_0, S0=S.0[i,j], k=k, c=c_array[i,j], z=z, Area=natarea[i,j,]))
        for (q in 2:length(wanted)){
          S_t[i, j, q] <- min(S_t[i, j, q-1], 
                              interpolated.species.at.t(t=wanted[q]-year_0, S0=S.0[i,j], k=k, c=c_array[i,j], z=z, Area=natarea[i,j,]))
          #The min ensures there is no resurrection: a cell's endemic diversity cannot increase over time whatever happens to area
        }
      }
    }
  }
  
  if (all(c_array==1) & !any(is.na(c_array==1))){
    # Didn't know c so need to now rescale
    S_t_rescaled <- array(NA, dim=dim(S_t))
    for (i in 1:nrow(S_t)){
      for (j in 1:ncol(S_t)){
        S_t_rescaled[i,j,] <- S_t[i, j,]*rsr_2015[i,j,1]
      }
    }
    c_array <- S_t_rescaled[,,1]/natarea[,,1]^z #Calculate c_array
    image(t(c_array), main="c_array")
    
  }else{
    # c was provided so no need to rescale
    S_t_rescaled <- S_t
  }
  # Set attributes to help with traceability
  attr(S_t_rescaled, which = "when_lu_made") <- attr(natarea, which="when_lu_made")
  attr(S_t_rescaled, which = "lu_years_range") <- c(min(years), max(years))
  attr(S_t_rescaled, which = "St_years") <- wanted
  attr(S_t_rescaled, which = "k") <- k
  attr(S_t_rescaled, which = "z") <- z
  
  return(list("S_t" = S_t_rescaled, "c_array" = c_array))
}

percentage_plot <- function(A, S, A0=NULL, S0=NULL, year, main_stem){
  # Plots time series of effective primary habitat and endemic species remaining
  # as percentages of A0 and S0 respectively
  
  A_t <- apply(A, MARGIN=3, sum, na.rm=TRUE)
  S_t <- apply(S, MARGIN=3, sum, na.rm=TRUE)
  if (is.null(A0)) A0 <- A_t[1]
  if (is.null(S0)) S0 <- S_t[1]
  plot(100*A_t/A0 ~ year, type="l", lwd=2, col="green", ylim=c(0,100), 
       ylab="Percentage remaining", xlab="Year", cex.lab=1.3,
       main=paste(main_stem, "(habitat: green; endemic species: blue)"))
  points(100*S_t/S0 ~ year, type="l", lwd=2, col="blue")
}


extinctions_map <- function(S, year, between=c(2020,2050), main_stem){
  # Simple map of extinctions between the years between[1] and between[2]
  diff <- S[,,which(year==between[1])] - S[,,which(year==between[2])]
  image(t(diff), ylim=c(1,0), main=paste(main_stem, "extinctions between", between[1], "and", between[2]))
}


get_totals <- function(S, S0=NULL, year){
  # Calculate extinction totals for a scenario from the vector of S values that it produces
  # S0 can be provided, if the span of years starts later than historical_span[1]
  
  if (is.null(S0)) S0 <- S[1]
  totals <- list(original_richness=NA, before_2020=NA, from_2020_to_2050=NA, from_2020_to_2100=NA)
  totals$original_richness <- S0
  totals$before_2020 <- S0 - S[which(year==2020)]
  totals$from_2020_to_2050 <- S[which(year==2020)] - S[which(year==2050)]
  totals$from_2020_to_2100 <- S[which(year==2020)] - S[which(year==2100)]
  return(totals)
}

plot4 <- function(pri, sec, mf, sf, layer=13){
  # Plots the cell fractions, for the specified layer, of primary, secondary, mf and sf
  # Makes no attempt to consider WHEN recovery or loss took place
  # A useful sense-check on input data files
  # layer = 13 is for 2050
  par(mfrow=c(2,2))
  plot(raster(pri, layer), main=pri)
  plot(raster(sec, layer), main=sec)
  plot(calc(stack(mf), sum), main=mf)
  plot(calc(stack(sf), sum), main=sf)
  par(mfrow=c(1,1))
}

bluered.plot <- function(raster1, raster2, main=""){
  # Slightly kludgy but plots differences within the range -1 to +1, with zero being white
  mindiff <- round(cellStats(raster1-raster2, "min"), 1)
  maxdiff <- round(cellStats(raster1-raster2, "max"), 1)
  if (maxdiff - mindiff < 0.1){
    print("Differences are too small for plot to be meaningful")
  }else{
    mycols <- bluered(20)
    tooblue <- 10 + (mindiff * 10)
    toored <- (1 -maxdiff) * 10
    mycols <- mycols[1:(20-toored)]
    mycols <- mycols[-c(1:tooblue)]
    plot(raster1 - raster2, col=mycols, main=main)
  }
}

check_restoration <- function(pri, sec, mf, sf, layer=19, stem=""){
  yr <- v_years[layer]
  to_sum <- which(v_res_years <= v_years[layer]) # Only sum restored land to this point
  natural <- stack(pri) + bii_sec * stack(sec)
  cat(paste("Maximum fraction of cell that is natural =",max(cellStats(natural, "max")),"\n"))
  restored <- stack(mf) + stack(sf)
  cat(paste("Maximum fraction of cell that is restored =",max(cellStats(restored, "max")),"\n"))
  fraction <- clamp(add_res_fraction(natural=natural, restored=restored, biis=bii_res), 
                    lower=0, upper=1, useValues=TRUE)
  par(mfrow=c(2,2))
  plot(subset(natural, layer), main=paste(stem,"pri + bii_sec in", yr))
  plot(calc(subset(restored, to_sum), sum), main=paste(stem,"mf + sf by", yr))
  plot(subset(fraction, layer), main=paste(stem,"Effective natural fraction in", yr))
  plot(subset(fraction, layer) - subset(natural, layer), 
       main=paste(stem,"Fraction added by restoration,", yr))
  par(mfrow=c(1,1))
  effective_area_restored <- cell_land_area * (subset(fraction, layer) - subset(natural, layer))
  cat(paste("Pri + BIIsec in", yr, ":", 
            cellStats(subset(natural, layer) * cell_land_area, "sum"), "\n"))
  cat(paste("Area set aside for restoration:", 
            cellStats(calc(subset(restored, to_sum), sum) * cell_land_area, "sum"), "\n"))
  cat(paste("Effective area contributed by restoration as of", yr,
            cellStats(effective_area_restored, "sum"), "\n"))
  cat(paste("Average BII of restored land =", 
            cellStats(effective_area_restored, "sum")/cellStats(calc(subset(restored, to_sum), sum) * cell_land_area, "sum")), "\n")
  return(effective_area_restored)
  }

plot.vivid.pri.sec <- function(df, years=v_years){
  lims <- c(min(df), max(df))
  plot(base ~ years, data=df, ylim=lims, type="l", col="black", ylab="Global total area", xlab="Year", cex.lab=1.3)
  points(early ~ years, data=df, type="l", col="green")
  points(early125 ~ years, data=df, type="l", col="darkgreen")
  points(late29 ~ years, data=df, type="l", col="brown")
}

plot.vivid.mf.sf <- function(df, years=v_res_years){
  par(mfrow=c(1,2))
  lims <- c(min(df), max(df))
  plot(base ~ years, data=df, ylim=lims, type="l", col="black", ylab="Global total area converted in 5-year bin", xlab="Year", cex.lab=1.3)
  points(early ~ years, data=df, type="l", col="green")
  points(early125 ~ years, data=df, type="l", col="darkgreen")
  points(late29 ~ years, data=df, type="l", col="brown")
  lims <- c(min(colSums(df)), max(colSums(df)))
  plot(cumsum(base) ~ years, data=df, ylim=lims, type="l", col="black", ylab="Cumulative global total area", xlab="Year", cex.lab=1.3)
  points(cumsum(early) ~ years, data=df, type="l", col="green")
  points(cumsum(early125) ~ years, data=df, type="l", col="darkgreen")
  points(cumsum(late29) ~ years, data=df, type="l", col="brown")
  par(mfrow=c(1,1))
}

compare.scenarios <- function(year=2015, base, early, late){
  layer <- which(v_years==year)
  bf <- raster(base, layer)
  ef <- raster(early, layer)
  lf <- raster(late, layer)
  if(maxValue(ef-lf) > 0.05 | minValue(ef-lf) < -0.05){
    bluered.plot(raster1=ef, raster2=lf, main=paste(year, ": ", early, " - ", late, sep=""))
  }else{
    cat(paste("Early - Late: maximum = ", maxValue(ef-lf),", minimum = ", minValue(ef-lf), "\n", sep=""))
  }
  if(maxValue(ef-bf) > 0.05 | minValue(ef-bf) < -0.05){
    bluered.plot(raster1=ef, raster2=bf, main=paste(year, ": ", early, " - ", base, sep=""))
  }else{
    cat(paste("Early - Base: maximum = ", maxValue(ef-bf),", minimum = ", minValue(ef-bf), "\n", sep=""))
  }
}