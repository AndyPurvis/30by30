# Function to convert IAM data to effective natural area on a 2-degree grid 
# Andy Purvis
# 2021-02-11

convert_iam_data <- function(iam_path, which_bii, ena_path, pdf_path, yb=2010, ye=2100, ys=5, layers=1){
  
  biplot <- function(brick, LU){
    plot(brick, 1, ylim=c(-90,90), axes=FALSE, main = paste(LU, "in 2010"), sub = iam_path)
    plot(seq(yb,ye,ys), as.numeric(cellStats(brick, sum)), ylab = "Total area (Mha)", xlab="Year", main=LU, sub=which_bii)
  }
  
  require(raster)
  require(scales)
  require(abind)
  require(gplots)
  require(colorspace)
  require(ncdf4) # package for netcdf manipulation, suggested by https://rpubs.com/boyerag/297592
  require(rgdal) # package for geospatial analysis, suggested by https://rpubs.com/boyerag/297592
  require(ggplot2) # package for plotting, suggested by https://rpubs.com/boyerag/297592
  
  stopifnot(which_bii %in% c("BII_crop0", "BII_humdom0")) #These are the only permitted options
  
  # How many layers?
  nl <- nlayers(brick(iam_path, varname = "LC_area_share"))
  
  # Read in cell areas
  cell_areas <- raster(iam_path, varname = "pixel_area")
  
  # turn warnings off
  options(warn=-1)
  
  # Multiply cell land-use fractions by cell area before summing into 2-degree grids
  # Note that, by default, aggregate has na.rm=TRUE, which is what I want as land-free cells are NA
  crop_other <- aggregate(brick(iam_path, 
                                varname="LC_area_share", lvar = 3, nl=nl, level=1) * cell_areas, 
                          fact = 4, fun = sum)
  crop_2Gbioen <- aggregate(brick(iam_path, 
                                  varname="LC_area_share", lvar = 3, nl=nl, level=2) * cell_areas, 
                            fact = 4, fun = sum)
  grassland <- aggregate(brick(iam_path, 
                               varname="LC_area_share", lvar = 3, nl=nl, level=3) * cell_areas, 
                         fact = 4, fun = sum)
  forest_unmanaged <- aggregate(brick(iam_path, 
                                      varname="LC_area_share", lvar = 3, nl=nl, level=4) * cell_areas, 
                                fact = 4, fun = sum)
  forest_managed <- aggregate(brick(iam_path, 
                                    varname="LC_area_share", lvar = 3, nl=nl, level=5) * cell_areas, 
                              fact = 4, fun = sum)
  restored <- aggregate(brick(iam_path, 
                              varname="LC_area_share", lvar = 3, nl=nl, level=6) * cell_areas, 
                        fact = 4, fun = sum)
  other <- aggregate(brick(iam_path, 
                           varname="LC_area_share", lvar = 3, nl=nl, level=7) * cell_areas, 
                     fact = 4, fun = sum)
  built_up <- aggregate(brick(iam_path, 
                              varname="LC_area_share", lvar = 3, nl=nl, level=8) * cell_areas, 
                        fact = 4, fun = sum)
  abn_cropland_other <- aggregate(brick(iam_path, 
                                        varname="LC_area_share", lvar = 3, nl=nl, level=9) * cell_areas, 
                                  fact = 4, fun = sum)
  abn_cropland_2Gbioen <- aggregate(brick(iam_path, 
                                          varname="LC_area_share", lvar = 3, nl=nl, level=10) * cell_areas, 
                                    fact = 4, fun = sum)
  abn_grassland <- aggregate(brick(iam_path, 
                                   varname="LC_area_share", lvar = 3, nl=nl, level=11) * cell_areas, 
                             fact = 4, fun = sum)
  abn_forest_managed <- aggregate(brick(iam_path, 
                                        varname="LC_area_share", lvar = 3, nl=nl, level=12) * cell_areas, 
                                  fact = 4, fun = sum)
  
  # Read in file of BII values
  bii_coefficients <- read.csv("output/Rescaled_BII_values_for_IAM_classes.csv")
  
  if (which_bii == "BII_humdom0"){
    bii <- bii_coefficients$BII_humdom0
  }else{
    bii <- bii_coefficients$BII_crop0
  }
  
  # Multiply each raster by its BII and add the products together
  effective_natural_area <- crop_other * bii[1] +
    crop_2Gbioen * bii[2] +
    grassland * bii[3] +
    forest_unmanaged * bii[4] +
    forest_managed * bii[5] +
    restored * bii[6] +
    other * bii[7] +
    built_up * bii[8] +
    abn_cropland_other * bii[9] +
    abn_cropland_2Gbioen * bii[10] +
    abn_grassland * bii[11] +
    abn_forest_managed * bii[12]
  
  # turn warnings back on
  options(warn=0)
  
  # Write to file, including the choice of BII values in the file name
  if (require(ncdf4)) {	
    rnc <- writeRaster(effective_natural_area, filename=ena_path, format="CDF", overwrite=TRUE) 
  }
  
  if(!is.null(pdf_path)){
    pdf(pdf_path, height=10, width=8)
    par(mfrow=c(3,2))
    biplot(crop_other, "crop_other")
    biplot(crop_2Gbioen, "crop_2Gbioen")
    biplot(grassland, "grassland")
    biplot(forest_unmanaged, "forest_unmanaged")
    biplot(forest_managed, "forest_managed")
    biplot(restored, "restored")
    biplot(other, "other")
    biplot(built_up, "built_up")
    biplot(abn_cropland_other, "abn_cropland_other")
    biplot(abn_cropland_2Gbioen, "abn_cropland_2Gbioen")
    biplot(abn_grassland, "abn_grassland")
    biplot(abn_forest_managed, "abn_forest_managed")
    biplot(effective_natural_area, "Natural")
    par(mfrow=c(1,1))
    dev.off()
  }
return(subset(effective_natural_area, layers)) # Just return the requested rasters
}
