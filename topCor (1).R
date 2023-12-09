topCor=function(img, slopeRast, aspectRast, azimuth, elevation, method="cor", includeReg=FALSE){
  # Source: https://web.pdx.edu/~nauna/Topographic_Normalization.pdf 
  require(terra, pracma)
  # convert all angles into radians
  illuGeom=c(deg2rad(azimuth),deg2rad(90-elevation)) #azimuth, zenith
  aspectRad=lapp(aspect, fun=function(x) deg2rad(x))
  slopeRad=lapp(slope, fun=function(x) deg2rad(x))
  
  #produce illumination model
  illuRast=cos(illuGeom[2])*cos(slopeRad)+sin(illuGeom[2])*sin(slopeRad)*cos(illuGeom[1]-aspectRad)
  #plot(illuRast)
  if(method=="illu"){
    return(illuRast)
    } 
  
  if (method=="cosine"){
    imgCor=img*cos(illuGeom[1])/illuRast
    return(imgCor)
  }
  
  if (method=="cosineMod"){
    illuMean=mean(values(illuRast), na.rm=TRUE)
    imgCor=img+((illuMean-illuRast)/illuMean)
    return(imgCor)
  }
  
  if(method=="C"){
    regression_summaries <- list() # stores the regression results for each band
    predicted_raster <- img # stores teh predicted values
    # Obtain regression summmaries for each band
    for (band_num in 1:nlyr(img)) {
      band <- subset(img, band_num)
        
      # Create a data frame with band and illumination values
      data_df <- data.frame(Band = values(band), illumVal = values(illuRast))
      names(data_df)=c("Band", "Illumination")
        
      # Fit the linear regression model for the band
      model <- lm(Band ~ Illumination, data = data_df)
      
      # Calculate the C parameter
      cval=model$coefficients[1]/model$coefficients[2]
      
      # Produce correction
      predicted_values <- values(predicted_raster[[i]])*(cos(illuGeom[1])+cval)/
        (values(illuRast)+cval)
      
      # Apply the linear regression model to predict the band values
      #predicted_values <- model$coefficients[1] + model$coefficients[2]*values(illuRast)
      
      # Set the predicted values in the corresponding band of the output raster
      predicted_raster[[band_num]]<- terra::setValues(predicted_raster[[band_num]],
                                                     predicted_values)
        
      # Print the regression summary for the current band
      cat("Prediction for Band", band_num, ":\n")
      #print(summary(model))
      
      # Store the summary of the regression model
      regression_summaries[[paste0("Band", band_num)]] <- summary(model)
    }
    #imgCor=img-predicted_raster
    #return(list(imgCor, regression_summaries))
    if(includeReg==TRUE) {
      return(list(predicted_raster, regression_summaries))
    } else {
      return(predicted_raster)
    }
  }
  
  if(method=="cor"){
    regression_summaries <- list() # stores the regression results for each band
    predicted_raster <- img # stores teh predicted values
    # Obtain regression summmaries for each band
    for (band_num in 1:nlyr(img)) {
      band <- subset(img, band_num)
      
      # Create a data frame with band and illumination values
      data_df <- data.frame(Band = values(band), illumVal = values(illuRast))
      names(data_df)=c("Band", "Illumination")
      
      # Fit the linear regression model for the band
      model <- lm(Band ~ Illumination, data = data_df)
      
      # Apply the linear regression model to predict the band values
      predicted_values <- model$coefficients[1] + model$coefficients[2]*values(illuRast)
      
      # Set the predicted values in the corresponding band of the output raster
      predicted_raster[[band_num]]<- terra::setValues(predicted_raster[[band_num]],
                                                      predicted_values)
      
      # Print the regression summary for the current band
      cat("Prediction for Band", band_num, ":\n")
      #print(summary(model))
      
      # Store the summary of the regression model
      regression_summaries[[paste0("Band", band_num)]] <- summary(model)
    }
    #slopeRadNorm=(slopeRad-min(values(slopeRad, na.rm=TRUE)))/(max(values(slopeRad, na.rm=TRUE))-min(values(slopeRad, na.rm=TRUE)))
    imgCor=img-(predicted_raster)#*slopeRadNorm)
    if(includeReg==TRUE) {
      return(list(imgCor, regression_summaries))
    } else {
      return(imgCor)
    }
  }
  
  if (method=="minnaert"){
    regression_summaries <- list() # stores the regression results for each band
    predicted_raster <- img # stores teh predicted values
    # Obtain regression summaries for each band
    for (band_num in 1:nlyr(img)) {
      band <- subset(img, band_num)
      
      # Create a data frame with band and illumination values
      data_df <- data.frame(Band = values(band), illumVal = values(illuRast))
      names(data_df)=c("Band", "Illumination")
      
      # Fit the linear regression model for the band
      #data_df2=data.frame(cbind(log(data_df$Band), log(data_df$Illumination/cos(illuGeom[1]))))
      model <- lm(log(Band) ~ log(Illumination/cos(illuGeom[1])), data = data_df)
      corrected_values=Band*cos(slopeRad)((cos(illuGeom[1])/(illuRast*cos(slopeRad)))^model$coefficients[2])
      
      predicted_raster[[band_num]]<- terra::setValues(predicted_raster[[band_num]],
                                                      corrected_values)
      # Store the summary of the regression model
      regression_summaries[[paste0("Band", band_num)]] <- summary(model)
    }
    return(list(predicted_raster, regression_summaries))
  }
}
