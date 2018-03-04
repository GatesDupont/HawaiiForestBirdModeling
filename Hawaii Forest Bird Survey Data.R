#----Libraries----
if(T){
  library(sp)
  library(leaflet)
  library(viridis)
  library(RColorBrewer)
  library(classInt)
  library(ggplot2)
  library(nlme)
  library(mgcv)
  library(gamm4)
  library(gstat)
  library(raster)
}

#----Read in the data----
if(T){
  df = read.csv("~/Desktop/HawaiiForestBirdSurveyData/HFBS.csv")
}

#----Re-projecting data; hi4/5;----
if(T){
  utm4 = df[df$Projection == "UTM Zone 4",]
  utm5 = df[df$Projection == "UTM Zone 5",]
  
  coords4 = cbind(utm4$EASTING, utm4$NORTHING)
  spdf4 = SpatialPointsDataFrame(coords4, utm4, proj4string=CRS("+proj=utm +zone=4 +datum=NAD83"))
  hi4 = spTransform(spdf4, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  coords5 = cbind(utm5$EASTING, utm5$NORTHING)
  spdf5 = SpatialPointsDataFrame(coords5, utm5, proj4string=CRS("+proj=utm +zone=5 +datum=NAD83"))
  hi5 = spTransform(spdf5, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
}

#----Leaflet of data by UTM----
if(F){
  leaflet() %>% addTiles() %>%
    addCircleMarkers(lng = hi4@coords[,1], lat = hi4@coords[,2], color = "blue") %>%
    addCircleMarkers(lng = hi5@coords[,1], lat = hi5@coords[,2], color = "red")
}

#----Back to normal; df = hi----
if(T){
  df4 = as.data.frame(hi4)
  df5 = as.data.frame(hi5)
  
  hi = rbind(df4,df5)
  
  hi = data.frame(hi$coords.x1, hi$coords.x2, hi$SpeciesName, hi$NumberBirdsDetected,
                  hi$Count, hi$SURVEY.AREA, hi$Year, hi$MONTH.., hi$OBSVDATE,
                  hi$TRANSECT, hi$STATION, hi$OBSCODE, hi$Origin, hi$Habitat)
  colnames(hi) = c('lng', 'lat', 'species', 'nseen', 
                   'countPeriod', 'survey', 'year', 'month', 'date',
                   'transect', 'station', 'observer', 'origin', 'habitat')
}

#----Species selection, sp.x----
if(T){
  sp.x = hi[hi$species=="Wild Turkey",]
}

#----Survey selection----
if(T){
  sp.x = sp.x[sp.x$survey == "Hamakua HFB Survey"|sp.x$survey == "Kau HFB Survey"|sp.x$survey == "Kipukas HFB Survey"|sp.x$survey == "Kohala HFB Survey"|sp.x$survey == "Kona HFB Survey"|sp.x$survey == "Mauna Kea HFB Survey"|sp.x$survey == "Puna HFB Survey",]
  #length(unique(sp.x$survey))
}

#----Plotting species counts----
if(T){
  myColors <- rev(heat.colors(12, alpha = 1))
  ggplot(sp.x, aes(lng, lat, col=factor(nseen))) + 
    geom_point() + scale_color_manual(values=myColors) +
    theme(panel.background = element_rect(fill = 'gray'),
          panel.grid.major = element_line(colour = 'gray'),
          panel.grid.minor = element_line(colour = 'gray'),
          plot.title = element_text(face = 'bold', hjust = 0.5, family = 'sans'))
}

#----Plotting by survey----
if(T){
  ggplot(sp.x, aes(lng, lat, col=factor(survey))) + 
    geom_point() + scale_color_discrete(name="Month") +
    theme(panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'white'),
          panel.grid.minor = element_line(colour = 'white'),
          plot.title = element_text(face = 'bold', hjust = 0.5, family = 'sans'))
}

#----GAM color pallette plot----
if(F){
  pal = brewer.pal(5, "Reds")
  q5 = classIntervals(sp.x$nseen, n=5, style = "quantile")
  q5Colours = findColours(q5,pal)
  plot(sp.x$lng, sp.x$lat, col = q5Colours, pch = 19, axes = T, cex = 0.3, main = "maxFlock")
  legend("topleft", fill = attr(q5Colours, "palette"),
         legend = names(attr(q5Colours, "table")),bty = "n")
}

#----First GAM----
if(F){
  sp.gam = gam(nseen~s(lat,lng),data = sp.x)
  summary(sp.gam)
}

#----Deviance smoothing----
if(F){
  dev.rss = c()
  kval = c()
  for(i in seq(120)){
    dev.rss = c(dev.rss, deviance(gam(nseen~s(lat,lng, k = i),data = sp.x)))
    kval = c(kval, i)
  }
  plot(kval, dev.rss, xlab = "Parameters", ylab = "Deviance (RSS)", 
       pch=15, main = "Smoothing Parameter Guide")
}

#----AIC smoothing----
if(F){
  aic.k = c()
  kval = c()
  for(i in seq(200)){
    aic.k = c(aic.k, AIC(gam(nseen~s(lat,lng, k = i),data = sp.x)))
    kval = c(kval, i)
  }
  plot(kval, aic.k, xlab = "Parameters", ylab = "AIC", pch=15, 
       main = "Smoothing Parameter Guide")
}

#----GAM k=70----
if(T){
  #sp.gam.xy = gam(nseen~s(lat,lng,k=70),data = sp.x)
  sp.gam.xy = gam(nseen~s(lat,lng,k=70)+month+year+observer,data = sp.x) 
  sp.gam.xy.pred = predict(sp.gam.xy,se.fit=T)
  summary(sp.gam.xy)
}

#----GAM k=70 predictions df----
if(T){
  sp.gam.xy.4pred = data.frame(
    x = sp.x$lng,
    y = sp.x$lat ,
    pred = fitted(sp.gam.xy))
  head(sp.gam.xy.4pred)
  coordinates(sp.gam.xy.4pred) = c("x","y")
}

#----GAM k=70 pred. plot----
if(T){
  pal = brewer.pal(5,"Reds")
  q5 = classIntervals(sp.gam.xy.4pred$pred, n=5, style = "quantile")
  q5Colours = findColours(q5, pal)
  plot(sp.gam.xy.4pred, col=q5Colours,pch=19,cex=0.7,axes=T,main="GAM k=120")
  legend("topleft", fill=attr(q5Colours, "palette"),
         legend = names(attr(q5Colours,"table")),cex=0.7,bty="n")
}

#----GAM model selection----
if(F){
  
  # Fundamental
  gam0 = gam(nseen~s(lat,lng,k=70),data = sp.x)
  
  # Month,Year
  gam1 = gam(nseen~s(lat,lng,k=70)+month+year,data = sp.x) # Runner up
  gam2 = gam(nseen~s(lat,lng,k=70)+year,data = sp.x)
  gam3 = gam(nseen~s(lat,lng,k=70)+month,data = sp.x) # Top model only by .2
  AIC(gam0,gam1,gam2,gam3)
  
  # Regular date
  gam4 = gam(nseen~s(lat,lng,k=70)+date,data = sp.x) # Better model but df??
  AIC(gam0,gam4)
  
  # Observer
  gam5 = gam(nseen~s(lat,lng,k=70)+observer,data = sp.x)
  #gam6 = gamm(nseen~s(lat,lng,k=70), random = list(observer=~1),data = sp.x) # Note random var format
  #gam7 = gamm(nseen~s(lat,lng,k=70),data = sp.x) # Note random var format
  AIC(gam0,gam5)
  #AIC(gam6)
  
  # Observer 2
  gam8 = gam(nseen~s(lat,lng,k=70)+month+year+observer,data = sp.x) 
  gam9 = gam(nseen~s(lat,lng,k=70)+month*year+observer,data = sp.x)
  AIC(gam8,gam9)
}

#----Variogram----
if(T){
  
  sp.x$logNseen = log(sp.x$nseen)
  
  Xloc = sp.x$lng
  Yloc = sp.x$lat
  
  sp.x.spdf <- SpatialPointsDataFrame(coords = sp.x[,c("lng", "lat")], data = sp.x,
                                             proj4string = CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  sp.x.spdf <- sp.x.spdf[-zerodist(sp.x.spdf)[,1],] 
  
  # Log transformed
  log.nseen.vario = variogram(logNseen~1,sp.x.spdf)
  plot(log.nseen.vario, pch=20,col=1,cex=2)
  logNseen.fit = fit.variogram(log.nseen.vario,
                            vgm(psill=0.35,"Sph",range=40,nugget=0.2))
  plot(log.nseen.vario,logNseen.fit,pch=20,col=2,cex=2,lwd=3,
       main=paste("Log nseen Variogram","- ", unique(sp.x$species)))
  
  # Original
  nseen.vario = variogram(nseen~1,sp.x.spdf)
  plot(nseen.vario, pch=20,col=1,cex=2)
  nseen.fit = fit.variogram(nseen.vario,
                               vgm(psill=1.4,"Sph",range=40,nugget=1.4)) #3, 2
  plot(nseen.vario,nseen.fit,pch=20,col=2,cex=2,lwd=3,
       main=paste("nseen Variogram","- ", unique(sp.x$species)))
}

#----Kriging----
if(T){
  # Set up polygon for prediction grid
  nseen.chull = chull(Xloc, Yloc)
  plot(Xloc, Yloc)
  lines(Xloc[nseen.chull],Yloc[nseen.chull])
  
  #Prediction grid
  library(geoR)
  nseen.grid = polygrid(
    xgrid=seq(min(Xloc),max(Xloc),length=20),
    ygrid=seq(min(Yloc),max(Yloc),length=20),
    cbind(
      Xloc[nseen.chull],
      Yloc[nseen.chull]))
  names(nseen.grid)=c("Xloc","Yloc")
  coordinates(nseen.grid)=c("Xloc","Yloc")
  nseen.grid = as(nseen.grid, "SpatialPixels")
  
  # Assure that data share projection
  proj4string(nseen.grid) = CRS(proj4string(sp.x.spdf))
  
  # Plot the sampling points with prediction grid
  plot(Yloc~Xloc,cex=1.2,pch=20,col=2)
  points(data.frame(nseen.grid)$Xloc,data.frame(nseen.grid)$Yloc,pch="+")	
  
  # Solve krige for predictions
  nseen.ok = krige(nseen~1, sp.x.spdf, nseen.grid, nseen.fit)
  
  # View range
  range(nseen.ok$var1.pred)
  
  # Produced kriged image
  image(nseen.ok["var1.pred"],col=heat.colors(10),axes=T)
  mtext(side=1,line=2,"Xloc")
  mtext(side=2,line=2,"Yloc")
  title(paste("Predicted", unique(sp.x$species), "Count - Big Island"))
  legend("topleft",legend=round(seq(min(nseen.ok$var1.pred),max(nseen.ok$var1.pred),length=10),2),fill=heat.colors(10),
         bty="n",title="Count")
}

#----Adding raster to leaflet----
if(T){
  # Fetching custom map tiles and adding citation
  custom_map = "https://api.mapbox.com/styles/v1/heliornis/cjboo3dac64w02srud7535eec/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1IjoiaGVsaW9ybmlzIiwiYSI6ImNqNGtjZjE1cjBoYmcycXAzbmlmNmJieGEifQ.ERz0DjQKEE1PBd7myLKwZA"
  mb_attribution <- "© <a href='https://www.mapbox.com/map-feedback/'>Mapbox</a> © <a href='https://www.sciencebase.gov/catalog/item/54c2db46e4b043905e0185cb'>Hawaii Forest Birds Survey Records</a> © <a href='https://www.gatesdupont.com/'>Gates Dupont</a>"
  
  r = raster(nseen.ok)
  pal = colorNumeric(c("#FFFF00", "#FFA500", "#FF0000"), values(r),
                     na.color = "transparent")
  
  leaflet() %>% addTiles(urlTemplate = custom_map, attribution = mb_attribution) %>%
    addRasterImage(r, colors = pal, opacity = 0.65) %>%
    addLegend(pal = pal, values = values(r),
              title = "Predicted count")
}
