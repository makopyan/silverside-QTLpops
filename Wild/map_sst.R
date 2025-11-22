library(ncdf4)
library(httr)
library(tidyverse)
library(reshape2)


junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v1_0_monthly.nc?analysed_sst[(2017-05-01T12:00:00Z):1:(2018-12-01T12:00:00Z)][(25):1:(50)][(275):1:(300)]', write_disk("sst.nc", overwrite=TRUE))

nc=nc_open('sst.nc')

names(nc$var)
v1=nc$var[[1]]
sst=ncvar_get(nc,v1)
dim(sst)
dates=as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') 
dates
lon=v1$dim[[1]]$vals 
lat=v1$dim[[2]]$vals

nc_close(nc) 
rm(junk,v1) 
file.remove('sst.nc')


h=hist(sst[,,1], 100, plot=FALSE) 
breaks=h$breaks 
n=length(breaks)-1


jet.colors <-colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))


c=jet.colors(n)


# Convert your SST array to a dataframe for ggplot
sst_df <- expand.grid(lon = lon, lat = lat)
sst_df$temp <- as.vector(sst[,,1])

ggplot() +
  geom_tile(data = sst_df, aes(x = lon, y = lat, fill = temp)) +
  scale_fill_gradientn(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0",
                                  "#FDDBC7","#F4A582","#D6604D","#B2182B"),
                       name = "SST", na.value = "white") +
  geom_point(data = pop, aes(x = long, y = lat), color = 'black', size = 2) +
  geom_text(data = pop, aes(x = long, y = lat, label = pop), 
            color = 'black', size = 3, vjust = -0.5) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
