## load packages (dependencies): 
library(CoordinateCleaner)
library(raster)
library(tidyverse)
library(rgdal)
library(maptools)
library(taxize)
library(rgbif)
library(data.table)
library(phytools)
library(sp)
library(phyloregion)
library(gdm)
library(terra)
library(Matrix)
library(dplyr)
library(canaper)
library(betapart)
library(sf)
library(graphics)
library(dunn.test)
library(multcompView)
library(igraph)
library(units)
library(fields)
library(betareg)
library(geosphere)
library(MuMIn)

# load in the occurrence dataset
load("FullyResolvedDist_Basic.RData")

# Read in the biome shapefile
# Reading this in with terra then converting it to an sf and making valid
# is the most effective for use in downstream sf functions
shp<-vect("Ecoregions2017/Ecoregions2017.shp")
shp<-st_as_sf(shp)
shp<-st_make_valid(shp)

# Read in Phylogenies
# Smith and Brown tree for seed plants
mega <- read.tree("Phylogeny/GBOTB.tre")
# Read in Nitta tree for Ferns
fern <- read.tree("Phylogeny/fern_megatree.tre")

# Read in rasters (predictor variables)
Rast1<-raster::brick("Rasters/Bioclim_5sec.tif") # 19 bioclim variables
elev<-raster("Rasters/elevation.tif")
fire<-raster("Rasters/FireFrequency.tif")
pstab<-raster("Rasters/precipStability.tif")
tstab<-raster("Rasters/tempStability.tif")

# Run phylogdm and taxgdm function scripts too

## Occurrence + biome intersection ----


# Run this in parallel. Even still it takes FOREVER

# Set up cluster

totalCores <- detectCores()
# Number of cores to run:
cluster <- makeCluster(totalCores[1]-38, outfile="biomes06142023")
registerDoParallel(cluster)
getDoParWorkers()
batch<-1:241926 #run the loop in batches to save incrementally and run faster
startTime <- Sys.time()

# run parallel loop to overlay the GBIF occurrence points with the biome shapefile
biomes <-  foreach(i =  batch, .combine = rbind) %dopar% {
  thisbatch<-Bif_PBasic[(((i-1)*1000)+1):(i*1000),]
  coords <- cbind(thisbatch$lon, thisbatch$lat)
  point_shp <- sp::SpatialPoints(coords,proj4string=sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
  pointsSPDF <- sp::SpatialPointsDataFrame(coords=point_shp,data=thisbatch)
  thisbatch_over<-sp::over(pointsSPDF, shp)
  biomes<-data.frame(cbind(thisbatch, thisbatch_over))
}
stopCluster(cluster)
endTime <- Sys.time()
print(endTime - startTime)


# For points near coastlines that did not intersect biome shapefile, get nearest ecoregion

# Which didn't intersect?
ecos_na<-biomes[which(is.na(biomes$ECO_NAME)),]
# Set up coordinates for input into for loop
coords <- st_as_sf(ecos_na, coords = c("lon", "lat"), crs = st_crs(shp))

nearPoly<-coords %>% 
  mutate(
    PolyNum = st_nearest_feature(geometry, shp)
  )

nearPoly_name<-c()
for(i in 1:nrow(nearPoly)){
  print(i)
  nearPoly_name[i]<-shp$ECO_NAME[nearPoly$PolyNum[i]]
}

# Okay great! Now merge the ecos_na df with the shp df
ecos_na$ECO_NAME<-nearPoly_name
head(ecos_na$ECO_NAME)
ecos_resolved<-ecos_na[,c(3:10,12)]
ecos_resolved<-merge(ecos_resolved, shp, by = "ECO_NAME")
ecos_resolved<-ecos_resolved[,-24] #remove "geometry" col
which(is.na(ecos_resolved$ECO_NAME)) #none!!!!! Yay!!!!!

# Remove the na occurrences from the biomes df
biomes<-biomes[-which(is.na(biomes$ECO_NAME)),-c(1:2)]
# Bind the resolved occs to the biomes df
biomes<-rbind(biomes, ecos_resolved)
write.csv(biomes, "biomes_full_2.csv")

## Remove tagged occurrences and seagrasses ----

# Remove any of the occurrences tagged with introduced, extinct, or location doubtful
biomes_clean<-biomes[-which(biomes$introduced > 0)]
biomes_clean<-biomes_clean[-which(biomes_clean$extinct > 0)]
biomes_clean<-biomes_clean[-which(biomes_clean$location_doubtful > 0)]
colnames(biomes_clean)
write.csv(biomes_clean, "biomes_clean_2.csv")
remove(biomes)

occs<-biomes_clean[,c("taxon_name", "lat", "lon", "BIOME_NAME")] # isolate just the columns we want
remove(biomes_clean) # we will just use the occs df from now on

# remove any "Rock and Ice" ecoregions (the biome is "N/A")
occs<-occs[-which(occs$BIOME_NAME == "N/A"),]

# re-name cols for input into points2comm in phyloregion
colnames(occs)<-c("species", "decimallatitude", "decimallongitude", "BIOME_NAME")

# remove those silly seagrasses
seagrass_genera<-c("Thalassodendron", "Amphibolis", "Oceana", "Cymodocea", 
                   "Syringodium", "Ruppia", "Posidonia", "Zostera", "Phyllospadix",
                   "Halodule", "Halophila", "Thalassia", "Enhalus")
occs<- occs %>% 
  filter(!sapply(strsplit(species, " "), `[`, 1) %in% seagrass_genera)

write.csv(occs, "occs_clean.csv", row.names = FALSE)
remove(seagrass_genera)

species<-unique(occs$taxon_name) # get a species list for updating phylogeny
write.csv(species, "Phylogeny/species.csv", row.names = FALSE)


## Bind together phylogenies and update phylogeny nomenclature ----

# Keep only angiosperms and gymnmosperms on Mega
mega <- extract.clade(mega, getMRCA(mega, c(which(mega$tip.label=="Ginkgo_biloba"),
                                            which(mega$tip.label=="Carex_grayi"))))
write.tree(mega, "Phylogeny/AngNGyms.tree")

# Remove angiosperm and gymnosperm tips from fern except for Sciadopitys_verticillata
fern <- drop.tip(fern, c(phangorn::Descendants(fern,
                                               getMRCA(fern, c(grep("Trithuria", fern$tip.label),
                                                               grep("Acorus", fern$tip.label))),
                                               "tips")[[1]],which(fern$tip.label=="Ginkgo_biloba")))

# Rename it in case that species occurs in both trees
fern$tip.label[which(fern$tip.label=="Sciadopitys_verticillata")] <- "temp_label"

# Attach it to the branch connecting the angiosperms to the gymnosperms on the fern tree
#at the position of the age of the angiosperm MRCA on the mega tree
combined_tree <- bind.tree(x = fern, y = mega, where = which(fern$tip.label=="temp_label"),
                           position = max(nodeHeights(mega)))
# Remove placeholder
combined_tree <- drop.tip(combined_tree, "temp_label")
write.tree(combined_tree, "Phylogeny/combined_tree.tree")
combined_tree<-read.tree("Phylogeny/combined_tree.tree")
grep("temp_label", combined_tree$tip.label) #none. yay!

remove(mega)
remove(fern)
str(combined_tree)
# There are 85,450 spp. in the tree


# How many are a direct match? 61,598 spp.
match<-data.frame(combined_tree$tip.label[which(combined_tree$tip.label %in% species$x)])

# These ones don't match:
unmatched<-data.frame(combined_tree$tip.label[-which(combined_tree$tip.label %in% species$x)]) #23852
length(combined_tree$tip.label)
unique(combined_tree$tip.label)

colnames(unmatched)[1]<-"search_names"

# Update unmatched tree names to match POWO names
current_name<-c()
accepted_name<-c()
unmatched$accepted<-NA
for(i in 1:nrow(unmatched)){
  print(i)
  current_name <- unmatched$search_names[i]
  print(current_name)
  # Use the `pow_search` function to search for the name in the POWO database
  # Limit the results to 1 for efficiency
  accepted_name <- pow_search(current_name, db = "wcsp", limit = 1)
  print(accepted_name)
  # If no data is returned (name not found in POWO), set the accepted name to NA
  if (is.null(accepted_name[["data"]])){
    unmatched$accepted[i] <- NA
  }
  # If the name is found and is an accepted name, update the 'accepted' column
  if (!is.null(accepted_name[["data"]]) && accepted_name[["data"]][["accepted"]] == TRUE) {
    unmatched$accepted[i]<- accepted_name[[2]][["name"]]
  }
  # If the name is found but is a synonym, update the 'accepted' column with its accepted name
  if (!is.null(accepted_name[["data"]]) && accepted_name[["data"]][["accepted"]] == FALSE){
    unmatched$accepted[i]<-accepted_name[["data"]][["synonymOf"]][["name"]]
  }
}

write.csv(unmatched, "Phylogeny/phy_resolved_names.csv", row.names = FALSE)
unmatched<-as.data.frame(fread("Phylogeny/phy_resolved_names.csv"))

# Remove any names in unmatched$accepted that are already in the tree to prevent duplicates
unmatched1<-unmatched[-which(unmatched$accepted %in% combined_tree$tip.label),]
# Okay that was a lot of them, now we hve 13,320
test <- unmatched1 %>% 
  group_by(accepted) %>% 
  count()

# Have to remove non-unique names, can't have duplicate names in phylogeny, 
# no way to know which is the one to use and which to drop
unmatched1<-unmatched1[which(unmatched1$accepted %in% test$accepted[which(test$n < 2)]),] 
SnB_updated<-combined_tree
SnB_updated$tip.label[1:5]
unmatched1$search_names[1:5]
for(i in 1:nrow(unmatched1)){
  print(i)
  SnB_updated$tip.label[which(unmatched1$search_names[i] == SnB_updated$tip.label)]<-unmatched1$accepted[i]
}

length(SnB_updated$tip.label)
unique(SnB_updated$tip.label)

which(SnB_updated$tip.label %in% species$x) # 68,923
write.tree(SnB_updated, "Phylogeny/SnB_updated.tre")
remove(combined_tree)

## Prepare rasters ----

# Align the extents of the rasters to those of the raster stack
elev<-resample(elev, Rast1)
fire<-resample(fire, Rast1)
pstab<-resample(pstab, Rast1)
tstab<-resample(tstab, Rast1)

# Join the layers into one raster
envRast<-addLayer(Rast1, elev, fire, pstab, tstab)
writeRaster(envRast, "envRast.tif", overwrite=TRUE)

remove(Rast1)
remove(elev)
remove(fire)
remove(pstab)
remove(tstab)


## Community grids + biome specific phylogenies + preliminary maps ----

# Make sure phylogeny tip labels match occs
SnB_updated$tip.label<-gsub("_", " ", SnB_updated$tip.label) 



###  "Tropical & Subtropical Dry Broadleaf Forests"

# Get biome-specific occurrence data for DBF
DBFOccs<-rbind(occs[which(occs$BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests"),])

# Only the columns we will need for points2comm() in phyloregion
DBFOccs_abbr<-DBFOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(DBFOccs_abbr, "DBF/DBFOccs_abbr.csv", row.names= FALSE)

# Great, now make the grid!
DBFGrid<-points2comm(dat = DBFOccs_abbr, res = 2, pol.grids = NULL) 
length(DBFGrid[["map"]]) #312
remove(DBFOccs, DBFOccs_abbr)

# We will need a df version of the matrix downstream
DBFMat<-sparse2dense(DBFGrid[["comm_dat"]])
str(DBFMat)
DBFMat<-as.data.frame(DBFMat)
rownames(DBFMat)[1:5]
colnames(DBFMat)[1:5]

#Make biome-specific phylogeny with only names in our new matrix
DBFPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(DBFMat)))==T)])
write.tree(DBFPhy, "DBF/DBFPhy.tree")

#Clean matrix to only include names that are in the phylogeny
DBFMat<-DBFMat[,which(names(DBFMat) %in% DBFPhy$tip.label)]
#Remove rows that no longer have data
DBFMat<-DBFMat[-which(rowSums(DBFMat)==0),]
#Convert matrix to presence-only
DBFMat[DBFMat>1]<-1
#Make a col for grid cells
DBFMat$grids<-row.names(DBFMat)
write.csv(DBFMat,"DBF/DBFMat.csv", row.names = FALSE)

#Save the map portion of DBFGrid with the appropriate CRS
DBFMap<-DBFGrid[["map"]][which(DBFGrid[["map"]]$grids %in% rownames(DBFMat)),]
length(unique(DBFMap$grids)) # 308
length(DBFMap$grids) #308
crs(DBFMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(DBFMap, "DBF/DBFMap.shp", overwrite = TRUE)
remove(DBFGrid)





### AWESOME! Now repeat this 14 times (13 more biomes + global): ###


### "Deserts & Xeric Shrublands"

DesOccs<-rbind(occs[which(occs$BIOME_NAME == "Deserts & Xeric Shrublands"),])

DesOccs_abbr<-DesOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(DesOccs_abbr, "Des/DesOccs_abbr.csv", row.names= FALSE)

DesGrid<-points2comm(dat = DesOccs_abbr, res = 2, pol.grids = NULL) 
length(DesGrid[["map"]]) #959
remove(DesOccs, DesOccs_abbr)

DesMat<-sparse2dense(DesGrid[["comm_dat"]])
str(DesMat)
DesMat<-as.data.frame(DesMat)
rownames(DesMat)[1:5]
colnames(DesMat)[1:5]

DesPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(DesMat)))==T)])
write.tree(DesPhy, "Des/DesPhy.tree")

DesMat<-DesMat[,which(names(DesMat) %in% DesPhy$tip.label)]
DesMat<-DesMat[-which(rowSums(DesMat)==0),]
DesMat[DesMat>1]<-1
DesMat$grids<-row.names(DesMat)
write.csv(DesMat, "Des/DesMat.csv", row.names=FALSE)

DesMap<-DesGrid[["map"]][which(DesGrid[["map"]]$grids %in% rownames(DesMat)),]
length(unique(DesMap$grids)) #947
length(DesMap$grids) #947
crs(DesMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(DesMap, "Des/DesMap.shp", overwrite = TRUE)
remove(DesGrid)




### "Flooded Grasslands & Savannas"

FloGOccs<-rbind(occs[which(occs$BIOME_NAME == "Flooded Grasslands & Savannas"),])

FloGOccs_abbr<-FloGOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(FloGOccs_abbr, "FloG/FloGOccs_abbr.csv", row.names= FALSE)

FloGGrid<-points2comm(dat = FloGOccs_abbr, res = 2, pol.grids = NULL) 
length(FloGGrid[["map"]]) #152
remove(FloGOccs, FloGOccs_abbr)

FloGMat<-sparse2dense(FloGGrid[["comm_dat"]])
str(FloGMat)
FloGMat<-as.data.frame(FloGMat)
rownames(FloGMat)[1:5]
colnames(FloGMat)[1:5]

FloGPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(FloGMat)))==T)])
write.tree(FloGPhy, "FloG/FloGPhy.tree")

FloGMat<-FloGMat[,which(names(FloGMat) %in% FloGPhy$tip.label)]
FloGMat<-FloGMat[-which(rowSums(FloGMat)==0),]
FloGMat[FloGMat>1]<-1
FloGMat$grids<-row.names(FloGMat)
write.csv(FloGMat, "FloG/FloGMat.csv", row.names = FALSE)

FloGMap<-FloGGrid[["map"]][which(FloGGrid[["map"]]$grids %in% rownames(FloGMat)),]
length(unique(FloGMap$grids)) #146
length(FloGMap$grids) #146
crs(FloGMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(FloGMap, "FloG/FloGMap.shp", overwrite = TRUE)
remove(FloGGrid)




### "Mangroves"

MangOccs<-rbind(occs[which(occs$BIOME_NAME == "Mangroves"),])

MangOccs_abbr<-MangOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(MangOccs_abbr, "Mang/MangOccs_abbr.csv", row.names= FALSE)

MangGrid<-points2comm(dat = MangOccs_abbr, res = 2, pol.grids = NULL) 
length(MangGrid[["map"]]) #283
remove(MangOccs, MangOccs_abbr)

MangMat<-sparse2dense(MangGrid[["comm_dat"]])
str(MangMat)
MangMat<-as.data.frame(MangMat)
rownames(MangMat)[1:5]
colnames(MangMat)[1:5]

MangPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(MangMat)))==T)])
write.tree(MangPhy, "Mang/MangPhy.tree")

MangMat<-MangMat[,which(names(MangMat) %in% MangPhy$tip.label)]
MangMat<-MangMat[-which(rowSums(MangMat)==0),]
MangMat[MangMat>1]<-1
MangMat$grids<-row.names(MangMat)
write.csv(MangMat, "Mang/MangMat.csv", row.names = FALSE)

MangMap<-MangGrid[["map"]][which(MangGrid[["map"]]$grids %in% rownames(MangMat)),]
length(unique(MangMap$grids)) #275
length(MangMap$grids) #275
crs(MangMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(MangMap, "Mang/MangMap.shp", overwrite = TRUE)
remove(MangGrid)




### "Tropical & Subtropical Moist Broadleaf Forests"

MBFOccs<-rbind(occs[which(occs$BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests"),])

MBFOccs_abbr<-MBFOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(MBFOccs_abbr, "MBF/MBFOccs_abbr.csv", row.names= FALSE)

MBFGrid<-points2comm(dat = MBFOccs_abbr, res = 2, pol.grids = NULL) 
length(MBFGrid[["map"]]) #1016
remove(MBFOccs, MBFOccs_abbr)

MBFMat<-sparse2dense(MBFGrid[["comm_dat"]])
str(MBFMat)
MBFMat<-as.data.frame(MBFMat)
rownames(MBFMat)[1:5]
colnames(MBFMat)[1:5]

MBFPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(MBFMat)))==T)])
write.tree(MBFPhy, "MBF/MBFPhy.tree")

MBFMat<-MBFMat[,which(names(MBFMat) %in% MBFPhy$tip.label)]
MBFMat<-MBFMat[-which(rowSums(MBFMat)==0),]
MBFMat[MBFMat>1]<-1
MBFMat$grids<-row.names(MBFMat)
write.csv(MBFMat, "MBF/MBFMat.csv", row.names = FALSE)

MBFMap<-MBFGrid[["map"]][which(MBFGrid[["map"]]$grids %in% rownames(MBFMat)),]
length(unique(MBFMap$grids)) #1004
length(MBFMap$grids) #1004
crs(MBFMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(MBFMap, "MBF/MBFMap.shp", overwrite = TRUE)
remove(MBFGrid)






### "Mediterranean Forests, Woodlands & Scrub"

MeditOccs<-rbind(occs[which(occs$BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub"),])

MeditOccs_abbr<-MeditOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(MeditOccs_abbr, "Medit/MeditOccs_abbr.csv", row.names= FALSE)

MeditGrid<-points2comm(dat = MeditOccs_abbr, res = 2, pol.grids = NULL) 
length(MeditGrid[["map"]]) #251
remove(MeditOccs, MeditOccs_abbr)

MeditMat<-sparse2dense(MeditGrid[["comm_dat"]])
str(MeditMat)
MeditMat<-as.data.frame(MeditMat)
rownames(MeditMat)[1:5]
colnames(MeditMat)[1:5]

MeditPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(MeditMat)))==T)])
write.tree(MeditPhy, "Medit/MeditPhy.tree")

MeditMat<-MeditMat[,which(names(MeditMat) %in% MeditPhy$tip.label)]
MeditMat<-MeditMat[-which(rowSums(MeditMat)==0),]
MeditMat[MeditMat>1]<-1
MeditMat$grids<-row.names(MeditMat)
write.csv(MeditMat, "Medit/MeditMat.csv", row.names = FALSE)

MeditMap<-MeditGrid[["map"]][which(MeditGrid[["map"]]$grids %in% rownames(MeditMat)),]
length(unique(MeditMap$grids)) #249
length(MeditMap$grids) #249
crs(MeditMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(MeditMap, "Medit/MeditMap.shp", overwrite = TRUE)
remove(MeditGrid)




### "Montane Grasslands & Shrublands"

MonGOccs<-rbind(occs[which(occs$BIOME_NAME == "Montane Grasslands & Shrublands"),])

MonGOccs_abbr<-MonGOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(MonGOccs_abbr, "MonG/MonGOccs_abbr.csv", row.names= FALSE)

MonGGrid<-points2comm(dat = MonGOccs_abbr, res = 2, pol.grids = NULL) 
length(MonGGrid[["map"]]) #342
remove(MonGOccs, MonGOccs_abbr)

MonGMat<-sparse2dense(MonGGrid[["comm_dat"]])
str(MonGMat)
MonGMat<-as.data.frame(MonGMat)
rownames(MonGMat)[1:5]
colnames(MonGMat)[1:5]

MonGPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(MonGMat)))==T)])
write.tree(MonGPhy, "MonG/MonGPhy.tree")

MonGMat<-MonGMat[,which(names(MonGMat) %in% MonGPhy$tip.label)]
MonGMat<-MonGMat[-which(rowSums(MonGMat)==0),]
MonGMat[MonGMat>1]<-1
MonGMat$grids<-row.names(MonGMat)
write.csv(MonGMat, "MonG/MonGMat.csv", row.names = FALSE)

MonGMap<-MonGGrid[["map"]][which(MonGGrid[["map"]]$grids %in% rownames(MonGMat)),]
length(unique(MonGMap$grids)) #335
length(MonGMap$grids) #335
crs(MonGMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(MonGMap, "MonG/MonGMap.shp", overwrite = TRUE)
remove(MonGGrid)




### "Boreal Forests/Taiga"

TaiOccs<-rbind(occs[which(occs$BIOME_NAME == "Boreal Forests/Taiga"),])

TaiOccs_abbr<-TaiOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TaiOccs_abbr, "Tai/TaiOccs_abbr.csv", row.names= FALSE)

TaiGrid<-points2comm(dat = TaiOccs_abbr, res = 2, pol.grids = NULL) 
length(TaiGrid[["map"]]) #942
remove(TaiOccs, TaiOccs_abbr)

TaiMat<-sparse2dense(TaiGrid[["comm_dat"]])
str(TaiMat)
TaiMat<-as.data.frame(TaiMat)
rownames(TaiMat)[1:5]
colnames(TaiMat)[1:5]

TaiPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TaiMat)))==T)])
write.tree(TaiPhy, "Tai/TaiPhy.tree")

TaiMat<-TaiMat[,which(names(TaiMat) %in% TaiPhy$tip.label)]
TaiMat<-TaiMat[-which(rowSums(TaiMat)==0),]
TaiMat[TaiMat>1]<-1
TaiMat$grids<-row.names(TaiMat)
write.csv(TaiMat, "Tai/TaiMat.csv", row.names = FALSE)

TaiMap<-TaiGrid[["map"]][which(TaiGrid[["map"]]$grids %in% rownames(TaiMat)),]
length(unique(TaiMap$grids)) #941
length(TaiMap$grids) #941
crs(TaiMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TaiMap, "Tai/TaiMap.shp", overwrite = TRUE)
remove(TaiGrid)





### "Temperate Broadleaf & Mixed Forests"

TemBFOccs<-rbind(occs[which(occs$BIOME_NAME == "Temperate Broadleaf & Mixed Forests"),])

TemBFOccs_abbr<-TemBFOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TemBFOccs_abbr, "TemBF/TemBFOccs_abbr.csv", row.names= FALSE)

TemBFGrid<-points2comm(dat = TemBFOccs_abbr, res = 2, pol.grids = NULL) 
length(TemBFGrid[["map"]]) #752
remove(TemBFOccs, TemBFOccs_abbr)

TemBFMat<-sparse2dense(TemBFGrid[["comm_dat"]])
str(TemBFMat)
TemBFMat<-as.data.frame(TemBFMat)
rownames(TemBFMat)[1:5]
colnames(TemBFMat)[1:5]

TemBFPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TemBFMat)))==T)])
write.tree(TemBFPhy, "TemBF/TemBFPhy.tree")

TemBFMat<-TemBFMat[,which(names(TemBFMat) %in% TemBFPhy$tip.label)]
TemBFMat<-TemBFMat[-which(rowSums(TemBFMat)==0),]
TemBFMat[TemBFMat>1]<-1
TemBFMat$grids<-row.names(TemBFMat)
write.csv(TemBFMat, "TemBF/TemBFMat.csv", row.names=FALSE)

TemBFMap<-TemBFGrid[["map"]][which(TemBFGrid[["map"]]$grids %in% rownames(TemBFMat)),]
length(unique(TemBFMap$grids)) #751
length(TemBFMap$grids) #751
crs(TemBFMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TemBFMap, "TemBF/TemBFMap.shp", overwrite = TRUE)
remove(TemBFGrid)





### "Temperate Conifer Forests"

TemCFOccs<-rbind(occs[which(occs$BIOME_NAME == "Temperate Conifer Forests"),])

TemCFOccs_abbr<-TemCFOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TemCFOccs_abbr, "TemCF/TemCFOccs_abbr.csv", row.names= FALSE)

TemCFGrid<-points2comm(dat = TemCFOccs_abbr, res = 2, pol.grids = NULL) 
length(TemCFGrid[["map"]]) #396
remove(TemCFOccs, TemCFOccs_abbr)

TemCFMat<-sparse2dense(TemCFGrid[["comm_dat"]])
str(TemCFMat)
TemCFMat<-as.data.frame(TemCFMat)
rownames(TemCFMat)[1:5]
colnames(TemCFMat)[1:5]

TemCFPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TemCFMat)))==T)])
write.tree(TemCFPhy, "TemCF/TemCFPhy.tree")

TemCFMat<-TemCFMat[,which(names(TemCFMat) %in% TemCFPhy$tip.label)]
TemCFMat<-TemCFMat[-which(rowSums(TemCFMat)==0),]
TemCFMat[TemCFMat>1]<-1
TemCFMat$grids<-row.names(TemCFMat)
write.csv(TemCFMat, "TemCF/TemCFMat.csv", row.names=FALSE)

TemCFMap<-TemCFGrid[["map"]][which(TemCFGrid[["map"]]$grids %in% rownames(TemCFMat)),]
length(unique(TemCFMap$grids)) #391
length(TemCFMap$grids) # 391
crs(TemCFMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TemCFMap, "TemCF/TemCFMap.shp", overwrite = TRUE)
remove(TemCFGrid)





### "Temperate Grasslands, Savannas & Shrublands"

TemGOccs<-rbind(occs[which(occs$BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands"),])

TemGOccs_abbr<-TemGOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TemGOccs_abbr, "TemG/TemGOccs_abbr.csv", row.names= FALSE)

TemGGrid<-points2comm(dat = TemGOccs_abbr, res = 2, pol.grids = NULL) 
length(TemGGrid[["map"]]) #571
remove(TemGOccs, TemGOccs_abbr)

TemGMat<-sparse2dense(TemGGrid[["comm_dat"]])
str(TemGMat)
TemGMat<-as.data.frame(TemGMat)
rownames(TemGMat)[1:5]
colnames(TemGMat)[1:5]

TemGPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TemGMat)))==T)])
write.tree(TemGPhy, "TemG/TemGPhy.tree")

TemGMat<-TemGMat[,which(names(TemGMat) %in% TemGPhy$tip.label)]
TemGMat<-TemGMat[-which(rowSums(TemGMat)==0),]
TemGMat[TemGMat>1]<-1
TemGMat$grids<-row.names(TemGMat)
write.csv(TemGMat, "TemG/TemGMat.csv", row.names = FALSE)

TemGMap<-TemGGrid[["map"]][which(TemGGrid[["map"]]$grids %in% rownames(TemGMat)),]
length(unique(TemGMap$grids)) #569
length(TemGMap$grids) # 569
crs(TemGMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TemGMap, "TemG/TemGMap.shp", overwrite = TRUE)
remove(TemGGrid)




### "Tropical & Subtropical Coniferous Forests"


TroCFOccs<-rbind(occs[which(occs$BIOME_NAME == "Tropical & Subtropical Coniferous Forests"),])

TroCFOccs_abbr<-TroCFOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TroCFOccs_abbr, "TroCF/TroCFOccs_abbr.csv", row.names= FALSE)

TroCFGrid<-points2comm(dat = TroCFOccs_abbr, res = 2, pol.grids = NULL) 
length(TroCFGrid[["map"]]) # 101
remove(TroCFOccs, TroCFOccs_abbr)

TroCFMat<-sparse2dense(TroCFGrid[["comm_dat"]])
str(TroCFMat)
TroCFMat<-as.data.frame(TroCFMat)
rownames(TroCFMat)[1:5]
colnames(TroCFMat)[1:5]

TroCFPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TroCFMat)))==T)])
write.tree(TroCFPhy, "TroCF/TroCFPhy.tree")

TroCFMat<-TroCFMat[,which(names(TroCFMat) %in% TroCFPhy$tip.label)]
TroCFMat<-TroCFMat[-which(rowSums(TroCFMat)==0),]
TroCFMat[TroCFMat>1]<-1
TroCFMat$grids<-row.names(TroCFMat)
write.csv(TroCFMat, "TroCF/TroCFMat.csv", row.names=FALSE)

TroCFMap<-TroCFGrid[["map"]][which(TroCFGrid[["map"]]$grids %in% rownames(TroCFMat)),]
length(unique(TroCFMap$grids)) #100
length(TroCFMap$grids) #100
crs(TroCFMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TroCFMap, "TroCF/TroCFMap.shp", overwrite = TRUE)
remove(TroCFGrid)




### "Tropical & Subtropical Grasslands, Savannas & Shrublands"

TroGOccs<-rbind(occs[which(occs$BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands"),])

TroGOccs_abbr<-TroGOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TroGOccs_abbr, "TroG/TroGOccs_abbr.csv", row.names= FALSE)

TroGGrid<-points2comm(dat = TroGOccs_abbr, res = 2, pol.grids = NULL) 
length(TroGGrid[["map"]]) #733
remove(TroGOccs, TroGOccs_abbr)

TroGMat<-sparse2dense(TroGGrid[["comm_dat"]])
str(TroGMat)
TroGMat<-as.data.frame(TroGMat)
rownames(TroGMat)[1:5]
colnames(TroGMat)[1:5]

TroGPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TroGMat)))==T)])
write.tree(TroGPhy, "TroG/TroGPhy.tree")

TroGMat<-TroGMat[,which(names(TroGMat) %in% TroGPhy$tip.label)]
TroGMat<-TroGMat[-which(rowSums(TroGMat)==0),]
TroGMat[TroGMat>1]<-1
TroGMat$grids<-row.names(TroGMat)
write.csv(TroGMat, "TroG/TroGMat.csv", row.names=FALSE)

TroGMap<-TroGGrid[["map"]][which(TroGGrid[["map"]]$grids %in% rownames(TroGMat)),]
length(unique(TroGMap$grids)) #726
length(TroGMap$grids) #726
crs(TroGMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TroGMap, "TroG/TroGMap.shp", overwrite = TRUE)
remove(TroGGrid)





### "Tundra"

TunOccs<-rbind(occs[which(occs$BIOME_NAME == "Tundra"),])

TunOccs_abbr<-TunOccs[,c("species", "decimallongitude", "decimallatitude")]
write.csv(TunOccs_abbr, "Tun/TunOccs_abbr.csv", row.names= FALSE)

TunGrid<-points2comm(dat = TunOccs_abbr, res = 2, pol.grids = NULL) 
length(TunGrid[["map"]]) #1024
remove(TunOccs, TunOccs_abbr)

TunMat<-sparse2dense(TunGrid[["comm_dat"]])
str(TunMat)
TunMat<-as.data.frame(TunMat)
rownames(TunMat)[1:5]
colnames(TunMat)[1:5]

TunPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, names(TunMat)))==T)])
write.tree(TunPhy, "Tun/TunPhy.tree")

TunMat<-TunMat[,which(names(TunMat) %in% TunPhy$tip.label)]
TunMat<-TunMat[-which(rowSums(TunMat)==0),]
TunMat[TunMat>1]<-1
TunMat$grids<-row.names(TunMat)
write.csv(TunMat, "Tun/TunMat.csv", row.names=FALSE)

TunMap<-TunGrid[["map"]][which(TunGrid[["map"]]$grids %in% rownames(TunMat)),]
length(unique(TunMap$grids)) #1007
length(TunMap$grids) #1007
crs(TunMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(TunMap, "Tun/TunMap.shp", overwrite = TRUE)
remove(TunGrid)





### Global

AllOccs_abbr<-occs[,c("species", "decimallongitude", "decimallatitude")]

AllGrid<-points2comm(dat = AllOccs_abbr, res = 2, pol.grids = NULL)
length(AllGrid[["map"]])

AllMat<-sparse2dense(AllGrid[["comm_dat"]])
str(AllMat)
AllMat<-as.data.frame(AllMat)
rownames(AllMat)[1:5]
colnames(AllMat)[1:5]

AllPhy<-drop.tip(SnB_updated, SnB_updated$tip.label[
  which(is.na(match(SnB_updated$tip.label, colnames(AllMat)))==T)])
write.tree(AllPhy, "All/AllPhy.tree")

AllMat<-AllMat[,which(colnames(AllMat) %in% AllPhy$tip.label)]
AllMat<-AllMat[-which(rowSums(AllMat)==0),]
AllMat[AllMat>1]<-1
AllMat$grids<-row.names(AllDBFMat)
write.csv(AllMat, "All/AllMat.csv", row.names= FALSE)


AllMap<-AllGrid[["map"]][which(AllGrid[["map"]]$grids %in% rownames(AllMat)),]
length(unique(AllMap$grids))
length(AllMap$grids)
crs(AllMap)<-"epsg:4326" # same as biome shapefile, and as rasters
writeVector(AllMap, "All/AllMap.shp", overwrite = TRUE)


## Spatial grid formatting (clip, centroid coordinates, etc.) ----

###  "Tropical & Subtropical Dry Broadleaf Forests"

# Set up grid-level map for clipping to biome boundaries
DBFMapsf<-DBFMap
DBFMapsf
DBFMapsf<-st_make_valid(DBFMapsf)

# Clip gridded biome map to original Dinerstein et al. (2017) biome map
DBFMap_clipped<-st_intersection(DBFMapsf, shp[shp$BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests",])
DBFMap_clipped<-DBFMap_clipped[,-c(4:18)]

# Remove weird geometry
DBFMap_clipped<-DBFMap_clipped[-which(st_geometry_type(DBFMap_clipped) == "GEOMETRYCOLLECTION"),]
# Cast as polygon
DBFMap_clipped<-st_cast(DBFMap_clipped, "POLYGON") # this keeps giving an error but
# t is not messing with any of the polygons so, shrug
st_write(DBFMap_clipped, "DBF/DBFMap_clipped.shp", append = FALSE)

# After clipping, some grid cells become multipolygons. We want one polygon for GDM input, so
# Assign the grid cell to the polygon with the largest area and get the area for each poly
DBFMap_1poly <- DBFMap_clipped %>% 
  mutate(area = st_area(geometry))
# Filter sf
DBFMap_1poly <- DBFMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()

# Double check that they're all polys
DBFMap_1poly<-st_cast(DBFMap_1poly, "POLYGON")
unique(st_geometry_type(DBFMap_1poly))

st_write(DBFMap_1poly, "DBF/DBFMap_1poly.shp", append = FALSE)

# Get Lat/Long of grid cell centroids for GDM input
DBFMap_1poly<-vect(DBFMap_1poly)
DBFCentz<-centroids(DBFMap_1poly, inside = TRUE)
writeVector(DBFCentz, "DBF/DBF_centroids.shp", overwrite=TRUE)

#convert to df for merging with comm mat
DBFCoordz<-as.data.frame(crds(DBFCentz))
colnames(DBFCoordz)<-c("lon", "lat")
DBFCoordz$grids<-DBFMap_1poly$grids
remove(DBFMapsf)

# Join coordz and comm mat
DBFMat<-as.data.frame(fread("DBF/DBFMat.csv"))
DBFMatCoordz<-left_join(DBFCoordz, DBFMat, by = "grids")
DBFMatCoordz[1:5, 1:15]
write.csv(DBFMatCoordz, "DBF/DBFMatCoords.csv", row.names = FALSE)




#REPEAT THIS PROCESS 14 TIMES (13 MORE BIOMES + 1 GLOBAL)



### "Deserts & Xeric Shrublands"

DesMapsf<-read_sf("Des/DesMap.shp")
DesMapsf
DesMapsf<-st_make_valid(DesMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
DesMap_clipped<-st_intersection(DesMapsf, shp[shp$BIOME_NAME == "Deserts & Xeric Shrublands",])
DesMap_clipped<-DesMap_clipped[,-c(4:18)]
DesMap_clipped<-st_cast(DesMap_clipped, "POLYGON")
st_write(DesMap_clipped, "Des/DesMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
DesMap_1poly <- DesMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
DesMap_1poly <- DesMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(DesMap_1poly))
DesMap_1poly<-st_cast(DesMap_1poly, "POLYGON")
st_write(DesMap_1poly, "Des/DesMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
DesMap_1poly<-vect(DesMap_1poly)
DesCentz<-centroids(DesMap_1poly, inside = TRUE)
writeVector(DesCentz, "Des/Des_centroids.shp", overwrite=TRUE)

DesCoordz<-as.data.frame(crds(DesCentz))
colnames(DesCoordz)<-c("lon", "lat")
DesCoordz$grids<-DesMap_1poly$grids
remove(DesMapsf)

# Join coordz + submat
DesMat<-as.data.frame(fread("Des/DesMat.csv"))
DesMatCoordz<-left_join(DesCoordz, DesMat, by = "grids")
DesMatCoordz[1:5, 1:15]
write.csv(DesMatCoordz, "Des/DesMatCoords.csv", row.names = FALSE)





### "Flooded Grasslands & Savannas"

FloGMapsf<-read_sf("FloG/FloGMap.shp")
FloGMapsf
FloGMapsf<-st_make_valid(FloGMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
FloGMap_clipped<-st_intersection(FloGMapsf, shp[shp$BIOME_NAME == "Flooded Grasslands & Savannas",])
FloGMap_clipped<-FloGMap_clipped[,-c(4:18)]
FloGMap_clipped<-st_cast(FloGMap_clipped, "POLYGON")
st_write(FloGMap_clipped, "FloG/FloGMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
FloGMap_1poly <- FloGMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
FloGMap_1poly <- FloGMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(FloGMap_1poly))
FloGMap_1poly<-st_cast(FloGMap_1poly, "POLYGON")
st_write(FloGMap_1poly, "FloG/FloGMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
FloGMap_1poly<-vect(FloGMap_1poly)
FloGCentz<-centroids(FloGMap_1poly, inside = TRUE)
writeVector(FloGCentz, "FloG/FloG_centroids.shp", overwrite=TRUE)

FloGCoordz<-as.data.frame(crds(FloGCentz))
colnames(FloGCoordz)<-c("lon", "lat")
FloGCoordz$grids<-FloGMap_1poly$grids
remove(FloGMapsf)

# Join coordz + submat
FloGMat<-as.data.frame(fread("FloG/FloGMat.csv"))
FloGMatCoordz<-left_join(FloGCoordz, FloGMat, by = "grids")
FloGMatCoordz[1:5, 1:15]
write.csv(FloGMatCoordz, "FloG/FloGMatCoords.csv", row.names = FALSE)






### "Mangroves

MangMapsf<-read_sf("Mang/MangMap.shp")
MangMapsf
MangMapsf<-st_make_valid(MangMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
MangMap_clipped<-st_intersection(MangMapsf, shp[shp$BIOME_NAME == "Mangroves",])
MangMap_clipped<-MangMap_clipped[,-c(4:18)]
MangMap_clipped<-st_cast(MangMap_clipped, "POLYGON")
st_write(MangMap_clipped, "Mang/MangMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
MangMap_1poly <- MangMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
MangMap_1poly <- MangMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(MangMap_1poly))
MangMap_1poly<-st_cast(MangMap_1poly, "POLYGON")
st_write(MangMap_1poly, "Mang/MangMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
MangMap_1poly<-vect(MangMap_1poly)
MangCentz<-centroids(MangMap_1poly, inside = TRUE)
writeVector(MangCentz, "Mang/Mang_centroids.shp", overwrite=TRUE)

MangCoordz<-as.data.frame(crds(MangCentz))
colnames(MangCoordz)<-c("lon", "lat")
MangCoordz$grids<-MangMap_1poly$grids
remove(MangMapsf)

# Join coordz + submat
MangMat<-as.data.frame(fread("Mang/MangMat.csv"))
MangMatCoordz<-left_join(MangCoordz, MangMat, by = "grids")
MangMatCoordz[1:5, 1:15]
write.csv(MangMatCoordz, "Mang/MangMatCoords.csv", row.names = FALSE)






### "Tropical & Subtropical Moist Broadleaf Forests"

MBFMapsf<-read_sf("MBF/MBFMap.shp")
MBFMapsf
MBFMapsf<-st_make_valid(MBFMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
MBFMap_clipped<-st_intersection(MBFMapsf, shp[shp$BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests",])
MBFMap_clipped<-MBFMap_clipped[,-c(4:18)]
MBFMap_clipped<-st_make_valid(MBFMap_clipped)
unique(st_geometry_type(MBFMap_clipped))
MBFMap_clipped<-MBFMap_clipped[-which(st_geometry_type(MBFMap_clipped) == "GEOMETRYCOLLECTION"),]
MBFMap_clipped<-st_cast(MBFMap_clipped, "POLYGON")
st_write(MBFMap_clipped, "MBF/MBFMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
MBFMap_1poly <- MBFMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
MBFMap_1poly <- MBFMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
MBFMap_1poly<-st_cast(MBFMap_1poly, "POLYGON")
st_write(MBFMap_1poly, "MBF/MBFMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
MBFMap_1poly<-vect(MBFMap_1poly)
MBFCentz<-centroids(MBFMap_1poly, inside = TRUE)
writeVector(MBFCentz, "MBF/MBF_centroids.shp", overwrite=TRUE)

MBFCoordz<-as.data.frame(crds(MBFCentz))
colnames(MBFCoordz)<-c("lon", "lat")
MBFCoordz$grids<-MBFMap_1poly$grids
remove(MBFMapsf)

# Join coordz + submat
MBFMat<-as.data.frame(fread("MBF/MBFMat.csv"))
MBFMatCoordz<-left_join(MBFCoordz, MBFMat, by = "grids")
MBFMatCoordz[1:5, 1:15]
write.csv(MBFMatCoordz, "MBF/MBFMatCoords.csv", row.names = FALSE)





### "Mediterranean Forests, Woodlands & Scrub"

MeditMapsf<-read_sf("Medit/MeditMap.shp")
MeditMapsf
MeditMapsf<-st_make_valid(MeditMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
MeditMap_clipped<-st_intersection(MeditMapsf, shp[shp$BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub",])
MeditMap_clipped<-MeditMap_clipped[,-c(4:18)]
MeditMap_clipped<-st_cast(MeditMap_clipped, "POLYGON")
st_write(MeditMap_clipped, "Medit/MeditMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
MeditMap_1poly <- MeditMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
MeditMap_1poly <- MeditMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(MeditMap_1poly))
MeditMap_1poly<-st_cast(MeditMap_1poly, "POLYGON")
st_write(MeditMap_1poly, "Medit/MeditMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
MeditMap_1poly<-vect(MeditMap_1poly)
MeditCentz<-centroids(MeditMap_1poly, inside = TRUE)
writeVector(MeditCentz, "Medit/Medit_centroids.shp", overwrite=TRUE)

MeditCoordz<-as.data.frame(crds(MeditCentz))
colnames(MeditCoordz)<-c("lon", "lat")
MeditCoordz$grids<-MeditMap_1poly$grids
remove(MeditMapsf)

# Join coordz + submat
MeditMat<-as.data.frame(fread("Medit/MeditMat.csv"))
MeditMatCoordz<-left_join(MeditCoordz, MeditMat, by = "grids")
MeditMatCoordz[1:5, 1:15]
write.csv(MeditMatCoordz, "Medit/MeditMatCoords.csv", row.names = FALSE)





### "Montane Grasslands & Shrublands"

MonGMapsf<-read_sf("MonG/MonGMap.shp")
MonGMapsf
MonGMapsf<-st_make_valid(MonGMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
MonGMap_clipped<-st_intersection(MonGMapsf, shp[shp$BIOME_NAME == "Montane Grasslands & Shrublands",])
MonGMap_clipped<-MonGMap_clipped[,-c(4:18)]
MonGMap_clipped<-st_cast(MonGMap_clipped, "POLYGON")
st_write(MonGMap_clipped, "MonG/MonGMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
MonGMap_1poly <- MonGMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
MonGMap_1poly <- MonGMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(MonGMap_1poly))
MonGMap_1poly<-st_cast(MonGMap_1poly, "POLYGON")
st_write(MonGMap_1poly, "MonG/MonGMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
MonGMap_1poly<-vect(MonGMap_1poly)
MonGCentz<-centroids(MonGMap_1poly, inside = TRUE)
writeVector(MonGCentz, "MonG/MonG_centroids.shp", overwrite=TRUE)

MonGCoordz<-as.data.frame(crds(MonGCentz))
colnames(MonGCoordz)<-c("lon", "lat")
MonGCoordz$grids<-MonGMap_1poly$grids
remove(MonGMapsf)

# Join coordz + submat
MonGMat<-as.data.frame(fread("MonG/MonGMat.csv"))
MonGMatCoordz<-left_join(MonGCoordz, MonGMat, by = "grids")
MonGMatCoordz[1:5, 1:15]
write.csv(MonGMatCoordz, "MonG/MonGMatCoords.csv", row.names = FALSE)





### "Boreal Forests/Taiga"

TaiMapsf<-read_sf("Tai/TaiMap.shp")
TaiMapsf
TaiMapsf<-st_make_valid(TaiMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TaiMap_clipped<-st_intersection(TaiMapsf, shp[shp$BIOME_NAME == "Boreal Forests/Taiga",])
TaiMap_clipped<-TaiMap_clipped[,-c(4:18)]
TaiMap_clipped<-st_cast(TaiMap_clipped, "POLYGON")
st_write(TaiMap_clipped, "Tai/TaiMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TaiMap_1poly <- TaiMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TaiMap_1poly <- TaiMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TaiMap_1poly))
TaiMap_1poly<-st_cast(TaiMap_1poly, "POLYGON")
st_write(TaiMap_1poly, "Tai/TaiMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TaiMap_1poly<-vect(TaiMap_1poly)
TaiCentz<-centroids(TaiMap_1poly, inside = TRUE)
writeVector(TaiCentz, "Tai/Tai_centroids.shp", overwrite=TRUE)

TaiCoordz<-as.data.frame(crds(TaiCentz))
colnames(TaiCoordz)<-c("lon", "lat")
TaiCoordz$grids<-TaiMap_1poly$grids
remove(TaiMapsf)

# Join coordz + submat
TaiMat<-as.data.frame(fread("Tai/TaiMat.csv"))
TaiMatCoordz<-left_join(TaiCoordz, TaiMat, by = "grids")
TaiMatCoordz[1:5, 1:15]
write.csv(TaiMatCoordz, "Tai/TaiMatCoords.csv", row.names = FALSE)





### "Temperate Broadleaf & Mixed Forests"

TemBFMapsf<-read_sf("TemBF/TemBFMap.shp")
TemBFMapsf
TemBFMapsf<-st_make_valid(TemBFMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TemBFMap_clipped<-st_intersection(TemBFMapsf, shp[shp$BIOME_NAME == "Temperate Broadleaf & Mixed Forests",])
TemBFMap_clipped<-TemBFMap_clipped[,-c(4:18)]
TemBFMap_clipped<-st_cast(TemBFMap_clipped, "POLYGON")
st_write(TemBFMap_clipped, "TemBF/TemBFMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TemBFMap_1poly <- TemBFMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TemBFMap_1poly <- TemBFMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TemBFMap_1poly))
TemBFMap_1poly<-st_cast(TemBFMap_1poly, "POLYGON")
st_write(TemBFMap_1poly, "TemBF/TemBFMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TemBFMap_1poly<-vect(TemBFMap_1poly)
TemBFCentz<-centroids(TemBFMap_1poly, inside = TRUE)
writeVector(TemBFCentz, "TemBF/TemBF_centroids.shp", overwrite=TRUE)

TemBFCoordz<-as.data.frame(crds(TemBFCentz))
colnames(TemBFCoordz)<-c("lon", "lat")
TemBFCoordz$grids<-TemBFMap_1poly$grids
remove(TemBFMapsf)

# Join coordz + submat
TemBFMat<-as.data.frame(fread("TemBF/TemBFMat.csv"))
TemBFMatCoordz<-left_join(TemBFCoordz, TemBFMat, by = "grids")
TemBFMatCoordz[1:5, 1:15]
write.csv(TemBFMatCoordz, "TemBF/TemBFMatCoords.csv", row.names = FALSE)





### "Temperate Conifer Forests"

TemCFMapsf<-read_sf("TemCF/TemCFMap.shp")
TemCFMapsf
TemCFMapsf<-st_make_valid(TemCFMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TemCFMap_clipped<-st_intersection(TemCFMapsf, shp[shp$BIOME_NAME == "Temperate Conifer Forests",])
TemCFMap_clipped<-TemCFMap_clipped[,-c(4:18)]
TemCFMap_clipped<-st_cast(TemCFMap_clipped, "POLYGON")
st_write(TemCFMap_clipped, "TemCF/TemCFMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TemCFMap_1poly <- TemCFMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TemCFMap_1poly <- TemCFMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TemCFMap_1poly))
TemCFMap_1poly<-st_cast(TemCFMap_1poly, "POLYGON")
st_write(TemCFMap_1poly, "TemCF/TemCFMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TemCFMap_1poly<-vect(TemCFMap_1poly)
TemCFCentz<-centroids(TemCFMap_1poly, inside = TRUE)
writeVector(TemCFCentz, "TemCF/TemCF_centroids.shp", overwrite=TRUE)

TemCFCoordz<-as.data.frame(crds(TemCFCentz))
colnames(TemCFCoordz)<-c("lon", "lat")
TemCFCoordz$grids<-TemCFMap_1poly$grids
remove(TemCFMapsf)

# Join coordz + submat
TemCFMat<-as.data.frame(fread("TemCF/TemCFMat.csv"))
TemCFMatCoordz<-left_join(TemCFCoordz, TemCFMat, by = "grids")
TemCFMatCoordz[1:5, 1:15]
write.csv(TemCFMatCoordz, "TemCF/TemCFMatCoords.csv", row.names = FALSE)






### "Temperate Grasslands, Savannas & Shrublands"

TemGMapsf<-read_sf("TemG/TemGMap.shp")
TemGMapsf
TemGMapsf<-st_make_valid(TemGMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TemGMap_clipped<-st_intersection(TemGMapsf, shp[shp$BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands",])
TemGMap_clipped<-TemGMap_clipped[,-c(4:18)]
TemGMap_clipped<-st_cast(TemGMap_clipped, "POLYGON")
st_write(TemGMap_clipped, "TemG/TemGMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TemGMap_1poly <- TemGMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TemGMap_1poly <- TemGMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TemGMap_1poly))
TemGMap_1poly<-st_cast(TemGMap_1poly, "POLYGON")
st_write(TemGMap_1poly, "TemG/TemGMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TemGMap_1poly<-vect(TemGMap_1poly)
TemGCentz<-centroids(TemGMap_1poly, inside = TRUE)
writeVector(TemGCentz, "TemG/TemG_centroids.shp", overwrite=TRUE)

TemGCoordz<-as.data.frame(crds(TemGCentz))
colnames(TemGCoordz)<-c("lon", "lat")
TemGCoordz$grids<-TemGMap_1poly$grids
remove(TemGMapsf)

# Join coordz + submat
TemGMat<-as.data.frame(fread("TemG/TemGMat.csv"))
TemGMatCoordz<-left_join(TemGCoordz, TemGMat, by = "grids")
TemGMatCoordz[1:5, 1:15]
write.csv(TemGMatCoordz, "TemG/TemGMatCoords.csv", row.names = FALSE)





### "Tropical & Subtropical Coniferous Forests"

TroCFMapsf<-read_sf("TroCF/TroCFMap.shp")
TroCFMapsf
TroCFMapsf<-st_make_valid(TroCFMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TroCFMap_clipped<-st_intersection(TroCFMapsf, shp[shp$BIOME_NAME == "Tropical & Subtropical Coniferous Forests",])
TroCFMap_clipped<-TroCFMap_clipped[,-c(4:18)]
TroCFMap_clipped<-st_cast(TroCFMap_clipped, "POLYGON")
st_write(TroCFMap_clipped, "TroCF/TroCFMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TroCFMap_1poly <- TroCFMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TroCFMap_1poly <- TroCFMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TroCFMap_1poly))
TroCFMap_1poly<-st_cast(TroCFMap_1poly, "POLYGON")
st_write(TroCFMap_1poly, "TroCF/TroCFMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TroCFMap_1poly<-vect(TroCFMap_1poly)
TroCFCentz<-centroids(TroCFMap_1poly, inside = TRUE)
writeVector(TroCFCentz, "TroCF/TroCF_centroids.shp", overwrite=TRUE)

TroCFCoordz<-as.data.frame(crds(TroCFCentz))
colnames(TroCFCoordz)<-c("lon", "lat")
TroCFCoordz$grids<-TroCFMap_1poly$grids
remove(TroCFMapsf)

# Join coordz + submat
TroCFMat<-as.data.frame(fread("TroCF/TroCFMat.csv"))
TroCFMatCoordz<-left_join(TroCFCoordz, TroCFMat, by = "grids")
TroCFMatCoordz[1:5, 1:15]
write.csv(TroCFMatCoordz, "TroCF/TroCFMatCoords.csv", row.names = FALSE)





### "Tropical & Subtropical Grasslands, Savannas & Shrublands"

TroGMapsf<-read_sf("TroG/TroGMap.shp")
TroGMapsf
TroGMapsf<-st_make_valid(TroGMapsf)
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TroGMap_clipped<-st_intersection(TroGMapsf, shp[shp$BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands",])
TroGMap_clipped<-TroGMap_clipped[,-c(4:18)]
TroGMap_clipped<-st_cast(TroGMap_clipped, "POLYGON")
st_write(TroGMap_clipped, "TroG/TroGMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TroGMap_1poly <- TroGMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TroGMap_1poly <- TroGMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TroGMap_1poly))
TroGMap_1poly<-st_cast(TroGMap_1poly, "POLYGON")
st_write(TroGMap_1poly, "TroG/TroGMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TroGMap_1poly<-vect(TroGMap_1poly)
TroGCentz<-centroids(TroGMap_1poly, inside = TRUE)
writeVector(TroGCentz, "TroG/TroG_centroids.shp", overwrite=TRUE)

TroGCoordz<-as.data.frame(crds(TroGCentz))
colnames(TroGCoordz)<-c("lon", "lat")
TroGCoordz$grids<-TroGMap_1poly$grids
remove(TroGMapsf)

# Join coordz + submat
TroGMat<-as.data.frame(fread("TroG/TroGMat.csv"))
TroGMatCoordz<-left_join(TroGCoordz, TroGMat, by = "grids")
TroGMatCoordz[1:5, 1:15]
write.csv(TroGMatCoordz, "TroG/TroGMatCoords.csv", row.names = FALSE)




### "Tundra"

TunMapsf<-read_sf("Tun/TunMap.shp")
TunMapsf
TunMapsf<-st_make_valid(TunMapsf)
TunMapsf<-TunMapsf[-which(st_is_valid(TunMapsf) == FALSE),]
#clip gridded biome map to original Dinerstein et al. (2017) biome map
TunMap_clipped<-st_intersection(TunMapsf, shp[shp$BIOME_NAME == "Tundra",])
TunMap_clipped<-TunMap_clipped[,-c(4:18)]
TunMap_clipped<-st_cast(TunMap_clipped, "POLYGON")
st_write(TunMap_clipped, "Tun/TunMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
TunMap_1poly <- TunMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
TunMap_1poly <- TunMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(TunMap_1poly))
TunMap_1poly<-st_cast(TunMap_1poly, "POLYGON")
st_write(TunMap_1poly, "Tun/TunMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
TunMap_1poly<-vect(TunMap_1poly)
TunCentz<-centroids(TunMap_1poly, inside = TRUE)
writeVector(TunCentz, "Tun/Tun_centroids.shp", overwrite=TRUE)

TunCoordz<-as.data.frame(crds(TunCentz))
colnames(TunCoordz)<-c("lon", "lat")
TunCoordz$grids<-TunMap_1poly$grids
remove(TunMapsf)

# Join coordz + submat
TunMat<-as.data.frame(fread("Tun/TunMat.csv"))
TunMatCoordz<-left_join(TunCoordz, TunMat, by = "grids")
TunMatCoordz[1:5, 1:15]
write.csv(TunMatCoordz, "Tun/TunMatCoords.csv", row.names = FALSE)




### Global

AllMapsf<-read_sf("All/AllMap.shp")
AllMapsf
AllMapsf<-st_make_valid(AllMapsf)
AllMapsf<-AllMapsf[-which(st_is_valid(AllMapsf)== FALSE),]
#clip gridded biome map to original Dinerstein et al. (2017) biome map
AllMap_clipped<-st_intersection(AllMapsf, shp)
AllMap_clipped<-AllMap_clipped[,-c(4:18)]
AllMap_clipped<-AllMap_clipped[-which(st_geometry_type(AllMap_clipped) == "GEOMETRYCOLLECTION"),]
AllMap_clipped<-AllMap_clipped[-which(st_geometry_type(AllMap_clipped) == "MULTILINESTRING"),]
AllMap_clipped<-st_cast(AllMap_clipped, "POLYGON")
AllMap_clipped<-st_make_valid(AllMap_clipped)
st_write(AllMap_clipped, "All/AllMap_clipped.shp", append = FALSE)
#after clipping, we are left with some grid cells being multipolygons
#assign the grid cell to the polygon with the largest area
#get area for each poly
AllMap_1poly <- AllMap_clipped %>% 
  mutate(area = st_area(geometry))
#filter sf
# none of these filter things actually worked i think....
AllMap_1poly <- AllMap_1poly %>%
  group_by(grids) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()
unique(st_geometry_type(AllMap_1poly))
AllMap_1poly<-st_cast(AllMap_1poly, "POLYGON")
st_write(AllMap_1poly, "All/AllMap_1poly.shp", append = FALSE)

# Get Lat/Long for GDM input
AllMap_1poly<-vect(AllMap_1poly)
AllCentz<-centroids(AllMap_1poly, inside = TRUE)
writeVector(AllCentz, "All/All_centroids.shp", overwrite=TRUE)

AllCoordz<-as.data.frame(crds(AllCentz))
colnames(AllCoordz)<-c("lon", "lat")
AllCoordz$grids<-AllMap_1poly$grids
remove(AllMapsf)

# Join coordz + submat
AllMat<-as.data.frame(fread("All/AllMat.csv"))
AllMatCoordz<-left_join(AllCoordz, AllMat, by = "grids")
AllMatCoordz[1:5, 1:15]
write.csv(AllMatCoordz, "All/AllMatCoords.csv", row.names = FALSE)





## Connectivity metrics ----

# Again, this is done by biome. See above or MS for biome abbreviations


# Start by making the clusters
# Get a vector with the grid cells that touch each other 
thistouches<-st_touches(DBFMap_clipped)

# Create an adjacency list
g<-graph.adjlist(thistouches)
# Get the number of clusters for each biome
c<-components(g)
# Assign a cluster column to the sf and write it
DBFMap_clipped$cluster<-c$membership
length(unique(DBFMap_clipped$cluster))
st_write(DBFMap_clipped, "DBF/DBFMap_clust.shp", append=FALSE)
# Make the clusters factors
DBFMap_clipped$cluster<-as.factor(DBFMap_clipped$cluster)

# Plot to double check
ggplot() +
  geom_sf(data = DBFMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

# Looks good! let's get some connectivitry measurements using the clusters
# Make a new df with Number of Clusters (NC) and area of each cluster
DBFCon <- DBFMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

# Minimum distance -
min_dist<-st_distance(DBFCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

# Add min_dist back into our connectivity df
DBFCon$min_dist<-NA
for(i in 1:nrow(DBFCon)){
  print(i)
  DBFCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

#CPL - Mean of shortest paths between the clusters
DBFCon$CPL<-mean(DBFCon$min_dist)

# SLC – Area of the largest cluster
DBFCon$SLC<-max(DBFCon$total_area)

# MSC – Mean area of clusters 
DBFCon$MSC<-mean(DBFCon$total_area)

# AWF – Area-weighted Flux. Evaluates the flow, weighted by area, between clusters
k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * DBFCon$total_area[i] *
        DBFCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
# Create a long matrix for summarizing
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(DBFCon), each=nrow(DBFCon)))
# Get the total AWF for each cluster
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
# Add each cluster's AWF back into connectivity df
DBFCon$AWF<-flux_tot$total_flux
# Get a mean, biome-scale AWF
DBFCon$meanAWF<-mean(DBFCon$AWF)

# Plot by cluster for funsies
ggplot() +
  geom_sf(data = DBFCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

DBFCon<-as.data.frame(DBFCon)
DBFCon<-as.data.frame(DBFCon[,-which(colnames(DBFCon) == "geometry")])

write.csv(DBFCon, "DBF/DBFCon.csv", row.names = FALSE )



# REPEAT THIS 14 TIMES 





DesMap_clipped<-read_sf("Des/DesMap_clipped.shp")
thistouches<-st_touches(DesMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
DesMap_clipped$cluster<-c$membership
length(unique(DesMap_clipped$cluster))
st_write(DesMap_clipped, "Des/DesMap_clust.shp", append=FALSE)
DesMap_clipped$cluster<-as.factor(DesMap_clipped$cluster)

ggplot() +
  geom_sf(data = DesMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

DesCon <- DesMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(DesCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

DesCon$min_dist<-NA
for(i in 1:nrow(DesCon)){
  print(i)
  DesCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

DesCon$CPL<-mean(DesCon$min_dist)

DesCon$SLC<-max(DesCon$total_area)

DesCon$MSC<-mean(DesCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * DesCon$total_area[i] *
        DesCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(DesCon), each=nrow(DesCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
DesCon$AWF<-flux_tot$total_flux
DesCon$meanAWF<-mean(DesCon$AWF)


ggplot() +
  geom_sf(data = DesCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

DesCon<-as.data.frame(DesCon)
DesCon<-as.data.frame(DesCon[,-which(colnames(DesCon) == "geometry")])

write.csv(DesCon, "Des/DesCon.csv", row.names = FALSE )






FloGMap_clipped<-read_sf("FloG/FloGMap_clipped.shp")
thistouches<-st_touches(FloGMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
FloGMap_clipped$cluster<-c$membership
length(unique(FloGMap_clipped$cluster))
st_write(FloGMap_clipped, "FloG/FloGMap_clust.shp", append=FALSE)
FloGMap_clipped$cluster<-as.factor(FloGMap_clipped$cluster)

ggplot() +
  geom_sf(data = FloGMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

FloGCon <- FloGMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(FloGCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

FloGCon$min_dist<-NA
for(i in 1:nrow(FloGCon)){
  print(i)
  FloGCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

FloGCon$CPL<-mean(FloGCon$min_dist)

FloGCon$SLC<-max(FloGCon$total_area)

FloGCon$MSC<-mean(FloGCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * FloGCon$total_area[i] *
        FloGCon$total_area[j]
    }
  }
}
colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(FloGCon), each=nrow(FloGCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
FloGCon$AWF<-flux_tot$total_flux
FloGCon$meanAWF<-mean(FloGCon$AWF)


ggplot() +
  geom_sf(data = FloGCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

FloGCon<-as.data.frame(FloGCon)
FloGCon<-as.data.frame(FloGCon[,-which(colnames(FloGCon) == "geometry")])

write.csv(FloGCon, "FloG/FloGCon.csv", row.names = FALSE )






MangMap_clipped<-read_sf("Mang/MangMap_clipped.shp")
thistouches<-st_touches(MangMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
MangMap_clipped$cluster<-c$membership
length(unique(MangMap_clipped$cluster))
st_write(MangMap_clipped, "Mang/MangMap_clust.shp", append=FALSE)
MangMap_clipped$cluster<-as.factor(MangMap_clipped$cluster)

ggplot() +
  geom_sf(data = MangMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

MangCon <- MangMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(MangCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

MangCon$min_dist<-NA
for(i in 1:nrow(MangCon)){
  print(i)
  MangCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

MangCon$CPL<-mean(MangCon$min_dist)

MangCon$SLC<-max(MangCon$total_area)

MangCon$MSC<-mean(MangCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * MangCon$total_area[i] *
        MangCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(MangCon), each=nrow(MangCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
MangCon$AWF<-flux_tot$total_flux
MangCon$meanAWF<-mean(MangCon$AWF)


ggplot() +
  geom_sf(data = MangCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

MangCon<-as.data.frame(MangCon)
MangCon<-as.data.frame(MangCon[,-which(colnames(MangCon) == "geometry")])

write.csv(MangCon, "Mang/MangCon.csv", row.names = FALSE )







MBFMap_clipped<-read_sf("MBF/MBFMap_clipped.shp")
thistouches<-st_touches(MBFMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
MBFMap_clipped$cluster<-c$membership
length(unique(MBFMap_clipped$cluster))
st_write(MBFMap_clipped, "MBF/MBFMap_clust.shp", append=FALSE)
MBFMap_clipped$cluster<-as.factor(MBFMap_clipped$cluster)

ggplot() +
  geom_sf(data = MBFMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

MBFCon <- MBFMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(MBFCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

MBFCon$min_dist<-NA
for(i in 1:nrow(MBFCon)){
  print(i)
  MBFCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

MBFCon$CPL<-mean(MBFCon$min_dist)

MBFCon$SLC<-max(MBFCon$total_area)

MBFCon$MSC<-mean(MBFCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 232:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * MBFCon$total_area[i] *
        MBFCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(MBFCon), each=nrow(MBFCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
MBFCon$AWF<-flux_tot$total_flux
MBFCon$meanAWF<-mean(MBFCon$AWF)


ggplot() +
  geom_sf(data = MBFCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

MBFCon<-as.data.frame(MBFCon)
MBFCon<-as.data.frame(MBFCon[,-which(colnames(MBFCon) == "geometry")])

write.csv(MBFCon, "MBF/MBFCon.csv", row.names = FALSE )







MeditMap_clipped<-read_sf("Medit/MeditMap_clipped.shp")
thistouches<-st_touches(MeditMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
MeditMap_clipped$cluster<-c$membership
length(unique(MeditMap_clipped$cluster))
st_write(MeditMap_clipped, "Medit/MeditMap_clust.shp", append=FALSE)
MeditMap_clipped$cluster<-as.factor(MeditMap_clipped$cluster)

ggplot() +
  geom_sf(data = MeditMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

MeditCon <- MeditMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(MeditCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

MeditCon$min_dist<-NA
for(i in 1:nrow(MeditCon)){
  print(i)
  MeditCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

MeditCon$CPL<-mean(MeditCon$min_dist)

MeditCon$SLC<-max(MeditCon$total_area)

MeditCon$MSC<-mean(MeditCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * MeditCon$total_area[i] *
        MeditCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(MeditCon), each=nrow(MeditCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
MeditCon$AWF<-flux_tot$total_flux
MeditCon$meanAWF<-mean(MeditCon$AWF)


ggplot() +
  geom_sf(data = MeditCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

MeditCon<-as.data.frame(MeditCon)
MeditCon<-as.data.frame(MeditCon[,-which(colnames(MeditCon) == "geometry")])

write.csv(MeditCon, "Medit/MeditCon.csv", row.names = FALSE )






MonGMap_clipped<-read_sf("MonG/MonGMap_clipped.shp")
thistouches<-st_touches(MonGMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
MonGMap_clipped$cluster<-c$membership
length(unique(MonGMap_clipped$cluster))
st_write(MonGMap_clipped, "MonG/MonGMap_clust.shp", append=FALSE)
MonGMap_clipped$cluster<-as.factor(MonGMap_clipped$cluster)

ggplot() +
  geom_sf(data = MonGMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

MonGCon <- MonGMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(MonGCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

MonGCon$min_dist<-NA
for(i in 1:nrow(MonGCon)){
  print(i)
  MonGCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

MonGCon$CPL<-mean(MonGCon$min_dist)

MonGCon$SLC<-max(MonGCon$total_area)

MonGCon$MSC<-mean(MonGCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * MonGCon$total_area[i] *
        MonGCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(MonGCon), each=nrow(MonGCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
MonGCon$AWF<-flux_tot$total_flux
MonGCon$meanAWF<-mean(MonGCon$AWF)


ggplot() +
  geom_sf(data = MonGCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

MonGCon<-as.data.frame(MonGCon)
MonGCon<-as.data.frame(MonGCon[,-which(colnames(MonGCon) == "geometry")])

write.csv(MonGCon, "MonG/MonGCon.csv", row.names = FALSE )







TaiMap_clipped<-read_sf("Tai/TaiMap_clipped.shp")
thistouches<-st_touches(TaiMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TaiMap_clipped$cluster<-c$membership
length(unique(TaiMap_clipped$cluster))
st_write(TaiMap_clipped, "Tai/TaiMap_clust.shp", append=FALSE)
TaiMap_clipped$cluster<-as.factor(TaiMap_clipped$cluster)

ggplot() +
  geom_sf(data = TaiMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TaiCon <- TaiMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TaiCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TaiCon$min_dist<-NA
for(i in 1:nrow(TaiCon)){
  print(i)
  TaiCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TaiCon$CPL<-mean(TaiCon$min_dist)

TaiCon$SLC<-max(TaiCon$total_area)

TaiCon$MSC<-mean(TaiCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TaiCon$total_area[i] *
        TaiCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TaiCon), each=nrow(TaiCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TaiCon$AWF<-flux_tot$total_flux
TaiCon$meanAWF<-mean(TaiCon$AWF)


ggplot() +
  geom_sf(data = TaiCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TaiCon<-as.data.frame(TaiCon)
TaiCon<-as.data.frame(TaiCon[,-which(colnames(TaiCon) == "geometry")])

write.csv(TaiCon, "Tai/TaiCon.csv", row.names = FALSE )






TemBFMap_clipped<-read_sf("TemBF/TemBFMap_clipped.shp")
thistouches<-st_touches(TemBFMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TemBFMap_clipped$cluster<-c$membership
length(unique(TemBFMap_clipped$cluster))
st_write(TemBFMap_clipped, "TemBF/TemBFMap_clust.shp", append=FALSE)
TemBFMap_clipped$cluster<-as.factor(TemBFMap_clipped$cluster)

ggplot() +
  geom_sf(data = TemBFMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TemBFCon <- TemBFMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TemBFCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TemBFCon$min_dist<-NA
for(i in 1:nrow(TemBFCon)){
  print(i)
  TemBFCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TemBFCon$CPL<-mean(TemBFCon$min_dist)

TemBFCon$SLC<-max(TemBFCon$total_area)

TemBFCon$MSC<-mean(TemBFCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TemBFCon$total_area[i] *
        TemBFCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TemBFCon), each=nrow(TemBFCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TemBFCon$AWF<-flux_tot$total_flux
TemBFCon$meanAWF<-mean(TemBFCon$AWF)


ggplot() +
  geom_sf(data = TemBFCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TemBFCon<-as.data.frame(TemBFCon)
TemBFCon<-as.data.frame(TemBFCon[,-which(colnames(TemBFCon) == "geometry")])

write.csv(TemBFCon, "TemBF/TemBFCon.csv", row.names = FALSE )






TemCFMap_clipped<-read_sf("TemCF/TemCFMap_clipped.shp")
thistouches<-st_touches(TemCFMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TemCFMap_clipped$cluster<-c$membership
length(unique(TemCFMap_clipped$cluster))
st_write(TemCFMap_clipped, "TemCF/TemCFMap_clust.shp", append=FALSE)
TemCFMap_clipped$cluster<-as.factor(TemCFMap_clipped$cluster)

ggplot() +
  geom_sf(data = TemCFMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TemCFCon <- TemCFMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TemCFCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TemCFCon$min_dist<-NA
for(i in 1:nrow(TemCFCon)){
  print(i)
  TemCFCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TemCFCon$CPL<-mean(TemCFCon$min_dist)

TemCFCon$SLC<-max(TemCFCon$total_area)

TemCFCon$MSC<-mean(TemCFCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TemCFCon$total_area[i] *
        TemCFCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TemCFCon), each=nrow(TemCFCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TemCFCon$AWF<-flux_tot$total_flux
TemCFCon$meanAWF<-mean(TemCFCon$AWF)


ggplot() +
  geom_sf(data = TemCFCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TemCFCon<-as.data.frame(TemCFCon)
TemCFCon<-as.data.frame(TemCFCon[,-which(colnames(TemCFCon) == "geometry")])

write.csv(TemCFCon, "TemCF/TemCFCon.csv", row.names = FALSE )






TemGMap_clipped<-read_sf("TemG/TemGMap_clipped.shp")
thistouches<-st_touches(TemGMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TemGMap_clipped$cluster<-c$membership
length(unique(TemGMap_clipped$cluster))
st_write(TemGMap_clipped, "TemG/TemGMap_clust.shp", append=FALSE)
TemGMap_clipped$cluster<-as.factor(TemGMap_clipped$cluster)

ggplot() +
  geom_sf(data = TemGMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TemGCon <- TemGMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TemGCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TemGCon$min_dist<-NA
for(i in 1:nrow(TemGCon)){
  print(i)
  TemGCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TemGCon$CPL<-mean(TemGCon$min_dist)

TemGCon$SLC<-max(TemGCon$total_area)

TemGCon$MSC<-mean(TemGCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TemGCon$total_area[i] *
        TemGCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TemGCon), each=nrow(TemGCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TemGCon$AWF<-flux_tot$total_flux
TemGCon$meanAWF<-mean(TemGCon$AWF)


ggplot() +
  geom_sf(data = TemGCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TemGCon<-as.data.frame(TemGCon)
TemGCon<-as.data.frame(TemGCon[,-which(colnames(TemGCon) == "geometry")])

write.csv(TemGCon, "TemG/TemGCon.csv", row.names = FALSE )







TroCFMap_clipped<-read_sf("TroCF/TroCFMap_clipped.shp")
thistouches<-st_touches(TroCFMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TroCFMap_clipped$cluster<-c$membership
length(unique(TroCFMap_clipped$cluster))
st_write(TroCFMap_clipped, "TroCF/TroCFMap_clust.shp", append=FALSE)
TroCFMap_clipped$cluster<-as.factor(TroCFMap_clipped$cluster)

ggplot() +
  geom_sf(data = TroCFMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TroCFCon <- TroCFMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TroCFCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TroCFCon$min_dist<-NA
for(i in 1:nrow(TroCFCon)){
  print(i)
  TroCFCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TroCFCon$CPL<-mean(TroCFCon$min_dist)

TroCFCon$SLC<-max(TroCFCon$total_area)

TroCFCon$MSC<-mean(TroCFCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TroCFCon$total_area[i] *
        TroCFCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TroCFCon), each=nrow(TroCFCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TroCFCon$AWF<-flux_tot$total_flux
TroCFCon$meanAWF<-mean(TroCFCon$AWF)


ggplot() +
  geom_sf(data = TroCFCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TroCFCon<-as.data.frame(TroCFCon)
TroCFCon<-as.data.frame(TroCFCon[,-which(colnames(TroCFCon) == "geometry")])

write.csv(TroCFCon, "TroCF/TroCFCon.csv", row.names = FALSE )







TroGMap_clipped<-read_sf("TroG/TroGMap_clipped.shp")
thistouches<-st_touches(TroGMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TroGMap_clipped$cluster<-c$membership
length(unique(TroGMap_clipped$cluster))
st_write(TroGMap_clipped, "TroG/TroGMap_clust.shp", append=FALSE)
TroGMap_clipped$cluster<-as.factor(TroGMap_clipped$cluster)

ggplot() +
  geom_sf(data = TroGMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TroGCon <- TroGMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TroGCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TroGCon$min_dist<-NA
for(i in 1:nrow(TroGCon)){
  print(i)
  TroGCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TroGCon$CPL<-mean(TroGCon$min_dist)

TroGCon$SLC<-max(TroGCon$total_area)

TroGCon$MSC<-mean(TroGCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TroGCon$total_area[i] *
        TroGCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TroGCon), each=nrow(TroGCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TroGCon$AWF<-flux_tot$total_flux
TroGCon$meanAWF<-mean(TroGCon$AWF)


ggplot() +
  geom_sf(data = TroGCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TroGCon<-as.data.frame(TroGCon)
TroGCon<-as.data.frame(TroGCon[,-which(colnames(TroGCon) == "geometry")])

write.csv(TroGCon, "TroG/TroGCon.csv", row.names = FALSE )







TunMap_clipped<-read_sf("Tun/TunMap_clipped.shp")
thistouches<-st_touches(TunMap_clipped)

g<-graph.adjlist(thistouches)
c<-components(g)
TunMap_clipped$cluster<-c$membership
length(unique(TunMap_clipped$cluster))
st_write(TunMap_clipped, "Tun/TunMap_clust.shp", append=FALSE)
TunMap_clipped$cluster<-as.factor(TunMap_clipped$cluster)

ggplot() +
  geom_sf(data = TunMap_clipped, aes(fill = as.factor(cluster))) +
  labs(fill = "Cluster") +
  theme_minimal()

TunCon <- TunMap_clipped %>%
  mutate(area = st_area(geometry)) %>%
  group_by(cluster) %>%
  summarise(total_area = sum(area)) %>%
  ungroup()

min_dist<-st_distance(TunCon$geometry)
min_dist<-as.data.frame(as.matrix(drop_units(min_dist)))
min_dist[min_dist == 0]<-NA
str(min_dist)

TunCon$min_dist<-NA
for(i in 1:nrow(TunCon)){
  print(i)
  TunCon$min_dist[i]<-min(min_dist[i,], na.rm = TRUE)
}

TunCon$CPL<-mean(TunCon$min_dist)

TunCon$SLC<-max(TunCon$total_area)

TunCon$MSC<-mean(TunCon$total_area)

k<-(-log(0.5))/(max(min_dist, na.rm = TRUE)/2)

AWF_mat<-matrix(NA, nrow = nrow(min_dist), ncol = ncol(min_dist))
for(i in 1:nrow(min_dist)){
  print(i)
  for(j in 1:nrow(min_dist)){
    if (i != j) {
      AWF_mat[i,j]<-exp(-k*(as.matrix(min_dist)[i,j])) * TunCon$total_area[i] *
        TunCon$total_area[j]
    }
  }
}

colnames(AWF_mat)<-1:nrow(AWF_mat)
AWF_mat<-as.data.frame(AWF_mat) %>%
  pivot_longer(cols = everything(), names_to =  "patch_j", values_to = "flux")
AWF_mat<-AWF_mat %>%
  mutate(patch_i = rep(1:nrow(TunCon), each=nrow(TunCon)))
flux_tot<- AWF_mat %>%
  group_by(patch_i) %>%
  summarise(total_flux = sum(flux, na.rm=TRUE))
TunCon$AWF<-flux_tot$total_flux
TunCon$meanAWF<-mean(TunCon$AWF)


ggplot() +
  geom_sf(data = TunCon, aes(fill = min_dist)) +
  labs(fill = "Cluster") +
  theme_minimal()

TunCon<-as.data.frame(TunCon)
TunCon<-as.data.frame(TunCon[,-which(colnames(TunCon) == "geometry")])

write.csv(TunCon, "Tun/TunCon.csv", row.names = FALSE )





## Sorensen vs Simpson beta diversities ----

# Code for the beta diversity calculations and plotting is unnecessarily long
# and makes way too many objects but it works

# Format matrix and phylogeny for each biome beginning with DBF
row.names(DBFMat)<-DBFMat$grids
DBFMat<-DBFMat[,-"grids"]
DBFMat<-as.matrix(as.data.frame(lapply(DBFMat, as.numeric)))
DBFMat<-dense2sparse(DBFMat)

DBFPhy$tip.label<-gsub("_", ".", DBFPhy$tip.label)

DBFMat<-DBFMat[,which(colnames(DBFMat)%in%DBFPhy$tip.label)] # double check, this shouldn't change DBFMat

#DBF Sor + Sim comparison using betapart
DBF_taxbeta<-beta.pair(DBFMat, index.family = "sorensen")
DBF_taxsor<-as.data.frame(as.matrix(DBF_taxbeta$beta.sor))
DBF_taxsor$grids<-row.names(DBF_taxsor)
write.csv(DBF_taxsor, "DBF/DBF_taxsor.csv", row.names = FALSE)
DBF_taxsim<-as.data.frame(as.matrix(DBF_taxbeta$beta.sim))
DBF_taxsim$grids<-row.names(DBF_taxsim)
write.csv(DBF_taxsim, "DBF/DBF_taxsim.csv", row.names = FALSE)

DBF_phylobeta<-phylobeta(DBFMat, DBFPhy, index.family = "sorensen")
DBF_physor<-as.data.frame(as.matrix(DBF_phylobeta$phylo.beta.sor))
DBF_physor$grids<-row.names(DBF_physor)
write.csv(DBF_physor, "DBF/DBF_physor.csv", row.names = FALSE)
DBF_physim<-as.data.frame(as.matrix(DBF_phylobeta$phylo.beta.sim))
DBF_physim$grids<-row.names(DBF_physim)
write.csv(DBF_physim, "DBF/DBF_physim.csv", row.names = FALSE)


#### Repeat this for the other 13 biomes


row.names(DesMat)<-DesMat$grids
DesMat<-DesMat[,-"grids"]
DesMat<-as.matrix(as.data.frame(lapply(DesMat, as.numeric)))
DesMat<-dense2sparse(DesMat)

DesPhy$tip.label<-gsub("_", ".", DesPhy$tip.label)

DesMat<-DesMat[,which(colnames(DesMat)%in%DesPhy$tip.label)]

#Des Sor + Sim comparison using betapart
Des_taxbeta<-beta.pair(DesMat, index.family = "sorensen")
Des_taxsor<-as.data.frame(as.matrix(Des_taxbeta$beta.sor))
Des_taxsor$grids<-row.names(Des_taxsor)
write.csv(Des_taxsor, "Des/Des_taxsor.csv", row.names = FALSE)
Des_taxsim<-as.data.frame(as.matrix(Des_taxbeta$beta.sim))
Des_taxsim$grids<-row.names(Des_taxsim)
write.csv(Des_taxsim, "Des/Des_taxsim.csv", row.names = FALSE)

Des_phylobeta<-phylobeta(DesMat, DesPhy, index.family = "sorensen")
Des_physor<-as.data.frame(as.matrix(Des_phylobeta$phylo.beta.sor))
Des_physor$grids<-row.names(Des_physor)
write.csv(Des_physor, "Des/Des_physor.csv", row.names = FALSE)
Des_physim<-as.data.frame(as.matrix(Des_phylobeta$phylo.beta.sim))
Des_physim$grids<-row.names(Des_physim)
write.csv(Des_physim, "Des/Des_physim.csv", row.names = FALSE)




row.names(FloGMat)<-FloGMat$grids
FloGMat<-FloGMat[,-"grids"]
FloGMat<-as.matrix(as.data.frame(lapply(FloGMat, as.numeric)))
FloGMat<-dense2sparse(FloGMat)

FloGPhy$tip.label<-gsub("_", ".", FloGPhy$tip.label)

FloGMat<-FloGMat[,which(colnames(FloGMat)%in%FloGPhy$tip.label)]

#FloG Sor + Sim comparison using betapart
FloG_taxbeta<-beta.pair(FloGMat, index.family = "sorensen")
FloG_taxsor<-as.data.frame(as.matrix(FloG_taxbeta$beta.sor))
FloG_taxsor$grids<-row.names(FloG_taxsor)
write.csv(FloG_taxsor, "FloG/FloG_taxsor.csv", row.names = FALSE)
FloG_taxsim<-as.data.frame(as.matrix(FloG_taxbeta$beta.sim))
FloG_taxsim$grids<-row.names(FloG_taxsim)
write.csv(FloG_taxsim, "FloG/FloG_taxsim.csv", row.names = FALSE)

FloG_phylobeta<-phylobeta(FloGMat, FloGPhy, index.family = "sorensen")
FloG_physor<-as.data.frame(as.matrix(FloG_phylobeta$phylo.beta.sor))
FloG_physor$grids<-row.names(FloG_physor)
write.csv(FloG_physor, "FloG/FloG_physor.csv", row.names = FALSE)
FloG_physim<-as.data.frame(as.matrix(FloG_phylobeta$phylo.beta.sim))
FloG_physim$grids<-row.names(FloG_physim)
write.csv(FloG_physim, "FloG/FloG_physim.csv", row.names = FALSE)







row.names(MangMat)<-MangMat$grids
MangMat<-MangMat[,-"grids"]
MangMat<-as.matrix(as.data.frame(lapply(MangMat, as.numeric)))
MangMat<-dense2sparse(MangMat)

MangPhy$tip.label<-gsub("_", ".", MangPhy$tip.label)

MangMat<-MangMat[,which(colnames(MangMat)%in%MangPhy$tip.label)]

#Mang Sor + Sim comparison using betapart
Mang_taxbeta<-beta.pair(MangMat, index.family = "sorensen")
Mang_taxsor<-as.data.frame(as.matrix(Mang_taxbeta$beta.sor))
Mang_taxsor$grids<-row.names(Mang_taxsor)
write.csv(Mang_taxsor, "Mang/Mang_taxsor.csv", row.names = FALSE)
Mang_taxsim<-as.data.frame(as.matrix(Mang_taxbeta$beta.sim))
Mang_taxsim$grids<-row.names(Mang_taxsim)
write.csv(Mang_taxsim, "Mang/Mang_taxsim.csv", row.names = FALSE)

Mang_phylobeta<-phylobeta(MangMat, MangPhy, index.family = "sorensen")
Mang_physor<-as.data.frame(as.matrix(Mang_phylobeta$phylo.beta.sor))
Mang_physor$grids<-row.names(Mang_physor)
write.csv(Mang_physor, "Mang/Mang_physor.csv", row.names = FALSE)
Mang_physim<-as.data.frame(as.matrix(Mang_phylobeta$phylo.beta.sim))
Mang_physim$grids<-row.names(Mang_physim)
write.csv(Mang_physim, "Mang/Mang_physim.csv", row.names = FALSE)







row.names(MBFMat)<-MBFMat$grids
MBFMat<-MBFMat[,-"grids"]
MBFMat<-as.matrix(as.data.frame(lapply(MBFMat, as.numeric)))
MBFMat<-dense2sparse(MBFMat)

MBFPhy$tip.label<-gsub("_", ".", MBFPhy$tip.label)

MBFMat<-MBFMat[,which(colnames(MBFMat)%in%MBFPhy$tip.label)]

#MBF Sor + Sim comparison using betapart
MBF_taxbeta<-beta.pair(MBFMat, index.family = "sorensen")
MBF_taxsor<-as.data.frame(as.matrix(MBF_taxbeta$beta.sor))
MBF_taxsor$grids<-row.names(MBF_taxsor)
write.csv(MBF_taxsor, "MBF/MBF_taxsor.csv", row.names = FALSE)
MBF_taxsim<-as.data.frame(as.matrix(MBF_taxbeta$beta.sim))
MBF_taxsim$grids<-row.names(MBF_taxsim)
write.csv(MBF_taxsim, "MBF/MBF_taxsim.csv", row.names = FALSE)

MBF_phylobeta<-phylobeta(MBFMat, MBFPhy, index.family = "sorensen")
MBF_physor<-as.data.frame(as.matrix(MBF_phylobeta$phylo.beta.sor))
MBF_physor$grids<-row.names(MBF_physor)
write.csv(MBF_physor, "MBF/MBF_physor.csv", row.names = FALSE)
MBF_physim<-as.data.frame(as.matrix(MBF_phylobeta$phylo.beta.sim))
MBF_physim$grids<-row.names(MBF_physim)
write.csv(MBF_physim, "MBF/MBF_physim.csv", row.names = FALSE)







row.names(MeditMat)<-MeditMat$grids
MeditMat<-MeditMat[,-"grids"]
MeditMat<-as.matrix(as.data.frame(lapply(MeditMat, as.numeric)))
MeditMat<-dense2sparse(MeditMat)

MeditPhy$tip.label<-gsub("_", ".", MeditPhy$tip.label)

MeditMat<-MeditMat[,which(colnames(MeditMat)%in%MeditPhy$tip.label)]

#Medit Sor + Sim comparison using betapart
Medit_taxbeta<-beta.pair(MeditMat, index.family = "sorensen")
Medit_taxsor<-as.data.frame(as.matrix(Medit_taxbeta$beta.sor))
Medit_taxsor$grids<-row.names(Medit_taxsor)
write.csv(Medit_taxsor, "Medit/Medit_taxsor.csv", row.names = FALSE)
Medit_taxsim<-as.data.frame(as.matrix(Medit_taxbeta$beta.sim))
Medit_taxsim$grids<-row.names(Medit_taxsim)
write.csv(Medit_taxsim, "Medit/Medit_taxsim.csv", row.names = FALSE)

Medit_phylobeta<-phylobeta(MeditMat, MeditPhy, index.family = "sorensen")
Medit_physor<-as.data.frame(as.matrix(Medit_phylobeta$phylo.beta.sor))
Medit_physor$grids<-row.names(Medit_physor)
write.csv(Medit_physor, "Medit/Medit_physor.csv", row.names = FALSE)
Medit_physim<-as.data.frame(as.matrix(Medit_phylobeta$phylo.beta.sim))
Medit_physim$grids<-row.names(Medit_physim)
write.csv(Medit_physim, "Medit/Medit_physim.csv", row.names = FALSE)






row.names(MonGMat)<-MonGMat$grids
MonGMat<-MonGMat[,-"grids"]
MonGMat<-as.matrix(as.data.frame(lapply(MonGMat, as.numeric)))
MonGMat<-dense2sparse(MonGMat)

MonGPhy$tip.label<-gsub("_", ".", MonGPhy$tip.label)

MonGMat<-MonGMat[,which(colnames(MonGMat)%in%MonGPhy$tip.label)]

#MonG Sor + Sim comparison using betapart
MonG_taxbeta<-beta.pair(MonGMat, index.family = "sorensen")
MonG_taxsor<-as.data.frame(as.matrix(MonG_taxbeta$beta.sor))
MonG_taxsor$grids<-row.names(MonG_taxsor)
write.csv(MonG_taxsor, "MonG/MonG_taxsor.csv", row.names = FALSE)
MonG_taxsim<-as.data.frame(as.matrix(MonG_taxbeta$beta.sim))
MonG_taxsim$grids<-row.names(MonG_taxsim)
write.csv(MonG_taxsim, "MonG/MonG_taxsim.csv", row.names = FALSE)

MonG_phylobeta<-phylobeta(MonGMat, MonGPhy, index.family = "sorensen")
MonG_physor<-as.data.frame(as.matrix(MonG_phylobeta$phylo.beta.sor))
MonG_physor$grids<-row.names(MonG_physor)
write.csv(MonG_physor, "MonG/MonG_physor.csv", row.names = FALSE)
MonG_physim<-as.data.frame(as.matrix(MonG_phylobeta$phylo.beta.sim))
MonG_physim$grids<-row.names(MonG_physim)
write.csv(MonG_physim, "MonG/MonG_physim.csv", row.names = FALSE)





row.names(TaiMat)<-TaiMat$grids
TaiMat<-TaiMat[,-"grids"]
TaiMat<-as.matrix(as.data.frame(lapply(TaiMat, as.numeric)))
TaiMat<-dense2sparse(TaiMat)

TaiPhy$tip.label<-gsub("_", ".", TaiPhy$tip.label)

TaiMat<-TaiMat[,which(colnames(TaiMat)%in%TaiPhy$tip.label)]

#Tai Sor + Sim comparison using betapart
Tai_taxbeta<-beta.pair(TaiMat, index.family = "sorensen")
Tai_taxsor<-as.data.frame(as.matrix(Tai_taxbeta$beta.sor))
Tai_taxsor$grids<-row.names(Tai_taxsor)
write.csv(Tai_taxsor, "Tai/Tai_taxsor.csv", row.names = FALSE)
Tai_taxsim<-as.data.frame(as.matrix(Tai_taxbeta$beta.sim))
Tai_taxsim$grids<-row.names(Tai_taxsim)
write.csv(Tai_taxsim, "Tai/Tai_taxsim.csv", row.names = FALSE)

Tai_phylobeta<-phylobeta(TaiMat, TaiPhy, index.family = "sorensen")
Tai_physor<-as.data.frame(as.matrix(Tai_phylobeta$phylo.beta.sor))
Tai_physor$grids<-row.names(Tai_physor)
write.csv(Tai_physor, "Tai/Tai_physor.csv", row.names = FALSE)
Tai_physim<-as.data.frame(as.matrix(Tai_phylobeta$phylo.beta.sim))
Tai_physim$grids<-row.names(Tai_physim)
write.csv(Tai_physim, "Tai/Tai_physim.csv", row.names = FALSE)




row.names(TemBFMat)<-TemBFMat$grids
TemBFMat<-TemBFMat[,-"grids"]
TemBFMat<-as.matrix(as.data.frame(lapply(TemBFMat, as.numeric)))
TemBFMat<-dense2sparse(TemBFMat)

TemBFPhy$tip.label<-gsub("_", ".", TemBFPhy$tip.label)

TemBFMat<-TemBFMat[,which(colnames(TemBFMat)%in%TemBFPhy$tip.label)]

#TemBF Sor + Sim comparison using betapart
TemBF_taxbeta<-beta.pair(TemBFMat, index.family = "sorensen")
TemBF_taxsor<-as.data.frame(as.matrix(TemBF_taxbeta$beta.sor))
TemBF_taxsor$grids<-row.names(TemBF_taxsor)
write.csv(TemBF_taxsor, "TemBF/TemBF_taxsor.csv", row.names = FALSE)
TemBF_taxsim<-as.data.frame(as.matrix(TemBF_taxbeta$beta.sim))
TemBF_taxsim$grids<-row.names(TemBF_taxsim)
write.csv(TemBF_taxsim, "TemBF/TemBF_taxsim.csv", row.names = FALSE)

TemBF_phylobeta<-phylobeta(TemBFMat, TemBFPhy, index.family = "sorensen")
TemBF_physor<-as.data.frame(as.matrix(TemBF_phylobeta$phylo.beta.sor))
TemBF_physor$grids<-row.names(TemBF_physor)
write.csv(TemBF_physor, "TemBF/TemBF_physor.csv", row.names = FALSE)
TemBF_physim<-as.data.frame(as.matrix(TemBF_phylobeta$phylo.beta.sim))
TemBF_physim$grids<-row.names(TemBF_physim)
write.csv(TemBF_physim, "TemBF/TemBF_physim.csv", row.names = FALSE)






row.names(TemCFMat)<-TemCFMat$grids
TemCFMat<-TemCFMat[,-"grids"]
TemCFMat<-as.matrix(as.data.frame(lapply(TemCFMat, as.numeric)))
TemCFMat<-dense2sparse(TemCFMat)

TemCFPhy$tip.label<-gsub("_", ".", TemCFPhy$tip.label)

TemCFMat<-TemCFMat[,which(colnames(TemCFMat)%in%TemCFPhy$tip.label)]

#TemCF Sor + Sim comparison using betapart
TemCF_taxbeta<-beta.pair(TemCFMat, index.family = "sorensen")
TemCF_taxsor<-as.data.frame(as.matrix(TemCF_taxbeta$beta.sor))
TemCF_taxsor$grids<-row.names(TemCF_taxsor)
write.csv(TemCF_taxsor, "TemCF/TemCF_taxsor.csv", row.names = FALSE)
TemCF_taxsim<-as.data.frame(as.matrix(TemCF_taxbeta$beta.sim))
TemCF_taxsim$grids<-row.names(TemCF_taxsim)
write.csv(TemCF_taxsim, "TemCF/TemCF_taxsim.csv", row.names = FALSE)

TemCF_phylobeta<-phylobeta(TemCFMat, TemCFPhy, index.family = "sorensen")
TemCF_physor<-as.data.frame(as.matrix(TemCF_phylobeta$phylo.beta.sor))
TemCF_physor$grids<-row.names(TemCF_physor)
write.csv(TemCF_physor, "TemCF/TemCF_physor.csv", row.names = FALSE)
TemCF_physim<-as.data.frame(as.matrix(TemCF_phylobeta$phylo.beta.sim))
TemCF_physim$grids<-row.names(TemCF_physim)
write.csv(TemCF_physim, "TemCF/TemCF_physim.csv", row.names = FALSE)




row.names(TemGMat)<-TemGMat$grids
TemGMat<-TemGMat[,-"grids"]
TemGMat<-as.matrix(as.data.frame(lapply(TemGMat, as.numeric)))
TemGMat<-dense2sparse(TemGMat)

TemGPhy$tip.label<-gsub("_", ".", TemGPhy$tip.label)

TemGMat<-TemGMat[,which(colnames(TemGMat)%in%TemGPhy$tip.label)]

#TemG Sor + Sim comparison using betapart
TemG_taxbeta<-beta.pair(TemGMat, index.family = "sorensen")
TemG_taxsor<-as.data.frame(as.matrix(TemG_taxbeta$beta.sor))
TemG_taxsor$grids<-row.names(TemG_taxsor)
write.csv(TemG_taxsor, "TemG/TemG_taxsor.csv", row.names = FALSE)
TemG_taxsim<-as.data.frame(as.matrix(TemG_taxbeta$beta.sim))
TemG_taxsim$grids<-row.names(TemG_taxsim)
write.csv(TemG_taxsim, "TemG/TemG_taxsim.csv", row.names = FALSE)

TemG_phylobeta<-phylobeta(TemGMat, TemGPhy, index.family = "sorensen")
TemG_physor<-as.data.frame(as.matrix(TemG_phylobeta$phylo.beta.sor))
TemG_physor$grids<-row.names(TemG_physor)
write.csv(TemG_physor, "TemG/TemG_physor.csv", row.names = FALSE)
TemG_physim<-as.data.frame(as.matrix(TemG_phylobeta$phylo.beta.sim))
TemG_physim$grids<-row.names(TemG_physim)
write.csv(TemG_physim, "TemG/TemG_physim.csv", row.names = FALSE)




row.names(TroCFMat)<-TroCFMat$grids
TroCFMat<-TroCFMat[,-"grids"]
TroCFMat<-as.matrix(as.data.frame(lapply(TroCFMat, as.numeric)))
TroCFMat<-dense2sparse(TroCFMat)

TroCFPhy$tip.label<-gsub("_", ".", TroCFPhy$tip.label)

TroCFMat<-TroCFMat[,which(colnames(TroCFMat)%in%TroCFPhy$tip.label)]

#TroCF Sor + Sim comparison using betapart
TroCF_taxbeta<-beta.pair(TroCFMat, index.family = "sorensen")
TroCF_taxsor<-as.data.frame(as.matrix(TroCF_taxbeta$beta.sor))
TroCF_taxsor$grids<-row.names(TroCF_taxsor)
write.csv(TroCF_taxsor, "TroCF/TroCF_taxsor.csv", row.names = FALSE)
TroCF_taxsim<-as.data.frame(as.matrix(TroCF_taxbeta$beta.sim))
TroCF_taxsim$grids<-row.names(TroCF_taxsim)
write.csv(TroCF_taxsim, "TroCF/TroCF_taxsim.csv", row.names = FALSE)

TroCF_phylobeta<-phylobeta(TroCFMat, TroCFPhy, index.family = "sorensen")
TroCF_physor<-as.data.frame(as.matrix(TroCF_phylobeta$phylo.beta.sor))
TroCF_physor$grids<-row.names(TroCF_physor)
write.csv(TroCF_physor, "TroCF/TroCF_physor.csv", row.names = FALSE)
TroCF_physim<-as.data.frame(as.matrix(TroCF_phylobeta$phylo.beta.sim))
TroCF_physim$grids<-row.names(TroCF_physim)
write.csv(TroCF_physim, "TroCF/TroCF_physim.csv", row.names = FALSE)




row.names(TroGMat)<-TroGMat$grids
TroGMat<-TroGMat[,-"grids"]
TroGMat<-as.matrix(as.data.frame(lapply(TroGMat, as.numeric)))
TroGMat<-dense2sparse(TroGMat)

TroGPhy$tip.label<-gsub("_", ".", TroGPhy$tip.label)

TroGMat<-TroGMat[,which(colnames(TroGMat)%in%TroGPhy$tip.label)]

#TroG Sor + Sim comparison using betapart
TroG_taxbeta<-beta.pair(TroGMat, index.family = "sorensen")
TroG_taxsor<-as.data.frame(as.matrix(TroG_taxbeta$beta.sor))
TroG_taxsor$grids<-row.names(TroG_taxsor)
write.csv(TroG_taxsor, "TroG/TroG_taxsor.csv", row.names = FALSE)
TroG_taxsim<-as.data.frame(as.matrix(TroG_taxbeta$beta.sim))
TroG_taxsim$grids<-row.names(TroG_taxsim)
write.csv(TroG_taxsim, "TroG/TroG_taxsim.csv", row.names = FALSE)

TroG_phylobeta<-phylobeta(TroGMat, TroGPhy, index.family = "sorensen")
TroG_physor<-as.data.frame(as.matrix(TroG_phylobeta$phylo.beta.sor))
TroG_physor$grids<-row.names(TroG_physor)
write.csv(TroG_physor, "TroG/TroG_physor.csv", row.names = FALSE)
TroG_physim<-as.data.frame(as.matrix(TroG_phylobeta$phylo.beta.sim))
TroG_physim$grids<-row.names(TroG_physim)
write.csv(TroG_physim, "TroG/TroG_physim.csv", row.names = FALSE)




row.names(TunMat)<-TunMat$grids
TunMat<-TunMat[,-"grids"]
TunMat<-as.matrix(as.data.frame(lapply(TunMat, as.numeric)))
TunMat<-dense2sparse(TunMat)

TunPhy$tip.label<-gsub("_", ".", TunPhy$tip.label)

TunMat<-TunMat[,which(colnames(TunMat)%in%TunPhy$tip.label)]

#Tun Sor + Sim comparison using betapart
Tun_taxbeta<-beta.pair(TunMat, index.family = "sorensen")
Tun_taxsor<-as.data.frame(as.matrix(Tun_taxbeta$beta.sor))
Tun_taxsor$grids<-row.names(Tun_taxsor)
write.csv(Tun_taxsor, "Tun/Tun_taxsor.csv", row.names = FALSE)
Tun_taxsim<-as.data.frame(as.matrix(Tun_taxbeta$beta.sim))
Tun_taxsim$grids<-row.names(Tun_taxsim)
write.csv(Tun_taxsim, "Tun/Tun_taxsim.csv", row.names = FALSE)

Tun_phylobeta<-phylobeta(TunMat, TunPhy, index.family = "sorensen")
Tun_physor<-as.data.frame(as.matrix(Tun_phylobeta$phylo.beta.sor))
Tun_physor$grids<-row.names(Tun_physor)
write.csv(Tun_physor, "Tun/Tun_physor.csv", row.names = FALSE)
Tun_physim<-as.data.frame(as.matrix(Tun_phylobeta$phylo.beta.sim))
Tun_physim$grids<-row.names(Tun_physim)
write.csv(Tun_physim, "Tun/Tun_physim.csv", row.names = FALSE)


## Plot Sor vs Sim beta diversities ----

# Format the turnover matrices to only keep the upper portion of each

DBF_taxsorUpper<-DBF_taxsor[,-which(colnames(DBF_taxsor) == "grids")]
DBF_taxsorUpper<-as.data.frame(DBF_taxsorUpper[upper.tri(DBF_taxsorUpper, diag = FALSE)])
DBF_taxsorUpper$biomeabbr<-"DBF"
colnames(DBF_taxsorUpper)[1]<-"taxsor"
DBF_taxsimUpper<-DBF_taxsim[,-which(colnames(DBF_taxsim) == "grids")]
DBF_taxsimUpper<-as.data.frame(DBF_taxsimUpper[upper.tri(DBF_taxsimUpper, diag = FALSE)])
DBF_taxsimUpper$biomeabbr<-"DBF"
colnames(DBF_taxsimUpper)[1]<-"taxsim"
DBF_physorUpper<-DBF_physor[,-which(colnames(DBF_physor) == "grids")]
DBF_physorUpper<-as.data.frame(DBF_physorUpper[upper.tri(DBF_physorUpper, diag = FALSE)])
DBF_physorUpper$biomeabbr<-"DBF"
colnames(DBF_physorUpper)[1]<-"physor"
DBF_physimUpper<-DBF_physim[,-which(colnames(DBF_physim) == "grids")]
DBF_physimUpper<-as.data.frame(DBF_physimUpper[upper.tri(DBF_physimUpper, diag = FALSE)])
DBF_physimUpper$biomeabbr<-"DBF"
colnames(DBF_physimUpper)[1]<-"physim"


Des_taxsorUpper<-Des_taxsor[,-which(colnames(Des_taxsor) == "grids")]
Des_taxsorUpper<-as.data.frame(Des_taxsorUpper[upper.tri(Des_taxsorUpper, diag = FALSE)])
Des_taxsorUpper$biomeabbr<-"Des"
colnames(Des_taxsorUpper)[1]<-"taxsor"
Des_taxsimUpper<-Des_taxsim[,-which(colnames(Des_taxsim) == "grids")]
Des_taxsimUpper<-as.data.frame(Des_taxsimUpper[upper.tri(Des_taxsimUpper, diag = FALSE)])
Des_taxsimUpper$biomeabbr<-"Des"
colnames(Des_taxsimUpper)[1]<-"taxsim"
Des_physorUpper<-Des_physor[,-which(colnames(Des_physor) == "grids")]
Des_physorUpper<-as.data.frame(Des_physorUpper[upper.tri(Des_physorUpper, diag = FALSE)])
Des_physorUpper$biomeabbr<-"Des"
colnames(Des_physorUpper)[1]<-"physor"
Des_physimUpper<-Des_physim[,-which(colnames(Des_physim) == "grids")]
Des_physimUpper<-as.data.frame(Des_physimUpper[upper.tri(Des_physimUpper, diag = FALSE)])
Des_physimUpper$biomeabbr<-"Des"
colnames(Des_physimUpper)[1]<-"physim"


FloG_taxsorUpper<-FloG_taxsor[,-which(colnames(FloG_taxsor) == "grids")]
FloG_taxsorUpper<-as.data.frame(FloG_taxsorUpper[upper.tri(FloG_taxsorUpper, diag = FALSE)])
FloG_taxsorUpper$biomeabbr<-"FloG"
colnames(FloG_taxsorUpper)[1]<-"taxsor"
FloG_taxsimUpper<-FloG_taxsim[,-which(colnames(FloG_taxsim) == "grids")]
FloG_taxsimUpper<-as.data.frame(FloG_taxsimUpper[upper.tri(FloG_taxsimUpper, diag = FALSE)])
FloG_taxsimUpper$biomeabbr<-"FloG"
colnames(FloG_taxsimUpper)[1]<-"taxsim"
FloG_physorUpper<-FloG_physor[,-which(colnames(FloG_physor) == "grids")]
FloG_physorUpper<-as.data.frame(FloG_physorUpper[upper.tri(FloG_physorUpper, diag = FALSE)])
FloG_physorUpper$biomeabbr<-"FloG"
colnames(FloG_physorUpper)[1]<-"physor"
FloG_physimUpper<-FloG_physim[,-which(colnames(FloG_physim) == "grids")]
FloG_physimUpper<-as.data.frame(FloG_physimUpper[upper.tri(FloG_physimUpper, diag = FALSE)])
FloG_physimUpper$biomeabbr<-"FloG"
colnames(FloG_physimUpper)[1]<-"physim"


Mang_taxsorUpper<-Mang_taxsor[,-which(colnames(Mang_taxsor) == "grids")]
Mang_taxsorUpper<-as.data.frame(Mang_taxsorUpper[upper.tri(Mang_taxsorUpper, diag = FALSE)])
Mang_taxsorUpper$biomeabbr<-"Mang"
colnames(Mang_taxsorUpper)[1]<-"taxsor"
Mang_taxsimUpper<-Mang_taxsim[,-which(colnames(Mang_taxsim) == "grids")]
Mang_taxsimUpper<-as.data.frame(Mang_taxsimUpper[upper.tri(Mang_taxsimUpper, diag = FALSE)])
Mang_taxsimUpper$biomeabbr<-"Mang"
colnames(Mang_taxsimUpper)[1]<-"taxsim"
Mang_physorUpper<-Mang_physor[,-which(colnames(Mang_physor) == "grids")]
Mang_physorUpper<-as.data.frame(Mang_physorUpper[upper.tri(Mang_physorUpper, diag = FALSE)])
Mang_physorUpper$biomeabbr<-"Mang"
colnames(Mang_physorUpper)[1]<-"physor"
Mang_physimUpper<-Mang_physim[,-which(colnames(Mang_physim) == "grids")]
Mang_physimUpper<-as.data.frame(Mang_physimUpper[upper.tri(Mang_physimUpper, diag = FALSE)])
Mang_physimUpper$biomeabbr<-"Mang"
colnames(Mang_physimUpper)[1]<-"physim"


MBF_taxsorUpper<-MBF_taxsor[,-which(colnames(MBF_taxsor) == "grids")]
MBF_taxsorUpper<-as.data.frame(MBF_taxsorUpper[upper.tri(MBF_taxsorUpper, diag = FALSE)])
MBF_taxsorUpper$biomeabbr<-"MBF"
colnames(MBF_taxsorUpper)[1]<-"taxsor"
MBF_taxsimUpper<-MBF_taxsim[,-which(colnames(MBF_taxsim) == "grids")]
MBF_taxsimUpper<-as.data.frame(MBF_taxsimUpper[upper.tri(MBF_taxsimUpper, diag = FALSE)])
MBF_taxsimUpper$biomeabbr<-"MBF"
colnames(MBF_taxsimUpper)[1]<-"taxsim"
MBF_physorUpper<-MBF_physor[,-which(colnames(MBF_physor) == "grids")]
MBF_physorUpper<-as.data.frame(MBF_physorUpper[upper.tri(MBF_physorUpper, diag = FALSE)])
MBF_physorUpper$biomeabbr<-"MBF"
colnames(MBF_physorUpper)[1]<-"physor"
MBF_physimUpper<-MBF_physim[,-which(colnames(MBF_physim) == "grids")]
MBF_physimUpper<-as.data.frame(MBF_physimUpper[upper.tri(MBF_physimUpper, diag = FALSE)])
MBF_physimUpper$biomeabbr<-"MBF"
colnames(MBF_physimUpper)[1]<-"physim"



Medit_taxsorUpper<-Medit_taxsor[,-which(colnames(Medit_taxsor) == "grids")]
Medit_taxsorUpper<-as.data.frame(Medit_taxsorUpper[upper.tri(Medit_taxsorUpper, diag = FALSE)])
Medit_taxsorUpper$biomeabbr<-"Medit"
colnames(Medit_taxsorUpper)[1]<-"taxsor"
Medit_taxsimUpper<-Medit_taxsim[,-which(colnames(Medit_taxsim) == "grids")]
Medit_taxsimUpper<-as.data.frame(Medit_taxsimUpper[upper.tri(Medit_taxsimUpper, diag = FALSE)])
Medit_taxsimUpper$biomeabbr<-"Medit"
colnames(Medit_taxsimUpper)[1]<-"taxsim"
Medit_physorUpper<-Medit_physor[,-which(colnames(Medit_physor) == "grids")]
Medit_physorUpper<-as.data.frame(Medit_physorUpper[upper.tri(Medit_physorUpper, diag = FALSE)])
Medit_physorUpper$biomeabbr<-"Medit"
colnames(Medit_physorUpper)[1]<-"physor"
Medit_physimUpper<-Medit_physim[,-which(colnames(Medit_physim) == "grids")]
Medit_physimUpper<-as.data.frame(Medit_physimUpper[upper.tri(Medit_physimUpper, diag = FALSE)])
Medit_physimUpper$biomeabbr<-"Medit"
colnames(Medit_physimUpper)[1]<-"physim"


MonG_taxsorUpper<-MonG_taxsor[,-which(colnames(MonG_taxsor) == "grids")]
MonG_taxsorUpper<-as.data.frame(MonG_taxsorUpper[upper.tri(MonG_taxsorUpper, diag = FALSE)])
MonG_taxsorUpper$biomeabbr<-"MonG"
colnames(MonG_taxsorUpper)[1]<-"taxsor"
MonG_taxsimUpper<-MonG_taxsim[,-which(colnames(MonG_taxsim) == "grids")]
MonG_taxsimUpper<-as.data.frame(MonG_taxsimUpper[upper.tri(MonG_taxsimUpper, diag = FALSE)])
MonG_taxsimUpper$biomeabbr<-"MonG"
colnames(MonG_taxsimUpper)[1]<-"taxsim"
MonG_physorUpper<-MonG_physor[,-which(colnames(MonG_physor) == "grids")]
MonG_physorUpper<-as.data.frame(MonG_physorUpper[upper.tri(MonG_physorUpper, diag = FALSE)])
MonG_physorUpper$biomeabbr<-"MonG"
colnames(MonG_physorUpper)[1]<-"physor"
MonG_physimUpper<-MonG_physim[,-which(colnames(MonG_physim) == "grids")]
MonG_physimUpper<-as.data.frame(MonG_physimUpper[upper.tri(MonG_physimUpper, diag = FALSE)])
MonG_physimUpper$biomeabbr<-"MonG"
colnames(MonG_physimUpper)[1]<-"physim"


Tai_taxsorUpper<-Tai_taxsor[,-which(colnames(Tai_taxsor) == "grids")]
Tai_taxsorUpper<-as.data.frame(Tai_taxsorUpper[upper.tri(Tai_taxsorUpper, diag = FALSE)])
Tai_taxsorUpper$biomeabbr<-"Tai"
colnames(Tai_taxsorUpper)[1]<-"taxsor"
Tai_taxsimUpper<-Tai_taxsim[,-which(colnames(Tai_taxsim) == "grids")]
Tai_taxsimUpper<-as.data.frame(Tai_taxsimUpper[upper.tri(Tai_taxsimUpper, diag = FALSE)])
Tai_taxsimUpper$biomeabbr<-"Tai"
colnames(Tai_taxsimUpper)[1]<-"taxsim"
Tai_physorUpper<-Tai_physor[,-which(colnames(Tai_physor) == "grids")]
Tai_physorUpper<-as.data.frame(Tai_physorUpper[upper.tri(Tai_physorUpper, diag = FALSE)])
Tai_physorUpper$biomeabbr<-"Tai"
colnames(Tai_physorUpper)[1]<-"physor"
Tai_physimUpper<-Tai_physim[,-which(colnames(Tai_physim) == "grids")]
Tai_physimUpper<-as.data.frame(Tai_physimUpper[upper.tri(Tai_physimUpper, diag = FALSE)])
Tai_physimUpper$biomeabbr<-"Tai"
colnames(Tai_physimUpper)[1]<-"physim"


TemBF_taxsorUpper<-TemBF_taxsor[,-which(colnames(TemBF_taxsor) == "grids")]
TemBF_taxsorUpper<-as.data.frame(TemBF_taxsorUpper[upper.tri(TemBF_taxsorUpper, diag = FALSE)])
TemBF_taxsorUpper$biomeabbr<-"TemBF"
colnames(TemBF_taxsorUpper)[1]<-"taxsor"
TemBF_taxsimUpper<-TemBF_taxsim[,-which(colnames(TemBF_taxsim) == "grids")]
TemBF_taxsimUpper<-as.data.frame(TemBF_taxsimUpper[upper.tri(TemBF_taxsimUpper, diag = FALSE)])
TemBF_taxsimUpper$biomeabbr<-"TemBF"
colnames(TemBF_taxsimUpper)[1]<-"taxsim"
TemBF_physorUpper<-TemBF_physor[,-which(colnames(TemBF_physor) == "grids")]
TemBF_physorUpper<-as.data.frame(TemBF_physorUpper[upper.tri(TemBF_physorUpper, diag = FALSE)])
TemBF_physorUpper$biomeabbr<-"TemBF"
colnames(TemBF_physorUpper)[1]<-"physor"
TemBF_physimUpper<-TemBF_physim[,-which(colnames(TemBF_physim) == "grids")]
TemBF_physimUpper<-as.data.frame(TemBF_physimUpper[upper.tri(TemBF_physimUpper, diag = FALSE)])
TemBF_physimUpper$biomeabbr<-"TemBF"
colnames(TemBF_physimUpper)[1]<-"physim"


TemCF_taxsorUpper<-TemCF_taxsor[,-which(colnames(TemCF_taxsor) == "grids")]
TemCF_taxsorUpper<-as.data.frame(TemCF_taxsorUpper[upper.tri(TemCF_taxsorUpper, diag = FALSE)])
TemCF_taxsorUpper$biomeabbr<-"TemCF"
colnames(TemCF_taxsorUpper)[1]<-"taxsor"
TemCF_taxsimUpper<-TemCF_taxsim[,-which(colnames(TemCF_taxsim) == "grids")]
TemCF_taxsimUpper<-as.data.frame(TemCF_taxsimUpper[upper.tri(TemCF_taxsimUpper, diag = FALSE)])
TemCF_taxsimUpper$biomeabbr<-"TemCF"
colnames(TemCF_taxsimUpper)[1]<-"taxsim"
TemCF_physorUpper<-TemCF_physor[,-which(colnames(TemCF_physor) == "grids")]
TemCF_physorUpper<-as.data.frame(TemCF_physorUpper[upper.tri(TemCF_physorUpper, diag = FALSE)])
TemCF_physorUpper$biomeabbr<-"TemCF"
colnames(TemCF_physorUpper)[1]<-"physor"
TemCF_physimUpper<-TemCF_physim[,-which(colnames(TemCF_physim) == "grids")]
TemCF_physimUpper<-as.data.frame(TemCF_physimUpper[upper.tri(TemCF_physimUpper, diag = FALSE)])
TemCF_physimUpper$biomeabbr<-"TemCF"
colnames(TemCF_physimUpper)[1]<-"physim"


TemG_taxsorUpper<-TemG_taxsor[,-which(colnames(TemG_taxsor) == "grids")]
TemG_taxsorUpper<-as.data.frame(TemG_taxsorUpper[upper.tri(TemG_taxsorUpper, diag = FALSE)])
TemG_taxsorUpper$biomeabbr<-"TemG"
colnames(TemG_taxsorUpper)[1]<-"taxsor"
TemG_taxsimUpper<-TemG_taxsim[,-which(colnames(TemG_taxsim) == "grids")]
TemG_taxsimUpper<-as.data.frame(TemG_taxsimUpper[upper.tri(TemG_taxsimUpper, diag = FALSE)])
TemG_taxsimUpper$biomeabbr<-"TemG"
colnames(TemG_taxsimUpper)[1]<-"taxsim"
TemG_physorUpper<-TemG_physor[,-which(colnames(TemG_physor) == "grids")]
TemG_physorUpper<-as.data.frame(TemG_physorUpper[upper.tri(TemG_physorUpper, diag = FALSE)])
TemG_physorUpper$biomeabbr<-"TemG"
colnames(TemG_physorUpper)[1]<-"physor"
TemG_physimUpper<-TemG_physim[,-which(colnames(TemG_physim) == "grids")]
TemG_physimUpper<-as.data.frame(TemG_physimUpper[upper.tri(TemG_physimUpper, diag = FALSE)])
TemG_physimUpper$biomeabbr<-"TemG"
colnames(TemG_physimUpper)[1]<-"physim"


TroCF_taxsorUpper<-TroCF_taxsor[,-which(colnames(TroCF_taxsor) == "grids")]
TroCF_taxsorUpper<-as.data.frame(TroCF_taxsorUpper[upper.tri(TroCF_taxsorUpper, diag = FALSE)])
TroCF_taxsorUpper$biomeabbr<-"TroCF"
colnames(TroCF_taxsorUpper)[1]<-"taxsor"
TroCF_taxsimUpper<-TroCF_taxsim[,-which(colnames(TroCF_taxsim) == "grids")]
TroCF_taxsimUpper<-as.data.frame(TroCF_taxsimUpper[upper.tri(TroCF_taxsimUpper, diag = FALSE)])
TroCF_taxsimUpper$biomeabbr<-"TroCF"
colnames(TroCF_taxsimUpper)[1]<-"taxsim"
TroCF_physorUpper<-TroCF_physor[,-which(colnames(TroCF_physor) == "grids")]
TroCF_physorUpper<-as.data.frame(TroCF_physorUpper[upper.tri(TroCF_physorUpper, diag = FALSE)])
TroCF_physorUpper$biomeabbr<-"TroCF"
colnames(TroCF_physorUpper)[1]<-"physor"
TroCF_physimUpper<-TroCF_physim[,-which(colnames(TroCF_physim) == "grids")]
TroCF_physimUpper<-as.data.frame(TroCF_physimUpper[upper.tri(TroCF_physimUpper, diag = FALSE)])
TroCF_physimUpper$biomeabbr<-"TroCF"
colnames(TroCF_physimUpper)[1]<-"physim"


TroG_taxsorUpper<-TroG_taxsor[,-which(colnames(TroG_taxsor) == "grids")]
TroG_taxsorUpper<-as.data.frame(TroG_taxsorUpper[upper.tri(TroG_taxsorUpper, diag = FALSE)])
TroG_taxsorUpper$biomeabbr<-"TroG"
colnames(TroG_taxsorUpper)[1]<-"taxsor"
TroG_taxsimUpper<-TroG_taxsim[,-which(colnames(TroG_taxsim) == "grids")]
TroG_taxsimUpper<-as.data.frame(TroG_taxsimUpper[upper.tri(TroG_taxsimUpper, diag = FALSE)])
TroG_taxsimUpper$biomeabbr<-"TroG"
colnames(TroG_taxsimUpper)[1]<-"taxsim"
TroG_physorUpper<-TroG_physor[,-which(colnames(TroG_physor) == "grids")]
TroG_physorUpper<-as.data.frame(TroG_physorUpper[upper.tri(TroG_physorUpper, diag = FALSE)])
TroG_physorUpper$biomeabbr<-"TroG"
colnames(TroG_physorUpper)[1]<-"physor"
TroG_physimUpper<-TroG_physim[,-which(colnames(TroG_physim) == "grids")]
TroG_physimUpper<-as.data.frame(TroG_physimUpper[upper.tri(TroG_physimUpper, diag = FALSE)])
TroG_physimUpper$biomeabbr<-"TroG"
colnames(TroG_physimUpper)[1]<-"physim"


Tun_taxsorUpper<-Tun_taxsor[,-which(colnames(Tun_taxsor) == "grids")]
Tun_taxsorUpper<-as.data.frame(Tun_taxsorUpper[upper.tri(Tun_taxsorUpper, diag = FALSE)])
Tun_taxsorUpper$biomeabbr<-"Tun"
colnames(Tun_taxsorUpper)[1]<-"taxsor"
Tun_taxsimUpper<-Tun_taxsim[,-which(colnames(Tun_taxsim) == "grids")]
Tun_taxsimUpper<-as.data.frame(Tun_taxsimUpper[upper.tri(Tun_taxsimUpper, diag = FALSE)])
Tun_taxsimUpper$biomeabbr<-"Tun"
colnames(Tun_taxsimUpper)[1]<-"taxsim"
Tun_physorUpper<-Tun_physor[,-which(colnames(Tun_physor) == "grids")]
Tun_physorUpper<-as.data.frame(Tun_physorUpper[upper.tri(Tun_physorUpper, diag = FALSE)])
Tun_physorUpper$biomeabbr<-"Tun"
colnames(Tun_physorUpper)[1]<-"physor"
Tun_physimUpper<-Tun_physim[,-which(colnames(Tun_physim) == "grids")]
Tun_physimUpper<-as.data.frame(Tun_physimUpper[upper.tri(Tun_physimUpper, diag = FALSE)])
Tun_physimUpper$biomeabbr<-"Tun"
colnames(Tun_physimUpper)[1]<-"physim"





# Bind the matrices into a single df for each biome
DBF_beta<-cbind(DBF_taxsorUpper, DBF_taxsimUpper, DBF_physorUpper, DBF_physimUpper)
DBF_beta<-DBF_beta[,-c(2,4,6)]
Des_beta<-cbind(Des_taxsorUpper, Des_taxsimUpper, Des_physorUpper, Des_physimUpper)
Des_beta<-Des_beta[,-c(2,4,6)]
FloG_beta<-cbind(FloG_taxsorUpper, FloG_taxsimUpper, FloG_physorUpper, FloG_physimUpper)
FloG_beta<-FloG_beta[,-c(2,4,6)]
Mang_beta<-cbind(Mang_taxsorUpper, Mang_taxsimUpper, Mang_physorUpper, Mang_physimUpper)
Mang_beta<-Mang_beta[,-c(2,4,6)]
MBF_beta<-cbind(MBF_taxsorUpper, MBF_taxsimUpper, MBF_physorUpper, MBF_physimUpper)
MBF_beta<-MBF_beta[,-c(2,4,6)]
Medit_beta<-cbind(Medit_taxsorUpper, Medit_taxsimUpper, Medit_physorUpper, Medit_physimUpper)
Medit_beta<-Medit_beta[,-c(2,4,6)]
MonG_beta<-cbind(MonG_taxsorUpper, MonG_taxsimUpper, MonG_physorUpper, MonG_physimUpper)
MonG_beta<-MonG_beta[,-c(2,4,6)]
Tai_beta<-cbind(Tai_taxsorUpper, Tai_taxsimUpper, Tai_physorUpper, Tai_physimUpper)
Tai_beta<-Tai_beta[,-c(2,4,6)]
TemBF_beta<-cbind(TemBF_taxsorUpper, TemBF_taxsimUpper, TemBF_physorUpper, TemBF_physimUpper)
TemBF_beta<-TemBF_beta[,-c(2,4,6)]
TemCF_beta<-cbind(TemCF_taxsorUpper, TemCF_taxsimUpper, TemCF_physorUpper, TemCF_physimUpper)
TemCF_beta<-TemCF_beta[,-c(2,4,6)]
TemG_beta<-cbind(TemG_taxsorUpper, TemG_taxsimUpper, TemG_physorUpper, TemG_physimUpper)
TemG_beta<-TemG_beta[,-c(2,4,6)]
TroCF_beta<-cbind(TroCF_taxsorUpper, TroCF_taxsimUpper, TroCF_physorUpper, TroCF_physimUpper)
TroCF_beta<-TroCF_beta[,-c(2,4,6)]
TroG_beta<-cbind(TroG_taxsorUpper, TroG_taxsimUpper, TroG_physorUpper, TroG_physimUpper)
TroG_beta<-TroG_beta[,-c(2,4,6)]
Tun_beta<-cbind(Tun_taxsorUpper, Tun_taxsimUpper, Tun_physorUpper, Tun_physimUpper)
Tun_beta<-Tun_beta[,-c(2,4,6)]

# Bind the matrices into a single df across all biomes
bd_summary<-rbind(DBF_beta, Des_beta, FloG_beta, Mang_beta, MBF_beta, Medit_beta,
                  MonG_beta, Tai_beta, TemBF_beta, TemCF_beta, TemG_beta, TroCF_beta,
                  TroG_beta, Tun_beta)
write.csv(bd_summary, "bd_summary.csv", row.names = FALSE)


# Make sure data classes are what we need for plotting
bd_summary$biomeabbr<-as.factor(bd_summary$biomeabbr)
bd_summary$taxsim<-as.numeric(bd_summary$taxsim)
bd_summary$taxsor<-as.numeric(bd_summary$taxsor)
bd_summary$physor<-as.numeric(bd_summary$physor)
bd_summary$physim<-as.numeric(bd_summary$physim)



# Plot Simpson's TBD
fit_taxsim<-kruskal.test(taxsim~biomeabbr, data = bd_summary)
print(fit_taxsim)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result_taxsim <- dunn.test(bd_summary$taxsim, bd_summary$biomeabbr,
                         method = "bonferroni")
print(dunn_result_taxsim)

# Assign labels for results that are significant
sigs<-dunn_result_taxsim$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result_taxsim$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)

# Plot
Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")

maxes<-bd_summary %>% group_by(biomeabbr) %>%
  summarise(max = max(taxsim, na.rm = TRUE))

quartz(w=3.375, h=2.5)

ggplot(bd_summary, aes(x = biomeabbr, y = taxsim, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(bd_summary$taxsim, na.rm = TRUE), 1.05 * max(bd_summary$taxsim, na.rm = TRUE))) +
  labs(y = "Taxonomic Simpson's", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.02*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)



# Plot Sorensen's TBD
fit_taxsor<-kruskal.test(taxsor~biomeabbr, data = bd_summary)
print(fit_taxsor)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result_taxsor <- dunn.test(bd_summary$taxsor, bd_summary$biomeabbr,
                         method = "bonferroni")
print(dunn_result_taxsor)

# Assign labels for results that are significant
sigs<-dunn_result_taxsor$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result_taxsor$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)

maxes<-bd_summary %>% group_by(biomeabbr) %>%
  summarise(max = max(taxsor, na.rm = TRUE))

# Plot
quartz(w=3.375, h=2.5)

ggplot(bd_summary, aes(x = biomeabbr, y = taxsor, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(bd_summary$taxsim), 1.05 * max(bd_summary$taxsim))) +
  labs(y = "Taxonomic Sørensen's", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.02*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)




# Plot Simpson's PBD
fit_physim<-kruskal.test(physim~biomeabbr, data = bd_summary)
print(fit_physim)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result_physim <- dunn.test(bd_summary$physim, bd_summary$biomeabbr,
                         method = "bonferroni")
print(dunn_result_physim)

# Assign labels for results that are significant
sigs<-dunn_result_physim$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result_physim$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)

maxes<-bd_summary %>% group_by(biomeabbr) %>%
  summarise(max = max(physim, na.rm = TRUE))

# Plot
quartz(w=3.375, h=2.5)

ggplot(bd_summary, aes(x = biomeabbr, y = physim, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(bd_summary$taxsim), 1.05 * max(bd_summary$taxsim))) +
  labs(y = "Phylogenetic Simpson's", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.02*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)


# Plot Sorensen's PBD
fit_physor<-kruskal.test(physor~biomeabbr, data = bd_summary)
print(fit_physor)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result_physor <- dunn.test(bd_summary$physor, bd_summary$biomeabbr,
                         method = "bonferroni")
print(dunn_result_physor)

# Assign labels for results that are significant
sigs<-dunn_result_physor$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result_physor$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)


maxes<-bd_summary %>% group_by(biomeabbr) %>%
  summarise(max = max(physor, na.rm = TRUE))

# Plot
quartz(w=3.375, h=2.5)

ggplot(bd_summary, aes(x = biomeabbr, y = physor, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(bd_summary$taxsim), 1.05 * max(bd_summary$taxsim))) +
  labs(y = "Phylogenetic Sørensen's", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.02*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)


## Format species matrices for GDM ----

# Again, format in chunks by biome 

# Make sure that our phylogeny and our matrix (including our coordinates) match
DBFPhy$tip.label[1:5]
colnames(DBFMatCoordz)[1:5]
df<-DBFMatCoordz[,c("lon", "lat", "grids")]
DBFPhy$tip.label<-gsub("_", ".", DBFPhy$tip.label)
DBFMatCoordz<-DBFMatCoordz[,which(colnames(DBFMatCoordz) %in% DBFPhy$tip.label)]
DBFPhy$tip.label<-DBFPhy$tip.label[which(DBFPhy$tip.label %in% colnames(DBFMatCoordz))]
length(DBFPhy$tip.label)
DBFMatCoordz<-cbind(DBFMatCoordz, df)
write.tree(DBFPhy, "DBF/DBFPhy.tre")
write.csv(DBFMatCoordz, "DBF/DBFMatCoords.csv")
remove(df)

# Set up data into data format 2 for input into formatsitepair() in gdm
DBF_region_location <- DBFMatCoordz[,c("grids", "lon", "lat")]
names(DBF_region_location) <- c("site","Long","Lat")
head(DBF_region_location)

# Only keep columns where the names match with DBFPhy$tip.label.
DBFMat1 <- DBFMatCoordz[,-(which(is.na(match(names(DBFMatCoordz), DBFPhy$tip.label))==T))]
rownames(DBFMat1) <- DBF_region_location[,1]
DBFMat1[1:5,1:5]
write.csv(DBFMat1, "DBF/DBFMat1.csv", row.names = FALSE)
write.csv(DBF_region_location, "DBF/DBF_region_location.csv", row.names=FALSE)

# Store intermediate step as a df
DBF_temp <- NULL
for(i in 1:length(DBFMat1)){
  print(i)
  # Find the rows (sites) where the value is 1 and create a new dataframe with the 
    # species name and corresponding site for these rows.
  DBF_temp <- rbind(DBF_temp,cbind(rep(names(DBFMat1)[i], length(which(DBFMat1[,i]==1))), 
                                   rownames(DBFMat1)[which(DBFMat1[,i]==1)]))
}

head(DBF_temp)
# Convert to a df
DBF_temp <- as.data.frame(DBF_temp)
names(DBF_temp) <- c("species","site")
head(DBF_temp)

DBF_temp2 <- merge(DBF_temp, DBF_region_location, by = "site")
head(DBF_temp2)

# Combine these dfs to create a lat, long, species list
DBF_sppTab<-DBF_temp2[, c("species", "site", "Lat", "Long")]
head(DBF_sppTab)
remove(DBF_temp2)
remove(DBF_temp)

write.csv(DBF_sppTab, "DBF/DBF_sppTab.csv", row.names = FALSE)



# Again, repeat this 14 times



DesMatCoordz<-read.csv("Des/DesMatCoords.csv")
DesPhy<-read.tree("Des/DesPhy.tree")

DesPhy$tip.label[1:5]
colnames(DesMatCoordz)[1:5]
df<-DesMatCoordz[,c("lon", "lat", "grids")]
DesPhy$tip.label<-gsub("_", ".", DesPhy$tip.label)
DesMatCoordz<-DesMatCoordz[,which(colnames(DesMatCoordz) %in% DesPhy$tip.label)]
DesPhy$tip.label<-DesPhy$tip.label[which(DesPhy$tip.label %in% colnames(DesMatCoordz))]
length(DesPhy$tip.label)
DesMatCoordz<-cbind(DesMatCoordz, df)
write.tree(DesPhy, "Des/DesPhy.tre")
write.csv(DesMatCoordz, "Des/DesMatCoords.csv")
remove(df)

Des_region_location <- DesMatCoordz[,c("grids", "lon", "lat")]
names(Des_region_location) <- c("site","Long","Lat")
head(Des_region_location)
DesMat1 <- DesMatCoordz[,-(which(is.na(match(names(DesMatCoordz), DesPhy$tip.label))==T))]
rownames(DesMat1) <- Des_region_location[,1]
DesMat1[1:5,1:5]
write.csv(DesMat1, "Des/DesMat1.csv", row.names = FALSE)
write.csv(Des_region_location, "Des/Des_region_location.csv", row.names=FALSE)

Des_temp <- NULL
for(i in 1:length(DesMat1)){
  print(i)
  Des_temp <- rbind(Des_temp,cbind(rep(names(DesMat1)[i], length(which(DesMat1[,i]==1))), 
                                   rownames(DesMat1)[which(DesMat1[,i]==1)]))
}

head(Des_temp)
Des_temp <- as.data.frame(Des_temp)
names(Des_temp) <- c("species","site")
head(Des_temp)

Des_temp2 <- merge(Des_temp, Des_region_location, by = "site")
head(Des_temp2)

Des_sppTab<-Des_temp2[, c("species", "site", "Lat", "Long")]
head(Des_sppTab)
remove(Des_temp2)
remove(Des_temp)

write.csv(Des_sppTab, "Des/Des_sppTab.csv", row.names = FALSE)






FloGMatCoordz<-read.csv("FloG/FloGMatCoords.csv")
FloGPhy<-read.tree("FloG/FloGPhy.tree")

FloGPhy$tip.label[1:5]
colnames(FloGMatCoordz)[1:5]
df<-FloGMatCoordz[,c("lon", "lat", "grids")]
FloGPhy$tip.label<-gsub("_", ".", FloGPhy$tip.label)
FloGMatCoordz<-FloGMatCoordz[,which(colnames(FloGMatCoordz) %in% FloGPhy$tip.label)]
FloGPhy$tip.label<-FloGPhy$tip.label[which(FloGPhy$tip.label %in% colnames(FloGMatCoordz))]
length(FloGPhy$tip.label)
FloGMatCoordz<-cbind(FloGMatCoordz, df)
write.tree(FloGPhy, "FloG/FloGPhy.tre")
write.csv(FloGMatCoordz, "FloG/FloGMatCoords.csv")
remove(df)


FloG_region_location <- FloGMatCoordz[,c("grids", "lon", "lat")]
names(FloG_region_location) <- c("site","Long","Lat")
head(FloG_region_location)
FloGMat1 <- FloGMatCoordz[,-(which(is.na(match(names(FloGMatCoordz), FloGPhy$tip.label))==T))]
rownames(FloGMat1) <- FloG_region_location[,1]
FloGMat1[1:5,1:5]
write.csv(FloGMat1, "FloG/FloGMat1.csv", row.names = FALSE)
write.csv(FloG_region_location, "FloG/FloG_region_location.csv", row.names=FALSE)

FloG_temp <- NULL
for(i in 1:length(FloGMat1)){
  print(i)
  FloG_temp <- rbind(FloG_temp,cbind(rep(names(FloGMat1)[i], length(which(FloGMat1[,i]==1))), 
                                   rownames(FloGMat1)[which(FloGMat1[,i]==1)]))
}

head(FloG_temp)
FloG_temp <- as.data.frame(FloG_temp)
names(FloG_temp) <- c("species","site")
head(FloG_temp)

FloG_temp2 <- merge(FloG_temp, FloG_region_location, by = "site")
head(FloG_temp2)

FloG_sppTab<-FloG_temp2[, c("species", "site", "Lat", "Long")]
head(FloG_sppTab)
remove(FloG_temp2)
remove(FloG_temp)

write.csv(FloG_sppTab, "FloG/FloG_sppTab.csv", row.names = FALSE)





MangMatCoordz<-read.csv("Mang/MangMatCoords.csv")
MangPhy<-read.tree("Mang/MangPhy.tree")

MangPhy$tip.label[1:5]
colnames(MangMatCoordz)[1:5]
df<-MangMatCoordz[,c("lon", "lat", "grids")]
MangPhy$tip.label<-gsub("_", ".", MangPhy$tip.label)
MangMatCoordz<-MangMatCoordz[,which(colnames(MangMatCoordz) %in% MangPhy$tip.label)]
MangPhy$tip.label<-MangPhy$tip.label[which(MangPhy$tip.label %in% colnames(MangMatCoordz))]
length(MangPhy$tip.label)
MangMatCoordz<-cbind(MangMatCoordz, df)
write.tree(MangPhy, "Mang/MangPhy.tre")
write.csv(MangMatCoordz, "Mang/MangMatCoords.csv")
remove(df)

Mang_region_location <- MangMatCoordz[,c("grids", "lon", "lat")]
names(Mang_region_location) <- c("site","Long","Lat")
head(Mang_region_location)
MangMat1 <- MangMatCoordz[,-(which(is.na(match(names(MangMatCoordz), MangPhy$tip.label))==T))]
rownames(MangMat1) <- Mang_region_location[,1]
MangMat1[1:5,1:5]
write.csv(MangMat1, "Mang/MangMat1.csv", row.names = FALSE)
write.csv(Mang_region_location, "Mang/Mang_region_location.csv", row.names=FALSE)

Mang_temp <- NULL
for(i in 1:length(MangMat1)){
  print(i)
  Mang_temp <- rbind(Mang_temp,cbind(rep(names(MangMat1)[i], length(which(MangMat1[,i]==1))), 
                                   rownames(MangMat1)[which(MangMat1[,i]==1)]))
}

head(Mang_temp)
Mang_temp <- as.data.frame(Mang_temp)
names(Mang_temp) <- c("species","site")
head(Mang_temp)

Mang_temp2 <- merge(Mang_temp, Mang_region_location, by = "site")
head(Mang_temp2)

Mang_sppTab<-Mang_temp2[, c("species", "site", "Lat", "Long")]
head(Mang_sppTab)
remove(Mang_temp2)
remove(Mang_temp)

write.csv(Mang_sppTab, "Mang/Mang_sppTab.csv", row.names = FALSE)






MBFMatCoordz<-read.csv("MBF/MBFMatCoords.csv")
MBFPhy<-read.tree("MBF/MBFPhy.tree")

MBFPhy$tip.label[1:5]
colnames(MBFMatCoordz)[1:5]
df<-MBFMatCoordz[,c("lon", "lat", "grids")]
MBFPhy$tip.label<-gsub("_", ".", MBFPhy$tip.label)
MBFMatCoordz<-MBFMatCoordz[,which(colnames(MBFMatCoordz) %in% MBFPhy$tip.label)]
MBFPhy$tip.label<-MBFPhy$tip.label[which(MBFPhy$tip.label %in% colnames(MBFMatCoordz))]
length(MBFPhy$tip.label)
MBFMatCoordz<-cbind(MBFMatCoordz, df)
write.tree(MBFPhy, "MBF/MBFPhy.tre")
write.csv(MBFMatCoordz, "MBF/MBFMatCoords.csv")
remove(df)

MBF_region_location <- MBFMatCoordz[,c("grids", "lon", "lat")]
names(MBF_region_location) <- c("site","Long","Lat")
head(MBF_region_location)
MBFMat1 <- MBFMatCoordz[,-(which(is.na(match(names(MBFMatCoordz), MBFPhy$tip.label))==T))]
rownames(MBFMat1) <- MBF_region_location[,1]
MBFMat1[1:5,1:5]
write.csv(MBFMat1, "MBF/MBFMat1.csv", row.names = FALSE)
write.csv(MBF_region_location, "MBF/MBF_region_location.csv", row.names=FALSE)

MBF_temp <- NULL
for(i in 1:length(MBFMat1)){
  print(i)
  MBF_temp <- rbind(MBF_temp,cbind(rep(names(MBFMat1)[i], length(which(MBFMat1[,i]==1))), 
                                   rownames(MBFMat1)[which(MBFMat1[,i]==1)]))
}

head(MBF_temp)
MBF_temp <- as.data.frame(MBF_temp)
names(MBF_temp) <- c("species","site")
head(MBF_temp)

MBF_temp2 <- merge(MBF_temp, MBF_region_location, by = "site")
head(MBF_temp2)

MBF_sppTab<-MBF_temp2[, c("species", "site", "Lat", "Long")]
head(MBF_sppTab)
remove(MBF_temp2)
remove(MBF_temp)

write.csv(MBF_sppTab, "MBF/MBF_sppTab.csv", row.names = FALSE)






MeditMatCoordz<-read.csv("Medit/MeditMatCoords.csv")
MeditPhy<-read.tree("Medit/MeditPhy.tree")

MeditPhy$tip.label[1:5]
colnames(MeditMatCoordz)[1:5]
df<-MeditMatCoordz[,c("lon", "lat", "grids")]
MeditPhy$tip.label<-gsub("_", ".", MeditPhy$tip.label)
MeditMatCoordz<-MeditMatCoordz[,which(colnames(MeditMatCoordz) %in% MeditPhy$tip.label)]
MeditPhy$tip.label<-MeditPhy$tip.label[which(MeditPhy$tip.label %in% colnames(MeditMatCoordz))]
length(MeditPhy$tip.label)
MeditMatCoordz<-cbind(MeditMatCoordz, df)
write.tree(MeditPhy, "Medit/MeditPhy.tre")
write.csv(MeditMatCoordz, "Medit/MeditMatCoords.csv")
remove(df)

Medit_region_location <- MeditMatCoordz[,c("grids", "lon", "lat")]
names(Medit_region_location) <- c("site","Long","Lat")
head(Medit_region_location)
MeditMat1 <- MeditMatCoordz[,-(which(is.na(match(names(MeditMatCoordz), MeditPhy$tip.label))==T))]
rownames(MeditMat1) <- Medit_region_location[,1]
MeditMat1[1:5,1:5]
write.csv(MeditMat1, "Medit/MeditMat1.csv", row.names = FALSE)
write.csv(Medit_region_location, "Medit/Medit_region_location.csv", row.names=FALSE)

Medit_temp <- NULL
for(i in 1:length(MeditMat1)){
  print(i)
  Medit_temp <- rbind(Medit_temp,cbind(rep(names(MeditMat1)[i], length(which(MeditMat1[,i]==1))), 
                                   rownames(MeditMat1)[which(MeditMat1[,i]==1)]))
}

head(Medit_temp)
Medit_temp <- as.data.frame(Medit_temp)
names(Medit_temp) <- c("species","site")
head(Medit_temp)

Medit_temp2 <- merge(Medit_temp, Medit_region_location, by = "site")
head(Medit_temp2)

Medit_sppTab<-Medit_temp2[, c("species", "site", "Lat", "Long")]
head(Medit_sppTab)
remove(Medit_temp2)
remove(Medit_temp)

write.csv(Medit_sppTab, "Medit/Medit_sppTab.csv", row.names = FALSE)






MonGMatCoordz<-read.csv("MonG/MonGMatCoords.csv")
MonGPhy<-read.tree("MonG/MonGPhy.tree")

MonGPhy$tip.label[1:5]
colnames(MonGMatCoordz)[1:5]
df<-MonGMatCoordz[,c("lon", "lat", "grids")]
MonGPhy$tip.label<-gsub("_", ".", MonGPhy$tip.label)
MonGMatCoordz<-MonGMatCoordz[,which(colnames(MonGMatCoordz) %in% MonGPhy$tip.label)]
MonGPhy$tip.label<-MonGPhy$tip.label[which(MonGPhy$tip.label %in% colnames(MonGMatCoordz))]
length(MonGPhy$tip.label)
MonGMatCoordz<-cbind(MonGMatCoordz, df)
write.tree(MonGPhy, "MonG/MonGPhy.tre")
write.csv(MonGMatCoordz, "MonG/MonGMatCoords.csv")
remove(df)

MonG_region_location <- MonGMatCoordz[,c("grids", "lon", "lat")]
names(MonG_region_location) <- c("site","Long","Lat")
head(MonG_region_location)
MonGMat1 <- MonGMatCoordz[,-(which(is.na(match(names(MonGMatCoordz), MonGPhy$tip.label))==T))]
rownames(MonGMat1) <- MonG_region_location[,1]
MonGMat1[1:5,1:5]
write.csv(MonGMat1, "MonG/MonGMat1.csv", row.names = FALSE)
write.csv(MonG_region_location, "MonG/MonG_region_location.csv", row.names=FALSE)

MonG_temp <- NULL
for(i in 1:length(MonGMat1)){
  print(i)
  MonG_temp <- rbind(MonG_temp,cbind(rep(names(MonGMat1)[i], length(which(MonGMat1[,i]==1))), 
                                   rownames(MonGMat1)[which(MonGMat1[,i]==1)]))
}

head(MonG_temp)
MonG_temp <- as.data.frame(MonG_temp)
names(MonG_temp) <- c("species","site")
head(MonG_temp)

MonG_temp2 <- merge(MonG_temp, MonG_region_location, by = "site")
head(MonG_temp2)

MonG_sppTab<-MonG_temp2[, c("species", "site", "Lat", "Long")]
head(MonG_sppTab)
remove(MonG_temp2)
remove(MonG_temp)

write.csv(MonG_sppTab, "MonG/MonG_sppTab.csv", row.names = FALSE)




TaiMatCoordz<-read.csv("Tai/TaiMatCoords.csv")
TaiPhy<-read.tree("Tai/TaiPhy.tree")

TaiPhy$tip.label[1:5]
colnames(TaiMatCoordz)[1:5]
df<-TaiMatCoordz[,c("lon", "lat", "grids")]
TaiPhy$tip.label<-gsub("_", ".", TaiPhy$tip.label)
TaiMatCoordz<-TaiMatCoordz[,which(colnames(TaiMatCoordz) %in% TaiPhy$tip.label)]
TaiPhy$tip.label<-TaiPhy$tip.label[which(TaiPhy$tip.label %in% colnames(TaiMatCoordz))]
length(TaiPhy$tip.label)
TaiMatCoordz<-cbind(TaiMatCoordz, df)
write.tree(TaiPhy, "Tai/TaiPhy.tre")
write.csv(TaiMatCoordz, "Tai/TaiMatCoords.csv")
remove(df)

Tai_region_location <- TaiMatCoordz[,c("grids", "lon", "lat")]
names(Tai_region_location) <- c("site","Long","Lat")
head(Tai_region_location)
TaiMat1 <- TaiMatCoordz[,-(which(is.na(match(names(TaiMatCoordz), TaiPhy$tip.label))==T))]
rownames(TaiMat1) <- Tai_region_location[,1]
TaiMat1[1:5,1:5]
write.csv(TaiMat1, "Tai/TaiMat1.csv", row.names = FALSE)
write.csv(Tai_region_location, "Tai/Tai_region_location.csv", row.names=FALSE)

Tai_temp <- NULL
for(i in 1:length(TaiMat1)){
  print(i)
  Tai_temp <- rbind(Tai_temp,cbind(rep(names(TaiMat1)[i], length(which(TaiMat1[,i]==1))), 
                                   rownames(TaiMat1)[which(TaiMat1[,i]==1)]))
}

head(Tai_temp)
Tai_temp <- as.data.frame(Tai_temp)
names(Tai_temp) <- c("species","site")
head(Tai_temp)

Tai_temp2 <- merge(Tai_temp, Tai_region_location, by = "site")
head(Tai_temp2)

Tai_sppTab<-Tai_temp2[, c("species", "site", "Lat", "Long")]
head(Tai_sppTab)
remove(Tai_temp2)
remove(Tai_temp)

write.csv(Tai_sppTab, "Tai/Tai_sppTab.csv", row.names = FALSE)





TemBFMatCoordz<-read.csv("TemBF/TemBFMatCoords.csv")
TemBFPhy<-read.tree("TemBF/TemBFPhy.tree")

TemBFPhy$tip.label[1:5]
colnames(TemBFMatCoordz)[1:5]
df<-TemBFMatCoordz[,c("lon", "lat", "grids")]
TemBFPhy$tip.label<-gsub("_", ".", TemBFPhy$tip.label)
TemBFMatCoordz<-TemBFMatCoordz[,which(colnames(TemBFMatCoordz) %in% TemBFPhy$tip.label)]
TemBFPhy$tip.label<-TemBFPhy$tip.label[which(TemBFPhy$tip.label %in% colnames(TemBFMatCoordz))]
length(TemBFPhy$tip.label)
TemBFMatCoordz<-cbind(TemBFMatCoordz, df)
write.tree(TemBFPhy, "TemBF/TemBFPhy.tre")
write.csv(TemBFMatCoordz, "TemBF/TemBFMatCoords.csv")
remove(df)

TemBF_region_location <- TemBFMatCoordz[,c("grids", "lon", "lat")]
names(TemBF_region_location) <- c("site","Long","Lat")
head(TemBF_region_location)
TemBFMat1 <- TemBFMatCoordz[,-(which(is.na(match(names(TemBFMatCoordz), TemBFPhy$tip.label))==T))]
rownames(TemBFMat1) <- TemBF_region_location[,1]
TemBFMat1[1:5,1:5]
write.csv(TemBFMat1, "TemBF/TemBFMat1.csv", row.names = FALSE)
write.csv(TemBF_region_location, "TemBF/TemBF_region_location.csv", row.names=FALSE)

TemBF_temp <- NULL
for(i in 1:length(TemBFMat1)){
  print(i)
  TemBF_temp <- rbind(TemBF_temp,cbind(rep(names(TemBFMat1)[i], length(which(TemBFMat1[,i]==1))), 
                                   rownames(TemBFMat1)[which(TemBFMat1[,i]==1)]))
}

head(TemBF_temp)
TemBF_temp <- as.data.frame(TemBF_temp)
names(TemBF_temp) <- c("species","site")
head(TemBF_temp)

TemBF_temp2 <- merge(TemBF_temp, TemBF_region_location, by = "site")
head(TemBF_temp2)

TemBF_sppTab<-TemBF_temp2[, c("species", "site", "Lat", "Long")]
head(TemBF_sppTab)
remove(TemBF_temp2)
remove(TemBF_temp)

write.csv(TemBF_sppTab, "TemBF/TemBF_sppTab.csv", row.names = FALSE)






TemCFMatCoordz<-read.csv("TemCF/TemCFMatCoords.csv")
TemCFPhy<-read.tree("TemCF/TemCFPhy.tree")

TemCFPhy$tip.label[1:5]
colnames(TemCFMatCoordz)[1:5]
df<-TemCFMatCoordz[,c("lon", "lat", "grids")]
TemCFPhy$tip.label<-gsub("_", ".", TemCFPhy$tip.label)
TemCFMatCoordz<-TemCFMatCoordz[,which(colnames(TemCFMatCoordz) %in% TemCFPhy$tip.label)]
TemCFPhy$tip.label<-TemCFPhy$tip.label[which(TemCFPhy$tip.label %in% colnames(TemCFMatCoordz))]
length(TemCFPhy$tip.label)
TemCFMatCoordz<-cbind(TemCFMatCoordz, df)
write.tree(TemCFPhy, "TemCF/TemCFPhy.tre")
write.csv(TemCFMatCoordz, "TemCF/TemCFMatCoords.csv")
remove(df)

TemCF_region_location <- TemCFMatCoordz[,c("grids", "lon", "lat")]
names(TemCF_region_location) <- c("site","Long","Lat")
head(TemCF_region_location)
TemCFMat1 <- TemCFMatCoordz[,-(which(is.na(match(names(TemCFMatCoordz), TemCFPhy$tip.label))==T))]
rownames(TemCFMat1) <- TemCF_region_location[,1]
TemCFMat1[1:5,1:5]
write.csv(TemCFMat1, "TemCF/TemCFMat1.csv", row.names = FALSE)
write.csv(TemCF_region_location, "TemCF/TemCF_region_location.csv", row.names=FALSE)

TemCF_temp <- NULL
for(i in 1:length(TemCFMat1)){
  print(i)
  TemCF_temp <- rbind(TemCF_temp,cbind(rep(names(TemCFMat1)[i], length(which(TemCFMat1[,i]==1))), 
                                   rownames(TemCFMat1)[which(TemCFMat1[,i]==1)]))
}

head(TemCF_temp)
TemCF_temp <- as.data.frame(TemCF_temp)
names(TemCF_temp) <- c("species","site")
head(TemCF_temp)

TemCF_temp2 <- merge(TemCF_temp, TemCF_region_location, by = "site")
head(TemCF_temp2)

TemCF_sppTab<-TemCF_temp2[, c("species", "site", "Lat", "Long")]
head(TemCF_sppTab)
remove(TemCF_temp2)
remove(TemCF_temp)

write.csv(TemCF_sppTab, "TemCF/TemCF_sppTab.csv", row.names = FALSE)






TemGMatCoordz<-read.csv("TemG/TemGMatCoords.csv")
TemGPhy<-read.tree("TemG/TemGPhy.tree")

TemGPhy$tip.label[1:5]
colnames(TemGMatCoordz)[1:5]
df<-TemGMatCoordz[,c("lon", "lat", "grids")]
TemGPhy$tip.label<-gsub("_", ".", TemGPhy$tip.label)
TemGMatCoordz<-TemGMatCoordz[,which(colnames(TemGMatCoordz) %in% TemGPhy$tip.label)]
TemGPhy$tip.label<-TemGPhy$tip.label[which(TemGPhy$tip.label %in% colnames(TemGMatCoordz))]
length(TemGPhy$tip.label)
TemGMatCoordz<-cbind(TemGMatCoordz, df)
write.tree(TemGPhy, "TemG/TemGPhy.tre")
write.csv(TemGMatCoordz, "TemG/TemGMatCoords.csv")
remove(df)

TemG_region_location <- TemGMatCoordz[,c("grids", "lon", "lat")]
names(TemG_region_location) <- c("site","Long","Lat")
head(TemG_region_location)
TemGMat1 <- TemGMatCoordz[,-(which(is.na(match(names(TemGMatCoordz), TemGPhy$tip.label))==T))]
rownames(TemGMat1) <- TemG_region_location[,1]
TemGMat1[1:5,1:5]
write.csv(TemGMat1, "TemG/TemGMat1.csv", row.names = FALSE)
write.csv(TemG_region_location, "TemG/TemG_region_location.csv", row.names=FALSE)

TemG_temp <- NULL
for(i in 1:length(TemGMat1)){
  print(i)
  TemG_temp <- rbind(TemG_temp,cbind(rep(names(TemGMat1)[i], length(which(TemGMat1[,i]==1))), 
                                   rownames(TemGMat1)[which(TemGMat1[,i]==1)]))
}

head(TemG_temp)
TemG_temp <- as.data.frame(TemG_temp)
names(TemG_temp) <- c("species","site")
head(TemG_temp)

TemG_temp2 <- merge(TemG_temp, TemG_region_location, by = "site")
head(TemG_temp2)

TemG_sppTab<-TemG_temp2[, c("species", "site", "Lat", "Long")]
head(TemG_sppTab)
remove(TemG_temp2)
remove(TemG_temp)

write.csv(TemG_sppTab, "TemG/TemG_sppTab.csv", row.names = FALSE)






TroCFMatCoordz<-read.csv("TroCF/TroCFMatCoords.csv")
TroCFPhy<-read.tree("TroCF/TroCFPhy.tree")

TroCFPhy$tip.label[1:5]
colnames(TroCFMatCoordz)[1:5]
df<-TroCFMatCoordz[,c("lon", "lat", "grids")]
TroCFPhy$tip.label<-gsub("_", ".", TroCFPhy$tip.label)
TroCFMatCoordz<-TroCFMatCoordz[,which(colnames(TroCFMatCoordz) %in% TroCFPhy$tip.label)]
TroCFPhy$tip.label<-TroCFPhy$tip.label[which(TroCFPhy$tip.label %in% colnames(TroCFMatCoordz))]
length(TroCFPhy$tip.label)
TroCFMatCoordz<-cbind(TroCFMatCoordz, df)
write.tree(TroCFPhy, "TroCF/TroCFPhy.tre")
write.csv(TroCFMatCoordz, "TroCF/TroCFMatCoords.csv")
remove(df)

TroCF_region_location <- TroCFMatCoordz[,c("grids", "lon", "lat")]
names(TroCF_region_location) <- c("site","Long","Lat")
head(TroCF_region_location)
TroCFMat1 <- TroCFMatCoordz[,-(which(is.na(match(names(TroCFMatCoordz), TroCFPhy$tip.label))==T))]
rownames(TroCFMat1) <- TroCF_region_location[,1]
TroCFMat1[1:5,1:5]
write.csv(TroCFMat1, "TroCF/TroCFMat1.csv", row.names = FALSE)
write.csv(TroCF_region_location, "TroCF/TroCF_region_location.csv", row.names=FALSE)

TroCF_temp <- NULL
for(i in 1:length(TroCFMat1)){
  print(i)
  TroCF_temp <- rbind(TroCF_temp,cbind(rep(names(TroCFMat1)[i], length(which(TroCFMat1[,i]==1))), 
                                   rownames(TroCFMat1)[which(TroCFMat1[,i]==1)]))
}

head(TroCF_temp)
TroCF_temp <- as.data.frame(TroCF_temp)
names(TroCF_temp) <- c("species","site")
head(TroCF_temp)

TroCF_temp2 <- merge(TroCF_temp, TroCF_region_location, by = "site")
head(TroCF_temp2)

TroCF_sppTab<-TroCF_temp2[, c("species", "site", "Lat", "Long")]
head(TroCF_sppTab)
remove(TroCF_temp2)
remove(TroCF_temp)

write.csv(TroCF_sppTab, "TroCF/TroCF_sppTab.csv", row.names = FALSE)






TroGMatCoordz<-read.csv("TroG/TroGMatCoords.csv")
TroGPhy<-read.tree("TroG/TroGPhy.tree")

TroGPhy$tip.label[1:5]
colnames(TroGMatCoordz)[1:5]
df<-TroGMatCoordz[,c("lon", "lat", "grids")]
TroGPhy$tip.label<-gsub("_", ".", TroGPhy$tip.label)
TroGMatCoordz<-TroGMatCoordz[,which(colnames(TroGMatCoordz) %in% TroGPhy$tip.label)]
TroGPhy$tip.label<-TroGPhy$tip.label[which(TroGPhy$tip.label %in% colnames(TroGMatCoordz))]
length(TroGPhy$tip.label)
TroGMatCoordz<-cbind(TroGMatCoordz, df)
write.tree(TroGPhy, "TroG/TroGPhy.tre")
write.csv(TroGMatCoordz, "TroG/TroGMatCoords.csv")
remove(df)

TroG_region_location <- TroGMatCoordz[,c("grids", "lon", "lat")]
names(TroG_region_location) <- c("site","Long","Lat")
head(TroG_region_location)
TroGMat1 <- TroGMatCoordz[,-(which(is.na(match(names(TroGMatCoordz), TroGPhy$tip.label))==T))]
rownames(TroGMat1) <- TroG_region_location[,1]
TroGMat1[1:5,1:5]
write.csv(TroGMat1, "TroG/TroGMat1.csv", row.names = FALSE)
write.csv(TroG_region_location, "TroG/TroG_region_location.csv", row.names=FALSE)

TroG_temp <- NULL
for(i in 1:length(TroGMat1)){
  print(i)
  TroG_temp <- rbind(TroG_temp,cbind(rep(names(TroGMat1)[i], length(which(TroGMat1[,i]==1))), 
                                   rownames(TroGMat1)[which(TroGMat1[,i]==1)]))
}

head(TroG_temp)
TroG_temp <- as.data.frame(TroG_temp)
names(TroG_temp) <- c("species","site")
head(TroG_temp)

TroG_temp2 <- merge(TroG_temp, TroG_region_location, by = "site")
head(TroG_temp2)

TroG_sppTab<-TroG_temp2[, c("species", "site", "Lat", "Long")]
head(TroG_sppTab)
remove(TroG_temp2)
remove(TroG_temp)

write.csv(TroG_sppTab, "TroG/TroG_sppTab.csv", row.names = FALSE)






TunMatCoordz<-read.csv("Tun/TunMatCoords.csv")
TunPhy<-read.tree("Tun/TunPhy.tree")

TunPhy$tip.label[1:5]
colnames(TunMatCoordz)[1:5]
df<-TunMatCoordz[,c("lon", "lat", "grids")]
TunPhy$tip.label<-gsub("_", ".", TunPhy$tip.label)
TunMatCoordz<-TunMatCoordz[,which(colnames(TunMatCoordz) %in% TunPhy$tip.label)]
TunPhy$tip.label<-TunPhy$tip.label[which(TunPhy$tip.label %in% colnames(TunMatCoordz))]
length(TunPhy$tip.label)
TunMatCoordz<-cbind(TunMatCoordz, df)
write.tree(TunPhy, "Tun/TunPhy.tre")
write.csv(TunMatCoordz, "Tun/TunMatCoords.csv")
remove(df)

Tun_region_location <- TunMatCoordz[,c("grids", "lon", "lat")]
names(Tun_region_location) <- c("site","Long","Lat")
head(Tun_region_location)
TunMat1 <- TunMatCoordz[,-(which(is.na(match(names(TunMatCoordz), TunPhy$tip.label))==T))]
rownames(TunMat1) <- Tun_region_location[,1]
TunMat1[1:5,1:5]
write.csv(TunMat1, "Tun/TunMat1.csv", row.names = FALSE)
write.csv(Tun_region_location, "Tun/Tun_region_location.csv", row.names=FALSE)

Tun_temp <- NULL
for(i in 1:length(TunMat1)){
  print(i)
  Tun_temp <- rbind(Tun_temp,cbind(rep(names(TunMat1)[i], length(which(TunMat1[,i]==1))), 
                                   rownames(TunMat1)[which(TunMat1[,i]==1)]))
}

head(Tun_temp)
Tun_temp <- as.data.frame(Tun_temp)
names(Tun_temp) <- c("species","site")
head(Tun_temp)

Tun_temp2 <- merge(Tun_temp, Tun_region_location, by = "site")
head(Tun_temp2)

Tun_sppTab<-Tun_temp2[, c("species", "site", "Lat", "Long")]
head(Tun_sppTab)
remove(Tun_temp2)
remove(Tun_temp)

write.csv(Tun_sppTab, "Tun/Tun_sppTab.csv", row.names = FALSE)






AllMatCoordz<-read.csv("All/AllMatCoords.csv")
AllPhy<-read.tree("All/AllPhy.tre")

AllPhy$tip.label[1:5]
colnames(AllMatCoordz)[1:5]
df<-AllMatCoordz[,c("lon", "lat", "grids")]
AllPhy$tip.label<-gsub("_", ".", AllPhy$tip.label)
AllMatCoordz<-AllMatCoordz[,which(colnames(AllMatCoordz) %in% AllPhy$tip.label)]
AllPhy$tip.label<-AllPhy$tip.label[which(AllPhy$tip.label %in% colnames(AllMatCoordz))]
length(AllPhy$tip.label)
AllMatCoordz<-cbind(AllMatCoordz, df)
write.tree(AllPhy, "All/AllPhy.tre")
write.csv(AllMatCoordz, "All/AllMatCoords.csv")
remove(df)


All_region_location <- AllMatCoordz[,c("grids", "lon", "lat")]
names(All_region_location) <- c("site","Long","Lat")
head(All_region_location)
AllMat1 <- AllMatCoordz[,-c(68516,68517,68518)]
rownames(AllMat1) <- All_region_location$site
AllMat1[1:5,1:5]
write.csv(AllMat1, "All/AllMat1.csv", row.names = FALSE)
write.csv(All_region_location, "All/All_region_location.csv", row.names=FALSE)

All_temp <- NULL
for(i in 1:length(AllMat1)){
  print(i)
  All_temp <- rbind(All_temp,cbind(rep(names(AllMat1)[i], length(which(AllMat1[,i]==1))), 
                                   rownames(AllMat1)[which(AllMat1[,i]==1)]))
}

head(All_temp)
All_temp <- as.data.frame(All_temp)
names(All_temp) <- c("species","site")
head(All_temp)

All_temp2 <- merge(All_temp, All_region_location, by = "site")
head(All_temp2)

All_sppTab<-All_temp2[, c("species", "site", "Lat", "Long")]
head(All_sppTab)
remove(All_temp2)
remove(All_temp)

write.csv(All_sppTab, "All/All_sppTab.csv", row.names = FALSE)




## GDMs ----

# DBF GDM
DBFMat1<-DBFMat1[-270,] # this site is missing from our spp_Tab, gonna remove it
DBF_sppTab<-DBF_sppTab[which(DBF_sppTab$species %in% colnames(DBFMat1)),]
DBFMat1<-DBFMat1[,which(colnames(DBFMat1) %in% DBF_sppTab$species), drop = FALSE]
# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = DBFPhy
predData <- envRast
rasts <- envRast

#DBF phylogenetic gdm
DBF_phylogdmTab <- formatsitepair_phylo(DBF_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        # The following arguments are added in to 
                                        # run phylobeta in phyloregion
                                        community_data = DBFMat1,
                                        thetree = DBFPhy)

DBF_phylogdmTab <- na.omit(DBF_phylogdmTab) # went from 46056 to 44850
write.csv(DBF_phylogdmTab, "DBF/DBF_phylogdmTab.csv", row.names = FALSE)
DBF_gdm.phylo <- gdm(DBF_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))

# COPY THE SUMMARY AND PASTE in a text file
# Then, data were manually input into a csv file of biomes by raw predictor values
# These were scaled into proportions later, when plotting
summary(DBF_gdm.phylo) # This shows how much deviance is explained and the 
# contribution of each predictor variable


# DBF Taxonomic GDM
DBF_taxgdmTab<-formatsitepair_tax(DBF_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                               predData=envRast, siteColumn="site",sppColumn = "species")

DBF_taxgdmTab <- na.omit(DBF_taxgdmTab)
DBF_gdm.tax <- gdm(DBF_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))

# COPY THE SUMMARY AND PASTE in a text file
# Then, data were manually input into a csv file of biomes by raw predictor values
# These were scaled into proportions later, when plotting
summary(DBF_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(DBF_taxgdmTab, "DBF/DBF_taxgdmTab.csv", row.names=FALSE)






# Des GDM
Des_sppTab<-as.data.frame(fread("Des/Des_sppTab.csv"))
DesMat1<-as.data.frame(fread("Des/DesMat1.csv"))
DesPhy<-read.tree("Des/DesPhy.tre")
Des_sppTab<-Des_sppTab[which(Des_sppTab$species %in% colnames(DesMat1)),]
DesMat1<-DesMat1[,which(colnames(DesMat1) %in% Des_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = DesPhy
predData <- envRast
rasts <- envRast

#Des phylogenetic gdm
Des_phylogdmTab <- formatsitepair_phylo(Des_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = DesMat1,
                                        thetree = DesPhy)

Des_phylogdmTab <- na.omit(Des_phylogdmTab)
write.csv(Des_phylogdmTab, "Des/Des_phylogdmTab.csv", row.names = FALSE)
Des_gdm.phylo <- gdm(Des_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))

# COPY THE SUMMARY AND PASTE in a text file
# Then, data were manually input into a csv file of biomes by raw predictor values
# These were scaled into proportions later, when plotting
summary(Des_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# Des Taxonomic GDM
Des_taxgdmTab<-formatsitepair_tax(Des_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

Des_taxgdmTab <- na.omit(Des_taxgdmTab)
Des_gdm.tax <- gdm(Des_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))

# COPY THE SUMMARY AND PASTE in a text file
# Then, data were manually input into a csv file of biomes by raw predictor values
# These were scaled into proportions later, when plotting
summary(Des_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(Des_taxgdmTab, "Des/Des_taxgdmTab.csv", row.names=FALSE)






# FloG GDM
FloG_sppTab<-as.data.frame(fread("FloG/FloG_sppTab.csv"))
FloGMat1<-as.data.frame(fread("FloG/FloGMat1.csv"))
FloGPhy<-read.tree("FloG/FloGPhy.tre")
FloG_sppTab<-FloG_sppTab[which(FloG_sppTab$species %in% colnames(FloGMat1)),]
FloGMat1<-FloGMat1[,which(colnames(FloGMat1) %in% FloG_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = FloGPhy
predData <- envRast
rasts <- envRast

#FloG phylogenetic gdm
FloG_phylogdmTab <- formatsitepair_phylo(FloG_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = FloGMat1,
                                        thetree = FloGPhy)

FloG_phylogdmTab <- na.omit(FloG_phylogdmTab)
write.csv(FloG_phylogdmTab, "FloG/FloG_phylogdmTab.csv", row.names = FALSE)
FloG_gdm.phylo <- gdm(FloG_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(FloG_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# FloG Taxonomic GDM
FloG_taxgdmTab<-formatsitepair_tax(FloG_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

FloG_taxgdmTab <- na.omit(FloG_taxgdmTab)
FloG_gdm.tax <- gdm(FloG_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(FloG_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(FloG_taxgdmTab, "FloG/FloG_taxgdmTab.csv", row.names=FALSE)






# Mang GDM
Mang_sppTab<-as.data.frame(fread("Mang/Mang_sppTab.csv"))
MangMat1<-as.data.frame(fread("Mang/MangMat1.csv"))
MangPhy<-read.tree("Mang/MangPhy.tre")
Mang_sppTab<-Mang_sppTab[which(Mang_sppTab$species %in% colnames(MangMat1)),]
MangMat1<-MangMat1[,which(colnames(MangMat1) %in% Mang_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = MangPhy
predData <- envRast
rasts <- envRast

#Mang phylogenetic gdm
Mang_phylogdmTab <- formatsitepair_phylo(Mang_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = MangMat1,
                                        thetree = MangPhy)

Mang_phylogdmTab <- na.omit(Mang_phylogdmTab)
write.csv(Mang_phylogdmTab, "Mang/Mang_phylogdmTab.csv", row.names = FALSE)
Mang_gdm.phylo <- gdm(Mang_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(Mang_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# Mang Taxonomic GDM
Mang_taxgdmTab<-formatsitepair_tax(Mang_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

Mang_taxgdmTab <- na.omit(Mang_taxgdmTab)
Mang_gdm.tax <- gdm(Mang_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(Mang_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(Mang_taxgdmTab, "Mang/Mang_taxgdmTab.csv", row.names=FALSE)






# MBF GDM
MBF_sppTab<-as.data.frame(fread("MBF/MBF_sppTab.csv"))
MBFMat1<-as.data.frame(fread("MBF/MBFMat1.csv"))
MBFPhy<-read.tree("MBF/MBFPhy.tre")
MBF_sppTab<-MBF_sppTab[which(MBF_sppTab$species %in% colnames(MBFMat1)),]
MBFMat1<-MBFMat1[,which(colnames(MBFMat1) %in% MBF_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = MBFPhy
predData <- envRast
rasts <- envRast

#MBF phylogenetic gdm
MBF_phylogdmTab <- formatsitepair_phylo(MBF_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = MBFMat1,
                                        thetree = MBFPhy)

MBF_phylogdmTab <- na.omit(MBF_phylogdmTab)
write.csv(MBF_phylogdmTab, "MBF/MBF_phylogdmTab.csv", row.names = FALSE)
MBF_gdm.phylo <- gdm(MBF_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(MBF_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# MBF Taxonomic GDM
MBF_taxgdmTab<-formatsitepair_tax(MBF_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

MBF_taxgdmTab <- na.omit(MBF_taxgdmTab)
MBF_gdm.tax <- gdm(MBF_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(MBF_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(MBF_taxgdmTab, "MBF/MBF_taxgdmTab.csv", row.names=FALSE)






# Medit GDM
Medit_sppTab<-as.data.frame(fread("Medit/Medit_sppTab.csv"))
MeditMat1<-as.data.frame(fread("Medit/MeditMat1.csv"))
MeditPhy<-read.tree("Medit/MeditPhy.tre")
Medit_sppTab<-Medit_sppTab[which(Medit_sppTab$species %in% colnames(MeditMat1)),]
MeditMat1<-MeditMat1[,which(colnames(MeditMat1) %in% Medit_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = MeditPhy
predData <- envRast
rasts <- envRast

#Medit phylogenetic gdm
Medit_phylogdmTab <- formatsitepair_phylo(Medit_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = MeditMat1,
                                        thetree = MeditPhy)

Medit_phylogdmTab <- na.omit(Medit_phylogdmTab)
write.csv(Medit_phylogdmTab, "Medit/Medit_phylogdmTab.csv", row.names = FALSE)
Medit_gdm.phylo <- gdm(Medit_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(Medit_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# Medit Taxonomic GDM
Medit_taxgdmTab<-formatsitepair_tax(Medit_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

Medit_taxgdmTab <- na.omit(Medit_taxgdmTab)
Medit_gdm.tax <- gdm(Medit_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(Medit_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(Medit_taxgdmTab, "Medit/Medit_taxgdmTab.csv", row.names=FALSE)






# MonG GDM
MonG_sppTab<-as.data.frame(fread("MonG/MonG_sppTab.csv"))
MonGMat1<-as.data.frame(fread("MonG/MonGMat1.csv"))
MonGPhy<-read.tree("MonG/MonGPhy.tre")
MonG_sppTab<-MonG_sppTab[which(MonG_sppTab$species %in% colnames(MonGMat1)),]
MonGMat1<-MonGMat1[,which(colnames(MonGMat1) %in% MonG_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = MonGPhy
predData <- envRast
rasts <- envRast

#MonG phylogenetic gdm
MonG_phylogdmTab <- formatsitepair_phylo(MonG_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = MonGMat1,
                                        thetree = MonGPhy)

MonG_phylogdmTab <- na.omit(MonG_phylogdmTab)
write.csv(MonG_phylogdmTab, "MonG/MonG_phylogdmTab.csv", row.names = FALSE)
MonG_gdm.phylo <- gdm(MonG_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(MonG_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# MonG Taxonomic GDM
MonG_taxgdmTab<-formatsitepair_tax(MonG_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

MonG_taxgdmTab <- na.omit(MonG_taxgdmTab)
MonG_gdm.tax <- gdm(MonG_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(MonG_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(MonG_taxgdmTab, "MonG/MonG_taxgdmTab.csv", row.names=FALSE)






# Tai GDM
Tai_sppTab<-as.data.frame(fread("Tai/Tai_sppTab.csv"))
TaiMat1<-as.data.frame(fread("Tai/TaiMat1.csv"))
TaiPhy<-read.tree("Tai/TaiPhy.tre")
Tai_sppTab<-Tai_sppTab[which(Tai_sppTab$species %in% colnames(TaiMat1)),]
TaiMat1<-TaiMat1[,which(colnames(TaiMat1) %in% Tai_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TaiPhy
predData <- envRast
rasts <- envRast

#Tai phylogenetic gdm
Tai_phylogdmTab <- formatsitepair_phylo(Tai_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TaiMat1,
                                        thetree = TaiPhy)

Tai_phylogdmTab <- na.omit(Tai_phylogdmTab)
write.csv(Tai_phylogdmTab, "Tai/Tai_phylogdmTab.csv", row.names = FALSE)
Tai_gdm.phylo <- gdm(Tai_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(Tai_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# Tai Taxonomic GDM
Tai_taxgdmTab<-formatsitepair_tax(Tai_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

Tai_taxgdmTab <- na.omit(Tai_taxgdmTab)
Tai_gdm.tax <- gdm(Tai_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(Tai_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(Tai_taxgdmTab, "Tai/Tai_taxgdmTab.csv", row.names=FALSE)






# TemBF GDM
TemBF_sppTab<-as.data.frame(fread("TemBF/TemBF_sppTab.csv"))
TemBFMat1<-as.data.frame(fread("TemBF/TemBFMat1.csv"))
TemBFPhy<-read.tree("TemBF/TemBFPhy.tre")
TemBF_sppTab<-TemBF_sppTab[which(TemBF_sppTab$species %in% colnames(TemBFMat1)),]
TemBFMat1<-TemBFMat1[,which(colnames(TemBFMat1) %in% TemBF_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TemBFPhy
predData <- envRast
rasts <- envRast

#TemBF phylogenetic gdm
TemBF_phylogdmTab <- formatsitepair_phylo(TemBF_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TemBFMat1,
                                        thetree = TemBFPhy)

TemBF_phylogdmTab <- na.omit(TemBF_phylogdmTab)
write.csv(TemBF_phylogdmTab, "TemBF/TemBF_phylogdmTab.csv", row.names = FALSE)
TemBF_gdm.phylo <- gdm(TemBF_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(TemBF_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# TemBF Taxonomic GDM
TemBF_taxgdmTab<-formatsitepair_tax(TemBF_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

TemBF_taxgdmTab <- na.omit(TemBF_taxgdmTab)
TemBF_gdm.tax <- gdm(TemBF_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(TemBF_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(TemBF_taxgdmTab, "TemBF/TemBF_taxgdmTab.csv", row.names=FALSE)






# TemCF GDM
TemCF_sppTab<-as.data.frame(fread("TemCF/TemCF_sppTab.csv"))
TemCFMat1<-as.data.frame(fread("TemCF/TemCFMat1.csv"))
TemCFPhy<-read.tree("TemCF/TemCFPhy.tre")
TemCF_sppTab<-TemCF_sppTab[which(TemCF_sppTab$species %in% colnames(TemCFMat1)),]
TemCFMat1<-TemCFMat1[,which(colnames(TemCFMat1) %in% TemCF_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TemCFPhy
predData <- envRast
rasts <- envRast

#TemCF phylogenetic gdm
TemCF_phylogdmTab <- formatsitepair_phylo(TemCF_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TemCFMat1,
                                        thetree = TemCFPhy)

TemCF_phylogdmTab <- na.omit(TemCF_phylogdmTab)
write.csv(TemCF_phylogdmTab, "TemCF/TemCF_phylogdmTab.csv", row.names = FALSE)
TemCF_gdm.phylo <- gdm(TemCF_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(TemCF_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# TemCF Taxonomic GDM
TemCF_taxgdmTab<-formatsitepair_tax(TemCF_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

TemCF_taxgdmTab <- na.omit(TemCF_taxgdmTab)
TemCF_gdm.tax <- gdm(TemCF_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(TemCF_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(TemCF_taxgdmTab, "TemCF/TemCF_taxgdmTab.csv", row.names=FALSE)






# TemG GDM
TemG_sppTab<-as.data.frame(fread("TemG/TemG_sppTab.csv"))
TemGMat1<-as.data.frame(fread("TemG/TemGMat1.csv"))
TemGPhy<-read.tree("TemG/TemGPhy.tre")
TemG_sppTab<-TemG_sppTab[which(TemG_sppTab$species %in% colnames(TemGMat1)),]
TemGMat1<-TemGMat1[,which(colnames(TemGMat1) %in% TemG_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TemGPhy
predData <- envRast
rasts <- envRast

#TemG phylogenetic gdm
TemG_phylogdmTab <- formatsitepair_phylo(TemG_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TemGMat1,
                                        thetree = TemGPhy)

TemG_phylogdmTab <- na.omit(TemG_phylogdmTab)
write.csv(TemG_phylogdmTab, "TemG/TemG_phylogdmTab.csv", row.names = FALSE)
TemG_gdm.phylo <- gdm(TemG_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(TemG_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# TemG Taxonomic GDM
TemG_taxgdmTab<-formatsitepair_tax(TemG_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

TemG_taxgdmTab <- na.omit(TemG_taxgdmTab)
TemG_gdm.tax <- gdm(TemG_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(TemG_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(TemG_taxgdmTab, "TemG/TemG_taxgdmTab.csv", row.names=FALSE)






# TroCF GDM
TroCF_sppTab<-as.data.frame(fread("TroCF/TroCF_sppTab.csv"))
TroCFMat1<-as.data.frame(fread("TroCF/TroCFMat1.csv"))
TroCFPhy<-read.tree("TroCF/TroCFPhy.tre")
TroCF_sppTab<-TroCF_sppTab[which(TroCF_sppTab$species %in% colnames(TroCFMat1)),]
TroCFMat1<-TroCFMat1[,which(colnames(TroCFMat1) %in% TroCF_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TroCFPhy
predData <- envRast
rasts <- envRast

#TroCF phylogenetic gdm
TroCF_phylogdmTab <- formatsitepair_phylo(TroCF_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TroCFMat1,
                                        thetree = TroCFPhy)

TroCF_phylogdmTab <- na.omit(TroCF_phylogdmTab)
write.csv(TroCF_phylogdmTab, "TroCF/TroCF_phylogdmTab.csv", row.names = FALSE)
TroCF_gdm.phylo <- gdm(TroCF_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(TroCF_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# TroCF Taxonomic GDM
TroCF_taxgdmTab<-formatsitepair_tax(TroCF_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

TroCF_taxgdmTab <- na.omit(TroCF_taxgdmTab)
TroCF_gdm.tax <- gdm(TroCF_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(TroCF_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(TroCF_taxgdmTab, "TroCF/TroCF_taxgdmTab.csv", row.names=FALSE)






# TroG GDM
TroG_sppTab<-as.data.frame(fread("TroG/TroG_sppTab.csv"))
TroGMat1<-as.data.frame(fread("TroG/TroGMat1.csv"))
TroGPhy<-read.tree("TroG/TroGPhy.tre")
TroG_sppTab<-TroG_sppTab[which(TroG_sppTab$species %in% colnames(TroGMat1)),]
TroGMat1<-TroGMat1[,which(colnames(TroGMat1) %in% TroG_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TroGPhy
predData <- envRast
rasts <- envRast

#TroG phylogenetic gdm
TroG_phylogdmTab <- formatsitepair_phylo(TroG_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TroGMat1,
                                        thetree = TroGPhy)

TroG_phylogdmTab <- na.omit(TroG_phylogdmTab)
write.csv(TroG_phylogdmTab, "TroG/TroG_phylogdmTab.csv", row.names = FALSE)
TroG_gdm.phylo <- gdm(TroG_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(TroG_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# TroG Taxonomic GDM
TroG_taxgdmTab<-formatsitepair_tax(TroG_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

TroG_taxgdmTab <- na.omit(TroG_taxgdmTab)
TroG_gdm.tax <- gdm(TroG_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(TroG_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(TroG_taxgdmTab, "TroG/TroG_taxgdmTab.csv", row.names=FALSE)






# Tun GDM
Tun_sppTab<-as.data.frame(fread("Tun/Tun_sppTab.csv"))
TunMat1<-as.data.frame(fread("Tun/TunMat1.csv"))
TunPhy<-read.tree("Tun/TunPhy.tre")
Tun_sppTab<-Tun_sppTab[which(Tun_sppTab$species %in% colnames(TunMat1)),]
TunMat1<-TunMat1[,which(colnames(TunMat1) %in% Tun_sppTab$species), drop = FALSE]
TunMat1<-fread("Tun/TunMat1.csv")

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = TunPhy
predData <- envRast
rasts <- envRast

#Tun phylogenetic gdm
Tun_phylogdmTab <- formatsitepair_phylo(Tun_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = TunMat1,
                                        thetree = TunPhy)

Tun_phylogdmTab <- na.omit(Tun_phylogdmTab)
write.csv(Tun_phylogdmTab, "Tun/Tun_phylogdmTab.csv", row.names = FALSE)
Tun_gdm.phylo <- gdm(Tun_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(Tun_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# Tun Taxonomic GDM
Tun_taxgdmTab<-formatsitepair_tax(Tun_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

Tun_taxgdmTab <- na.omit(Tun_taxgdmTab)
Tun_gdm.tax <- gdm(Tun_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(Tun_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(Tun_taxgdmTab, "Tun/Tun_taxgdmTab.csv", row.names=FALSE)






# All GDM
All_sppTab<-as.data.frame(fread("All/All_sppTab.csv"))
AllMat1<-as.data.frame(fread("All/AllMat1.csv"))
AllPhy<-read.tree("All/AllPhy.tre")
All_sppTab<-All_sppTab[which(All_sppTab$species %in% colnames(AllMat1)),]
AllMat1<-AllMat1[,which(colnames(AllMat1) %in% All_sppTab$species), drop = FALSE]

# inputs for gdm
locs <- c("Long","Lat")
XColumn <- "Long"
YColumn <- "Lat"
siteColumn="site"
sppColumn = "species"
thetree = AllPhy
predData <- envRast
rasts <- envRast

#All phylogenetic gdm
All_phylogdmTab <- formatsitepair_phylo(All_sppTab, bioFormat = 2, XColumn="Long", 
                                        YColumn="Lat",predData=envRast, 
                                        siteColumn="site",sppColumn = "species",
                                        community_data = AllMat1,
                                        thetree = AllPhy)

All_phylogdmTab <- na.omit(All_phylogdmTab)
write.csv(All_phylogdmTab, "All/All_phylogdmTab.csv", row.names = FALSE)
All_gdm.phylo <- gdm(All_phylogdmTab, geo=T) 
#plot(gdm.phylo, plot.layout=c(4,3))
summary(All_gdm.phylo) #This shows how much deviance is explained and the contribution of each predictor variable


# All Taxonomic GDM
All_taxgdmTab<-formatsitepair_tax(All_sppTab, bioFormat = 2, XColumn="Long", YColumn="Lat", 
                              predData=envRast, siteColumn="site",sppColumn = "species")

All_taxgdmTab <- na.omit(All_taxgdmTab)
All_gdm.tax <- gdm(All_taxgdmTab, geo=T)
#plot(gdm.tax, plot.layout=c(4,3))
summary(All_gdm.tax) #This shows how much deviance is explained and the contribution of each predictor variable
write.csv(All_taxgdmTab, "All/All_taxgdmTab.csv", row.names=FALSE)





## Get lat/long of cell to match with gdm coordinates for regressions ----

# This code identifies the grid cell centroids that GDM assigns, so that we can 
  # match them back with the original grid cell centroid post GDM

DBF_cellNum <- as.data.frame(cellFromXY(envRast, DBF_region_location[-1]))
colnames(DBF_cellNum)<-"cellName"
head(DBF_cellNum)

# This is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location
DBF_cell_LatLong <- as.data.frame(xyFromCell(envRast, DBF_cellNum$cellName))
DBFLocs<-cbind(DBF_region_location, DBF_cell_LatLong)
colnames(DBFLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(DBFLocs)
remove(DBF_cellNum)
remove(DBF_cell_LatLong)
uniquez<-DBFLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(DBFLocs, "DBF/DBFLocs.csv", row.names = FALSE)


# Again, do this 14 more times.


Des_region_location<-read.csv("Des/Des_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
Des_cellNum <- as.data.frame(cellFromXY(envRast, Des_region_location[-1]))
colnames(Des_cellNum)<-"cellName"
head(Des_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

Des_cell_LatLong <- as.data.frame(xyFromCell(envRast, Des_cellNum$cellName))
DesLocs<-cbind(Des_region_location, Des_cell_LatLong)
colnames(DesLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(DesLocs)
remove(Des_cellNum)
remove(Des_cell_LatLong)
uniquez<-DesLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(DesLocs, "Des/DesLocs.csv", row.names = FALSE)





FloG_region_location<-read.csv("FloG/FloG_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
FloG_cellNum <- as.data.frame(cellFromXY(envRast, FloG_region_location[-1]))
colnames(FloG_cellNum)<-"cellName"
head(FloG_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

FloG_cell_LatLong <- as.data.frame(xyFromCell(envRast, FloG_cellNum$cellName))
FloGLocs<-cbind(FloG_region_location, FloG_cell_LatLong)
colnames(FloGLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(FloGLocs)
remove(FloG_cellNum)
remove(FloG_cell_LatLong)
uniquez<-FloGLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(FloGLocs, "FloG/FloGLocs.csv", row.names = FALSE)





Mang_region_location<-read.csv("Mang/Mang_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
Mang_cellNum <- as.data.frame(cellFromXY(envRast, Mang_region_location[-1]))
colnames(Mang_cellNum)<-"cellName"
head(Mang_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

Mang_cell_LatLong <- as.data.frame(xyFromCell(envRast, Mang_cellNum$cellName))
MangLocs<-cbind(Mang_region_location, Mang_cell_LatLong)
colnames(MangLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(MangLocs)
remove(Mang_cellNum)
remove(Mang_cell_LatLong)
uniquez<-MangLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(MangLocs, "Mang/MangLocs.csv", row.names = FALSE)





MBF_region_location<-read.csv("MBF/MBF_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
MBF_cellNum <- as.data.frame(cellFromXY(envRast, MBF_region_location[-1]))
colnames(MBF_cellNum)<-"cellName"
head(MBF_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

MBF_cell_LatLong <- as.data.frame(xyFromCell(envRast, MBF_cellNum$cellName))
MBFLocs<-cbind(MBF_region_location, MBF_cell_LatLong)
colnames(MBFLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(MBFLocs)
remove(MBF_cellNum)
remove(MBF_cell_LatLong)
uniquez<-MBFLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(MBFLocs, "MBF/MBFLocs.csv", row.names = FALSE)





Medit_region_location<-read.csv("Medit/Medit_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
Medit_cellNum <- as.data.frame(cellFromXY(envRast, Medit_region_location[-1]))
colnames(Medit_cellNum)<-"cellName"
head(Medit_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

Medit_cell_LatLong <- as.data.frame(xyFromCell(envRast, Medit_cellNum$cellName))
MeditLocs<-cbind(Medit_region_location, Medit_cell_LatLong)
colnames(MeditLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(MeditLocs)
remove(Medit_cellNum)
remove(Medit_cell_LatLong)
uniquez<-MeditLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(MeditLocs, "Medit/MeditLocs.csv", row.names = FALSE)





MonG_region_location<-read.csv("MonG/MonG_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
MonG_cellNum <- as.data.frame(cellFromXY(envRast, MonG_region_location[-1]))
colnames(MonG_cellNum)<-"cellName"
head(MonG_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

MonG_cell_LatLong <- as.data.frame(xyFromCell(envRast, MonG_cellNum$cellName))
MonGLocs<-cbind(MonG_region_location, MonG_cell_LatLong)
colnames(MonGLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(MonGLocs)
remove(MonG_cellNum)
remove(MonG_cell_LatLong)
uniquez<-MonGLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(MonGLocs, "MonG/MonGLocs.csv", row.names = FALSE)





Tai_region_location<-read.csv("Tai/Tai_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
Tai_cellNum <- as.data.frame(cellFromXY(envRast, Tai_region_location[-1]))
colnames(Tai_cellNum)<-"cellName"
head(Tai_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

Tai_cell_LatLong <- as.data.frame(xyFromCell(envRast, Tai_cellNum$cellName))
TaiLocs<-cbind(Tai_region_location, Tai_cell_LatLong)
colnames(TaiLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TaiLocs)
remove(Tai_cellNum)
remove(Tai_cell_LatLong)
uniquez<-TaiLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TaiLocs, "Tai/TaiLocs.csv", row.names = FALSE)





TemBF_region_location<-read.csv("TemBF/TemBF_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
TemBF_cellNum <- as.data.frame(cellFromXY(envRast, TemBF_region_location[-1]))
colnames(TemBF_cellNum)<-"cellName"
head(TemBF_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

TemBF_cell_LatLong <- as.data.frame(xyFromCell(envRast, TemBF_cellNum$cellName))
TemBFLocs<-cbind(TemBF_region_location, TemBF_cell_LatLong)
colnames(TemBFLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TemBFLocs)
remove(TemBF_cellNum)
remove(TemBF_cell_LatLong)
uniquez<-TemBFLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TemBFLocs, "TemBF/TemBFLocs.csv", row.names = FALSE)





TemCF_region_location<-read.csv("TemCF/TemCF_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
TemCF_cellNum <- as.data.frame(cellFromXY(envRast, TemCF_region_location[-1]))
colnames(TemCF_cellNum)<-"cellName"
head(TemCF_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

TemCF_cell_LatLong <- as.data.frame(xyFromCell(envRast, TemCF_cellNum$cellName))
TemCFLocs<-cbind(TemCF_region_location, TemCF_cell_LatLong)
colnames(TemCFLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TemCFLocs)
remove(TemCF_cellNum)
remove(TemCF_cell_LatLong)
uniquez<-TemCFLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TemCFLocs, "TemCF/TemCFLocs.csv", row.names = FALSE)





TemG_region_location<-read.csv("TemG/TemG_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
TemG_cellNum <- as.data.frame(cellFromXY(envRast, TemG_region_location[-1]))
colnames(TemG_cellNum)<-"cellName"
head(TemG_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

TemG_cell_LatLong <- as.data.frame(xyFromCell(envRast, TemG_cellNum$cellName))
TemGLocs<-cbind(TemG_region_location, TemG_cell_LatLong)
colnames(TemGLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TemGLocs)
remove(TemG_cellNum)
remove(TemG_cell_LatLong)
uniquez<-TemGLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TemGLocs, "TemG/TemGLocs.csv", row.names = FALSE)





TroCF_region_location<-read.csv("TroCF/TroCF_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
TroCF_cellNum <- as.data.frame(cellFromXY(envRast, TroCF_region_location[-1]))
colnames(TroCF_cellNum)<-"cellName"
head(TroCF_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

TroCF_cell_LatLong <- as.data.frame(xyFromCell(envRast, TroCF_cellNum$cellName))
TroCFLocs<-cbind(TroCF_region_location, TroCF_cell_LatLong)
colnames(TroCFLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TroCFLocs)
remove(TroCF_cellNum)
remove(TroCF_cell_LatLong)
uniquez<-TroCFLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TroCFLocs, "TroCF/TroCFLocs.csv", row.names = FALSE)





TroG_region_location<-read.csv("TroG/TroG_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
TroG_cellNum <- as.data.frame(cellFromXY(envRast, TroG_region_location[-1]))
colnames(TroG_cellNum)<-"cellName"
head(TroG_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

TroG_cell_LatLong <- as.data.frame(xyFromCell(envRast, TroG_cellNum$cellName))
TroGLocs<-cbind(TroG_region_location, TroG_cell_LatLong)
colnames(TroGLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TroGLocs)
remove(TroG_cellNum)
remove(TroG_cell_LatLong)
uniquez<-TroGLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TroGLocs, "TroG/TroGLocs.csv", row.names = FALSE)





Tun_region_location<-read.csv("Tun/Tun_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
Tun_cellNum <- as.data.frame(cellFromXY(envRast, Tun_region_location[-1]))
colnames(Tun_cellNum)<-"cellName"
head(Tun_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

Tun_cell_LatLong <- as.data.frame(xyFromCell(envRast, Tun_cellNum$cellName))
TunLocs<-cbind(Tun_region_location, Tun_cell_LatLong)
colnames(TunLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(TunLocs)
remove(Tun_cellNum)
remove(Tun_cell_LatLong)
uniquez<-TunLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(TunLocs, "Tun/TunLocs.csv", row.names = FALSE)




All_region_location<-read.csv("All/All_region_location.csv")
envRast<-raster::brick("Rasters/envRast.tif")
All_cellNum <- as.data.frame(cellFromXY(envRast, All_region_location[-1]))
colnames(All_cellNum)<-"cellName"
head(All_cellNum)
# this is the function gdm calls on to assign cells from rasters to the nearest
# centroid lat/long - input these into region_location

All_cell_LatLong <- as.data.frame(xyFromCell(envRast, All_cellNum$cellName))
AllLocs<-cbind(All_region_location, All_cell_LatLong)
colnames(AllLocs)[c(4,5)]<-c("cellLong", "cellLat")
head(AllLocs)
remove(All_cellNum)
remove(All_cell_LatLong)
uniquez<-AllLocs %>% distinct(cellLong, cellLat) #great, we have 941 unique combinations
# of lats + longs
write.csv(AllLocs, "All/AllLocs.csv", row.names = FALSE)









 
## Canaper PD and RPD + plots ----


# Bind together the matrices
MatBind<-bind_rows(DBFMat, DesMat, FloGMat, MangMat, MBFMat, MeditMat, MonGMat, 
                   TaiMat, TemBFMat, TemCFMat, TemGMat, TroCFMat, TroGMat, TunMat)

row_obj<-(c(paste("DBF", DBFMat$grids, sep = "_"), paste("Des", DesMat$grids, sep = "_"),
            paste("FloG", FloGMat$grids, sep = "_"), paste("Mang", MangMat$grids, sep = "_"),
            paste("MBF", MBFMat$grids, sep = "_"), paste("Medit", MeditMat$grids, sep = "_"),
            paste("MonG", MonGMat$grids, sep = "_"), paste("Tai", TaiMat$grids, sep = "_"),
            paste("TemBF", TemBFMat$grids, sep = "_"), paste("TemCF", TemCFMat$grids, sep = "_"),
            paste("TemG", TemGMat$grids, sep = "_"), paste("TroCF", TroCFMat$grids, sep = "_"),
            paste("TroG", TroGMat$grids, sep = "_"), paste("Tun", TunMat$grids, sep = "_")))

remove(DBFMat, DesMat, FloGMat, MangMat, MBFMat, MeditMat, MonGMat, TaiMat, 
       TemBFMat, TemCFMat, TemGMat, TroCFMat, TroGMat, TunMat)

row.names(MatBind)[1:5]
MatBind<-MatBind[,-15853] #change this to the "grids" col

MatBind<-as.data.frame(lapply(MatBind, as.numeric))
MatBind[is.na(MatBind)] <- 0 #bind_rows() fills in any species missing from a grid cell
#with NA when it binds, change those to 0
MatBind[1:5, 1:5]
write.csv(row_obj, "MatBind_rownames.csv", row.names= FALSE)
write.csv(MatBind, "MatBind.csv", row.names = FALSE)


rownames(MatBind) <- row_obj[,1]
MatBind[1:5,1:5]
SnB_updated <- drop.tip(SnB_updated, which(is.na(match(SnB_updated$tip.label, names(MatBind)))))
SnB_updated
MatBind <- MatBind[,-which(is.na(match(names(MatBind), SnB_updated$tip.label)))]
dim(MatBind)
MatBind <- MatBind[-which(rowSums(MatBind)<1),]

# We will run Canaper in parallel

# Set up the cluster
set.seed(12345)
availableCores()
# Configure parallel workers
usethis::edit_r_environ()
gc()
plan(multisession, workers = 30) # Using 30 workers
set.seed(071421)
# Max size of objects that can be passed between workers.
options(future.globals.maxSize= 8912896000)

tic() # Set a timer

# Run randomization test
rand_res <- cpr_rand_test(
  MatBind, SnB_updated,
  null_model = "curveball",
  n_reps = 100, n_iterations = 100000,
  tbl_out = TRUE,
  site_col = "site"
)
toc() # Stop the timer
plan(sequential) # End the parallel
head(rand_res)

# use cpr to classify the significance of each metric
mypd <-
  cpr_classify_endem(rand_res) |>
  cpr_classify_signif("pd") |>
  cpr_classify_signif("rpd") |>
  cpr_classify_signif("pe") |>
  cpr_classify_signif("rpe")

write.csv(mypd,"CANAPER_results.csv")



# Create a shapefile of canaper for mapping in QGIS

# Assign unique grid cell names for each biome, just like the row_obj
DBFMap$row_name <-paste("DBF", DBFMap$grids, sep = "_")
head(DBFMap$row_name)
DesMap$row_name <-paste("Des", DesMap$grids, sep = "_")
head(DesMap$row_name)
FloGMap$row_name <-paste("FloG", FloGMap$grids, sep = "_")
head(FloGMap$row_name)
MangMap$row_name <-paste("Mang", MangMap$grids, sep = "_")
head(MangMap$row_name)
MBFMap$row_name <-paste("MBF", MBFMap$grids, sep = "_")
head(MBFMap$row_name)
MeditMap$row_name <-paste("Medit", MeditMap$grids, sep = "_")
head(MeditMap$row_name)
MonGMap$row_name <-paste("MonG", MonGMap$grids, sep = "_")
head(MonGMap$row_name)
TaiMap$row_name <-paste("Tai", TaiMap$grids, sep = "_")
head(TaiMap$row_name)
TemBFMap$row_name <-paste("TemBF", TemBFMap$grids, sep = "_")
head(TemBFMap$row_name)
TemCFMap$row_name <-paste("TemCF", TemCFMap$grids, sep = "_")
head(TemCFMap$row_name)
TemGMap$row_name <-paste("TemG", TemGMap$grids, sep = "_")
head(TemGMap$row_name)
TroCFMap$row_name <-paste("TroCF", TroCFMap$grids, sep = "_")
head(TroCFMap$row_name)
TroGMap$row_name <-paste("TroG", TroGMap$grids, sep = "_")
head(TroGMap$row_name)
TunMap$row_name <-paste("Tun", TunMap$grids, sep = "_")
head(TunMap$row_name)

# Bind all of the maps together
BiomeMap<-rbind(DBFMap, DesMap, FloGMap, MangMap, MBFMap, MeditMap, MonGMap, 
                TaiMap, TemBFMap, TemCFMap, TemGMap, TroCFMap, TroGMap, TunMap)


BiomeMap<-BiomeMap[which(BiomeMap$row_name %in% mypd$site),]
mypd1<-mypd[which(mypd$site %in% BiomeMap$row_name),]

# Initiate the columns we want in our new shapefile
BiomeMap$pd_obs<-NA
BiomeMap$rpd_obs<-NA
BiomeMap$pe_obs<-NA
BiomeMap$rpe_obs<-NA
BiomeMap$pd_signif<-NA
BiomeMap$rpd_signif<-NA
BiomeMap$pe_signif<-NA
BiomeMap$rpe_signif<-NA

# Assign the values to each column for each grid cell
for(i in 1:nrow(BiomeMap)){
  print(i)
  BiomeMap$pd_obs[i]<-mypd1$pd_obs[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$rpd_obs[i]<-mypd1$rpd_obs[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$pe_obs[i]<-mypd1$pe_obs[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$rpe_obs[i]<-mypd1$rpe_obs[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$pd_signif[i]<-mypd1$pd_signif[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$rpd_signif[i]<-mypd1$rpd_signif[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$pe_signif[i]<-mypd1$pe_signif[which(mypd1$site == BiomeMap$row_name[i])]
  BiomeMap$rpe_signif[i]<-mypd1$rpe_signif[which(mypd1$site == BiomeMap$row_name[i])]
}

BiomeMap$biomeabbr<-sub("_.*", "", BiomeMap$row_name)
BiomeMap$biomeabbr[1:45]
write_sf(BiomeMap, "BiomeMap_1.shp", append = FALSE)



# Plot the PD, RPD, and significances
mypd<-mypd[,-1]
mypd$biomeabbr<-sub("_.*", "", mypd$site)

# PD plot
fit<-kruskal.test(pd_obs~biomeabbr, data = mypd)
print(fit)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result <- dunn.test(mypd$pd_obs, mypd$biomeabbr, method = "bonferroni")
print(dunn_result)

# Assign labels for results that are significant
sigs<-dunn_result$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)

# Plotting PD
Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")
maxes<-mypd %>% group_by(biomeabbr) %>%
  summarise(max = max(pd_obs))


quartz(w=3.375, h=2.5)

ggplot(mypd, aes(x = biomeabbr, y = pd_obs, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(mypd$pd_obs), 1.05 * max(mypd$pd_obs))) +
  labs(y = "PD", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.03*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)


#RPD

fit<-kruskal.test(rpd_obs~biomeabbr, data = mypd)
print(fit)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result <- dunn.test(mypd$rpd_obs, mypd$biomeabbr, method = "bonferroni")
print(dunn_result)

# Assign labels for results that are significant
sigs<-dunn_result$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)

# Plotting PD
Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")
maxes<-mypd %>% group_by(biomeabbr) %>%
  summarise(max = max(rpd_obs))


quartz(w=3.375, h=2.5)

ggplot(mypd, aes(x = biomeabbr, y = rpd_obs, fill = biomeabbr)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4) +
  scale_y_continuous(limits = c(min(mypd$rpd_obs), 1.05 * max(mypd$rpd_obs))) +
  labs(y = "RPD", x = "Biome Abbreviation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 8),
    axis.title.x = element_text(family = "Times"),
    axis.title.y = element_text(family = "Times"),
    panel.grid = element_blank(),     # Removes the gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Adds an outline
    axis.ticks = element_line(color = "black"),  # Keeps tick marks
    legend.position = "none"
  ) +
  scale_fill_manual(values = Hexxs) +
  geom_text(data = biome_labels, aes(x = biomeabbr, y = 1.03*maxes$max, label = label), 
            family = "Times", size = 3, vjust = 0)


#PD_signif
fit<-kruskal.test(pd_signif~biomeabbr, data = mypd)
print(fit)

# Significant result, p < 2.2e-16, so now do a post-hoc Dunn's test to compare all pairs of groups
dunn_result <- dunn.test(mypd$pd_signif, mypd$biomeabbr, method = "bonferroni")
print(dunn_result)

# Assign labels for results that are significant
sigs<-dunn_result$P.adjusted<0.05
Names<-gsub(" ", "", dunn_result$comparisons)
names(sigs)<-Names
LABELS<-multcompLetters(sigs)
biome_labels <- data.frame(biomeabbr = names(LABELS$Letters), label = LABELS$Letters)



# Plotting significance categories

# Plot PD significance
dat_pd<-data.frame(biomeabbr = mypd$biomeabbr, 
                category = mypd$pd_signif)

# Calculate the proportion of each PD significance category for each biome
prop_dat_pd<- dat_pd %>% 
  group_by(biomeabbr, category) %>% 
  summarise(count = n(), groups = "drop") %>% # Count each category for each biome
  group_by(biomeabbr) %>% 
  mutate(prop = count/sum(count)) %>%  # calculate the proportion of each category for each biome
  ungroup()

# Color palette for our categories
colorz<-c("#882255", "#AA4499", "#DDCC77", "#88CCEE", "#332288" )

quartz(w = 3.375, h = 2.5)

# Order the categories for display in plot
prop_dat_pd$category <- factor(prop_dat_pd$category, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))

# Create stacked bar plot showing proportion of PD signif. categories
ggplot(prop_dat_pd, aes(x = biomeabbr, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Biome Abbreviation", y = "Proportion of PD signif.", fill = "category") +
  scale_fill_manual(values = colorz) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Times"), 
        legend.position = "none")



# Repeat the above plot for RPD

dat_rpd<-data.frame(biomeabbr = mypd$biomeabbr, 
                category = mypd$rpd_signif)

prop_dat_rpd<- dat_rpd %>% 
  group_by(biomeabbr, category) %>% 
  summarise(count = n(), groups = "drop") %>% 
  group_by(biomeabbr) %>% 
  mutate(prop = count/sum(count)) %>% 
  ungroup()


quartz(w = 3.375, h = 2.5)

prop_dat_rpd$category <- factor(prop_dat_rpd$category, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))
ggplot(prop_dat_rpd, aes(x = biomeabbr, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Biome Abbreviation", y = "Proportion of RPD signif. category", fill = "category") +
  scale_fill_manual(values = colorz) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Times"), 
        legend.position = "none")






## Environmental matrices for logistic glms and slope regressions----

row.names(MatBind)<-row_obj$x
MatBind[1:5, 1:5]

# Update the phylogeny with the correct spp.
SnB_updated$tip.label<-gsub("_", ".", SnB_updated$tip.label)
SnB_updated$tip.label[1:5]
colnames(MatBind)[1:5]
MatBind1<-as.matrix(MatBind)
MatBind1<-dense2sparse(MatBind1)

# set up a cols object for simplicity
columns<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
           "s1.envRast_1", "s2.envRast_1", "s1.envRast_12", "s2.envRast_12", "s1.envRast_22",
           "s2.envRast_22", "s2.envRast_23", "s2.envRast_23")



# And begin with DBF

DBFEnv<-DBF_phylogdmTab[, columns]

# Save a df with a subset of our environmental variables
colnames(DBFEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
DBFEnv$biomeabbr<-"DBF"

# Here is when we have to re-assign the grid cell names, so we will use the 
  # GDM and matching biome_locs dfs to do that

#Make sure both dfs have the same coordinate decimals
DBFEnv$s1.xCoord<-round(DBFEnv$s1.xCoord, 2)
DBFEnv$s2.xCoord<-round(DBFEnv$s2.xCoord, 2)
DBFEnv$s1.yCoord<-round(DBFEnv$s1.yCoord, 2)
DBFEnv$s2.yCoord<-round(DBFEnv$s2.yCoord, 2)
DBF_locs$cellLat<-round(DBF_locs$cellLat, 2)
DBF_locs$cellLong<-round(DBF_locs$cellLong, 2)


# Write a for loop to assign the sites
DBFEnv$site1<-NA
DBFEnv$site2<-NA
for(i in 1:nrow(DBFEnv)){
  print(i)
  DBFEnv$site1[i]<-DBF_locs$site[(which((DBFEnv$s1.xCoord[i] == DBF_locs$cellLong) & 
                                          (DBFEnv$s1.yCoord[i] == DBF_locs$cellLat)))]
  DBFEnv$site2[i]<-DBF_locs$site[(which((DBFEnv$s2.xCoord[i] == DBF_locs$cellLong) & 
                                          (DBFEnv$s2.yCoord[i] == DBF_locs$cellLat)))]
}

# Assign the correct values for each of the site's predictor variables
DBFEnv$bio1Mean<-NA
DBFEnv$bio12Mean<-NA
DBFEnv$pStabMean<-NA
DBFEnv$tStabMean<-NA
for(i in 1:nrow(DBFEnv)){
  print(i)
  DBFEnv$bio1Mean[i]<-mean(c(DBFEnv$s1.bio1[i], DBFEnv$s2.bio1[i]))
  DBFEnv$bio12Mean[i]<-mean(c(DBFEnv$s1.bio12[i], DBFEnv$s2.bio12[i]))
  DBFEnv$pStabMean[i]<-mean(c(DBFEnv$s1.precipStability[i], DBFEnv$s2.precipStability[i]))
  DBFEnv$tStabMean[i]<-mean(c(DBFEnv$s1.tempStability[i], DBFEnv$s2.tempStability[i]))
}

DBFEnv$geoDist<-distHaversine(matrix(c(DBFEnv$s1.xCoord, DBFEnv$s1.yCoord), ncol = 2),
                                 matrix(c(DBFEnv$s2.xCoord, DBFEnv$s2.yCoord), ncol = 2))/1000

write.csv(DBFEnv, "DBF/DBFEnv.csv", row.names = FALSE)



# Repeat across other biomes



Des_phylogdmTab<-as.data.frame(fread("Des/Des_phylogdmTab.csv"))
DesEnv<-Des_phylogdmTab[, columns]
colnames(DesEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
DesEnv$biomeabbr<-"Des"
Des_locs<-as.data.frame(fread("Des/DesLocs.csv"))
DesEnv$s1.xCoord<-round(DesEnv$s1.xCoord, 2)
DesEnv$s2.xCoord<-round(DesEnv$s2.xCoord, 2)
DesEnv$s1.yCoord<-round(DesEnv$s1.yCoord, 2)
DesEnv$s2.yCoord<-round(DesEnv$s2.yCoord, 2)
Des_locs$cellLat<-round(Des_locs$cellLat, 2)
Des_locs$cellLong<-round(Des_locs$cellLong, 2)

DesEnv$site1<-NA
DesEnv$site2<-NA
for(i in 1:nrow(DesEnv)){
  print(i)
  DesEnv$site1[i]<-Des_locs$site[(which((DesEnv$s1.xCoord[i] == Des_locs$cellLong) & 
                                          (DesEnv$s1.yCoord[i] == Des_locs$cellLat)))]
  DesEnv$site2[i]<-Des_locs$site[(which((DesEnv$s2.xCoord[i] == Des_locs$cellLong) & 
                                          (DesEnv$s2.yCoord[i] == Des_locs$cellLat)))]
}

DesEnv$bio1Mean<-NA
DesEnv$bio12Mean<-NA
DesEnv$pStabMean<-NA
DesEnv$tStabMean<-NA
for(i in 1:nrow(DesEnv)){
  print(i)
  DesEnv$bio1Mean[i]<-mean(c(DesEnv$s1.bio1[i], DesEnv$s2.bio1[i]))
  DesEnv$bio12Mean[i]<-mean(c(DesEnv$s1.bio12[i], DesEnv$s2.bio12[i]))
  DesEnv$pStabMean[i]<-mean(c(DesEnv$s1.precipStability[i], DesEnv$s2.precipStability[i]))
  DesEnv$tStabMean[i]<-mean(c(DesEnv$s1.tempStability[i], DesEnv$s2.tempStability[i]))
}

DesEnv$geoDist<-distHaversine(matrix(c(DesEnv$s1.xCoord, DesEnv$s1.yCoord), ncol = 2),
                              matrix(c(DesEnv$s2.xCoord, DesEnv$s2.yCoord), ncol = 2))/1000

write.csv(DesEnv, "Des/DesEnv.csv", row.names = FALSE)







FloG_phylogdmTab<-as.data.frame(fread("FloG/FloG_phylogdmTab.csv"))
FloGEnv<-FloG_phylogdmTab[, columns]
colnames(FloGEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                     "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                     "s2.precipStability", "s1.tempStability", "s2.tempStability")
FloGEnv$biomeabbr<-"FloG"
FloG_locs<-as.data.frame(fread("FloG/FloGLocs.csv"))
FloGEnv$s1.xCoord<-round(FloGEnv$s1.xCoord, 2)
FloGEnv$s2.xCoord<-round(FloGEnv$s2.xCoord, 2)
FloGEnv$s1.yCoord<-round(FloGEnv$s1.yCoord, 2)
FloGEnv$s2.yCoord<-round(FloGEnv$s2.yCoord, 2)
FloG_locs$cellLat<-round(FloG_locs$cellLat, 2)
FloG_locs$cellLong<-round(FloG_locs$cellLong, 2)

FloGEnv$site1<-NA
FloGEnv$site2<-NA
for(i in 1:nrow(FloGEnv)){
  print(i)
  FloGEnv$site1[i]<-FloG_locs$site[(which((FloGEnv$s1.xCoord[i] == FloG_locs$cellLong) & 
                                            (FloGEnv$s1.yCoord[i] == FloG_locs$cellLat)))]
  FloGEnv$site2[i]<-FloG_locs$site[(which((FloGEnv$s2.xCoord[i] == FloG_locs$cellLong) & 
                                            (FloGEnv$s2.yCoord[i] == FloG_locs$cellLat)))]
}

FloGEnv$bio1Mean<-NA
FloGEnv$bio12Mean<-NA
FloGEnv$pStabMean<-NA
FloGEnv$tStabMean<-NA
for(i in 1:nrow(FloGEnv)){
  print(i)
  FloGEnv$bio1Mean[i]<-mean(c(FloGEnv$s1.bio1[i], FloGEnv$s2.bio1[i]))
  FloGEnv$bio12Mean[i]<-mean(c(FloGEnv$s1.bio12[i], FloGEnv$s2.bio12[i]))
  FloGEnv$pStabMean[i]<-mean(c(FloGEnv$s1.precipStability[i], FloGEnv$s2.precipStability[i]))
  FloGEnv$tStabMean[i]<-mean(c(FloGEnv$s1.tempStability[i], FloGEnv$s2.tempStability[i]))
}

FloGEnv$geoDist<-distHaversine(matrix(c(FloGEnv$s1.xCoord, FloGEnv$s1.yCoord), ncol = 2),
                              matrix(c(FloGEnv$s2.xCoord, FloGEnv$s2.yCoord), ncol = 2))/1000

write.csv(FloGEnv, "FloG/FloGEnv.csv", row.names = FALSE)








Mang_phylogdmTab<-as.data.frame(fread("Mang/Mang_phylogdmTab.csv"))
MangEnv<-Mang_phylogdmTab[, columns]
colnames(MangEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                     "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                     "s2.precipStability", "s1.tempStability", "s2.tempStability")
MangEnv$biomeabbr<-"Mang"
Mang_locs<-as.data.frame(fread("Mang/MangLocs.csv"))
MangEnv$s1.xCoord<-round(MangEnv$s1.xCoord, 2)
MangEnv$s2.xCoord<-round(MangEnv$s2.xCoord, 2)
MangEnv$s1.yCoord<-round(MangEnv$s1.yCoord, 2)
MangEnv$s2.yCoord<-round(MangEnv$s2.yCoord, 2)
Mang_locs$cellLat<-round(Mang_locs$cellLat, 2)
Mang_locs$cellLong<-round(Mang_locs$cellLong, 2)

MangEnv$site1<-NA
MangEnv$site2<-NA
for(i in 1:nrow(MangEnv)){
  print(i)
  MangEnv$site1[i]<-Mang_locs$site[(which((MangEnv$s1.xCoord[i] == Mang_locs$cellLong) & 
                                            (MangEnv$s1.yCoord[i] == Mang_locs$cellLat)))]
  MangEnv$site2[i]<-Mang_locs$site[(which((MangEnv$s2.xCoord[i] == Mang_locs$cellLong) & 
                                            (MangEnv$s2.yCoord[i] == Mang_locs$cellLat)))]
}

MangEnv$bio1Mean<-NA
MangEnv$bio12Mean<-NA
MangEnv$pStabMean<-NA
MangEnv$tStabMean<-NA
for(i in 1:nrow(MangEnv)){
  print(i)
  MangEnv$bio1Mean[i]<-mean(c(MangEnv$s1.bio1[i], MangEnv$s2.bio1[i]))
  MangEnv$bio12Mean[i]<-mean(c(MangEnv$s1.bio12[i], MangEnv$s2.bio12[i]))
  MangEnv$pStabMean[i]<-mean(c(MangEnv$s1.precipStability[i], MangEnv$s2.precipStability[i]))
  MangEnv$tStabMean[i]<-mean(c(MangEnv$s1.tempStability[i], MangEnv$s2.tempStability[i]))
}

MangEnv$geoDist<-distHaversine(matrix(c(MangEnv$s1.xCoord, MangEnv$s1.yCoord), ncol = 2),
                              matrix(c(MangEnv$s2.xCoord, MangEnv$s2.yCoord), ncol = 2))/1000

write.csv(MangEnv, "Mang/MangEnv.csv", row.names = FALSE)








MBF_phylogdmTab<-as.data.frame(fread("MBF/MBF_phylogdmTab.csv"))
MBFEnv<-MBF_phylogdmTab[, columns]
colnames(MBFEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
MBFEnv$biomeabbr<-"MBF"
MBF_locs<-as.data.frame(fread("MBF/MBFLocs.csv"))
MBFEnv$s1.xCoord<-round(MBFEnv$s1.xCoord, 2)
MBFEnv$s2.xCoord<-round(MBFEnv$s2.xCoord, 2)
MBFEnv$s1.yCoord<-round(MBFEnv$s1.yCoord, 2)
MBFEnv$s2.yCoord<-round(MBFEnv$s2.yCoord, 2)
MBF_locs$cellLat<-round(MBF_locs$cellLat, 2)
MBF_locs$cellLong<-round(MBF_locs$cellLong, 2)

MBFEnv$site1<-NA
MBFEnv$site2<-NA
for(i in 1:nrow(MBFEnv)){
  print(i)
  MBFEnv$site1[i]<-MBF_locs$site[(which((MBFEnv$s1.xCoord[i] == MBF_locs$cellLong) & 
                                          (MBFEnv$s1.yCoord[i] == MBF_locs$cellLat)))]
  MBFEnv$site2[i]<-MBF_locs$site[(which((MBFEnv$s2.xCoord[i] == MBF_locs$cellLong) & 
                                          (MBFEnv$s2.yCoord[i] == MBF_locs$cellLat)))]
}

MBFEnv$bio1Mean<-NA
MBFEnv$bio12Mean<-NA
MBFEnv$pStabMean<-NA
MBFEnv$tStabMean<-NA
for(i in 1:nrow(MBFEnv)){
  print(i)
  MBFEnv$bio1Mean[i]<-mean(c(MBFEnv$s1.bio1[i], MBFEnv$s2.bio1[i]))
  MBFEnv$bio12Mean[i]<-mean(c(MBFEnv$s1.bio12[i], MBFEnv$s2.bio12[i]))
  MBFEnv$pStabMean[i]<-mean(c(MBFEnv$s1.precipStability[i], MBFEnv$s2.precipStability[i]))
  MBFEnv$tStabMean[i]<-mean(c(MBFEnv$s1.tempStability[i], MBFEnv$s2.tempStability[i]))
}

MBFEnv$geoDist<-distHaversine(matrix(c(MBFEnv$s1.xCoord, MBFEnv$s1.yCoord), ncol = 2),
                              matrix(c(MBFEnv$s2.xCoord, MBFEnv$s2.yCoord), ncol = 2))/1000

write.csv(MBFEnv, "MBF/MBFEnv.csv", row.names = FALSE)












Medit_phylogdmTab<-as.data.frame(fread("Medit/Medit_phylogdmTab.csv"))
MeditEnv<-Medit_phylogdmTab[, columns]
colnames(MeditEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                      "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                      "s2.precipStability", "s1.tempStability", "s2.tempStability")
MeditEnv$biomeabbr<-"Medit"
Medit_locs<-as.data.frame(fread("Medit/MeditLocs.csv"))
MeditEnv$s1.xCoord<-round(MeditEnv$s1.xCoord, 2)
MeditEnv$s2.xCoord<-round(MeditEnv$s2.xCoord, 2)
MeditEnv$s1.yCoord<-round(MeditEnv$s1.yCoord, 2)
MeditEnv$s2.yCoord<-round(MeditEnv$s2.yCoord, 2)
Medit_locs$cellLat<-round(Medit_locs$cellLat, 2)
Medit_locs$cellLong<-round(Medit_locs$cellLong, 2)

MeditEnv$site1<-NA
MeditEnv$site2<-NA
for(i in 1:nrow(MeditEnv)){
  print(i)
  MeditEnv$site1[i]<-Medit_locs$site[(which((MeditEnv$s1.xCoord[i] == Medit_locs$cellLong) & 
                                              (MeditEnv$s1.yCoord[i] == Medit_locs$cellLat)))]
  MeditEnv$site2[i]<-Medit_locs$site[(which((MeditEnv$s2.xCoord[i] == Medit_locs$cellLong) & 
                                              (MeditEnv$s2.yCoord[i] == Medit_locs$cellLat)))]
}

MeditEnv$bio1Mean<-NA
MeditEnv$bio12Mean<-NA
MeditEnv$pStabMean<-NA
MeditEnv$tStabMean<-NA
for(i in 1:nrow(MeditEnv)){
  print(i)
  MeditEnv$bio1Mean[i]<-mean(c(MeditEnv$s1.bio1[i], MeditEnv$s2.bio1[i]))
  MeditEnv$bio12Mean[i]<-mean(c(MeditEnv$s1.bio12[i], MeditEnv$s2.bio12[i]))
  MeditEnv$pStabMean[i]<-mean(c(MeditEnv$s1.precipStability[i], MeditEnv$s2.precipStability[i]))
  MeditEnv$tStabMean[i]<-mean(c(MeditEnv$s1.tempStability[i], MeditEnv$s2.tempStability[i]))
}

MeditEnv$geoDist<-distHaversine(matrix(c(MeditEnv$s1.xCoord, MeditEnv$s1.yCoord), ncol = 2),
                              matrix(c(MeditEnv$s2.xCoord, MeditEnv$s2.yCoord), ncol = 2))/1000


write.csv(MeditEnv, "Medit/MeditEnv.csv", row.names = FALSE)










MonG_phylogdmTab<-as.data.frame(fread("MonG/MonG_phylogdmTab.csv"))
MonGEnv<-MonG_phylogdmTab[, columns]
colnames(MonGEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                     "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                     "s2.precipStability", "s1.tempStability", "s2.tempStability")
MonGEnv$biomeabbr<-"MonG"
MonG_locs<-as.data.frame(fread("MonG/MonGLocs.csv"))
MonGEnv$s1.xCoord<-round(MonGEnv$s1.xCoord, 2)
MonGEnv$s2.xCoord<-round(MonGEnv$s2.xCoord, 2)
MonGEnv$s1.yCoord<-round(MonGEnv$s1.yCoord, 2)
MonGEnv$s2.yCoord<-round(MonGEnv$s2.yCoord, 2)
MonG_locs$cellLat<-round(MonG_locs$cellLat, 2)
MonG_locs$cellLong<-round(MonG_locs$cellLong, 2)

MonGEnv$site1<-NA
MonGEnv$site2<-NA
for(i in 1:nrow(MonGEnv)){
  print(i)
  MonGEnv$site1[i]<-MonG_locs$site[(which((MonGEnv$s1.xCoord[i] == MonG_locs$cellLong) & 
                                            (MonGEnv$s1.yCoord[i] == MonG_locs$cellLat)))]
  MonGEnv$site2[i]<-MonG_locs$site[(which((MonGEnv$s2.xCoord[i] == MonG_locs$cellLong) & 
                                            (MonGEnv$s2.yCoord[i] == MonG_locs$cellLat)))]
}

MonGEnv$bio1Mean<-NA
MonGEnv$bio12Mean<-NA
MonGEnv$pStabMean<-NA
MonGEnv$tStabMean<-NA
for(i in 1:nrow(MonGEnv)){
  print(i)
  MonGEnv$bio1Mean[i]<-mean(c(MonGEnv$s1.bio1[i], MonGEnv$s2.bio1[i]))
  MonGEnv$bio12Mean[i]<-mean(c(MonGEnv$s1.bio12[i], MonGEnv$s2.bio12[i]))
  MonGEnv$pStabMean[i]<-mean(c(MonGEnv$s1.precipStability[i], MonGEnv$s2.precipStability[i]))
  MonGEnv$tStabMean[i]<-mean(c(MonGEnv$s1.tempStability[i], MonGEnv$s2.tempStability[i]))
}

MonGEnv$geoDist<-distHaversine(matrix(c(MonGEnv$s1.xCoord, MonGEnv$s1.yCoord), ncol = 2),
                              matrix(c(MonGEnv$s2.xCoord, MonGEnv$s2.yCoord), ncol = 2))/1000


write.csv(MonGEnv, "MonG/MonGEnv.csv", row.names = FALSE)









Tai_phylogdmTab<-as.data.frame(fread("Tai/Tai_phylogdmTab.csv"))
TaiEnv<-Tai_phylogdmTab[, columns]
colnames(TaiEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
TaiEnv$biomeabbr<-"Tai"
Tai_locs<-as.data.frame(fread("Tai/TaiLocs.csv"))
TaiEnv$s1.xCoord<-round(TaiEnv$s1.xCoord, 2)
TaiEnv$s2.xCoord<-round(TaiEnv$s2.xCoord, 2)
TaiEnv$s1.yCoord<-round(TaiEnv$s1.yCoord, 2)
TaiEnv$s2.yCoord<-round(TaiEnv$s2.yCoord, 2)
Tai_locs$cellLat<-round(Tai_locs$cellLat, 2)
Tai_locs$cellLong<-round(Tai_locs$cellLong, 2)

TaiEnv$site1<-NA
TaiEnv$site2<-NA
for(i in 1:nrow(TaiEnv)){
  print(i)
  TaiEnv$site1[i]<-Tai_locs$site[(which((TaiEnv$s1.xCoord[i] == Tai_locs$cellLong) & 
                                          (TaiEnv$s1.yCoord[i] == Tai_locs$cellLat)))]
  TaiEnv$site2[i]<-Tai_locs$site[(which((TaiEnv$s2.xCoord[i] == Tai_locs$cellLong) & 
                                          (TaiEnv$s2.yCoord[i] == Tai_locs$cellLat)))]
}

TaiEnv$bio1Mean<-NA
TaiEnv$bio12Mean<-NA
TaiEnv$pStabMean<-NA
TaiEnv$tStabMean<-NA
for(i in 1:nrow(TaiEnv)){
  print(i)
  TaiEnv$bio1Mean[i]<-mean(c(TaiEnv$s1.bio1[i], TaiEnv$s2.bio1[i]))
  TaiEnv$bio12Mean[i]<-mean(c(TaiEnv$s1.bio12[i], TaiEnv$s2.bio12[i]))
  TaiEnv$pStabMean[i]<-mean(c(TaiEnv$s1.precipStability[i], TaiEnv$s2.precipStability[i]))
  TaiEnv$tStabMean[i]<-mean(c(TaiEnv$s1.tempStability[i], TaiEnv$s2.tempStability[i]))
}

TaiEnv$geoDist<-distHaversine(matrix(c(TaiEnv$s1.xCoord, TaiEnv$s1.yCoord), ncol = 2),
                              matrix(c(TaiEnv$s2.xCoord, TaiEnv$s2.yCoord), ncol = 2))/1000

write.csv(TaiEnv, "Tai/TaiEnv.csv", row.names = FALSE)










TemBF_phylogdmTab<-as.data.frame(fread("TemBF/TemBF_phylogdmTab.csv"))
TemBFEnv<-TemBF_phylogdmTab[, columns]
colnames(TemBFEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                      "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                      "s2.precipStability", "s1.tempStability", "s2.tempStability")
TemBFEnv$biomeabbr<-"TemBF"
TemBF_locs<-as.data.frame(fread("TemBF/TemBFLocs.csv"))
TemBFEnv$s1.xCoord<-round(TemBFEnv$s1.xCoord, 2)
TemBFEnv$s2.xCoord<-round(TemBFEnv$s2.xCoord, 2)
TemBFEnv$s1.yCoord<-round(TemBFEnv$s1.yCoord, 2)
TemBFEnv$s2.yCoord<-round(TemBFEnv$s2.yCoord, 2)
TemBF_locs$cellLat<-round(TemBF_locs$cellLat, 2)
TemBF_locs$cellLong<-round(TemBF_locs$cellLong, 2)

TemBFEnv$site1<-NA
TemBFEnv$site2<-NA
for(i in 1:nrow(TemBFEnv)){
  print(i)
  TemBFEnv$site1[i]<-TemBF_locs$site[(which((TemBFEnv$s1.xCoord[i] == TemBF_locs$cellLong) & 
                                              (TemBFEnv$s1.yCoord[i] == TemBF_locs$cellLat)))]
  TemBFEnv$site2[i]<-TemBF_locs$site[(which((TemBFEnv$s2.xCoord[i] == TemBF_locs$cellLong) & 
                                              (TemBFEnv$s2.yCoord[i] == TemBF_locs$cellLat)))]
}

TemBFEnv$bio1Mean<-NA
TemBFEnv$bio12Mean<-NA
TemBFEnv$pStabMean<-NA
TemBFEnv$tStabMean<-NA
for(i in 1:nrow(TemBFEnv)){
  print(i)
  TemBFEnv$bio1Mean[i]<-mean(c(TemBFEnv$s1.bio1[i], TemBFEnv$s2.bio1[i]))
  TemBFEnv$bio12Mean[i]<-mean(c(TemBFEnv$s1.bio12[i], TemBFEnv$s2.bio12[i]))
  TemBFEnv$pStabMean[i]<-mean(c(TemBFEnv$s1.precipStability[i], TemBFEnv$s2.precipStability[i]))
  TemBFEnv$tStabMean[i]<-mean(c(TemBFEnv$s1.tempStability[i], TemBFEnv$s2.tempStability[i]))
}


TemBFEnv$geoDist<-distHaversine(matrix(c(TemBFEnv$s1.xCoord, TemBFEnv$s1.yCoord), ncol = 2),
                              matrix(c(TemBFEnv$s2.xCoord, TemBFEnv$s2.yCoord), ncol = 2))/1000
write.csv(TemBFEnv, "TemBF/TemBFEnv.csv", row.names = FALSE)









TemCF_phylogdmTab<-as.data.frame(fread("TemCF/TemCF_phylogdmTab.csv"))
TemCFEnv<-TemCF_phylogdmTab[, columns]
colnames(TemCFEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                      "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                      "s2.precipStability", "s1.tempStability", "s2.tempStability")
TemCFEnv$biomeabbr<-"TemCF"
TemCF_locs<-as.data.frame(fread("TemCF/TemCFLocs.csv"))
TemCFEnv$s1.xCoord<-round(TemCFEnv$s1.xCoord, 2)
TemCFEnv$s2.xCoord<-round(TemCFEnv$s2.xCoord, 2)
TemCFEnv$s1.yCoord<-round(TemCFEnv$s1.yCoord, 2)
TemCFEnv$s2.yCoord<-round(TemCFEnv$s2.yCoord, 2)
TemCF_locs$cellLat<-round(TemCF_locs$cellLat, 2)
TemCF_locs$cellLong<-round(TemCF_locs$cellLong, 2)

TemCFEnv$site1<-NA
TemCFEnv$site2<-NA
for(i in 1:nrow(TemCFEnv)){
  print(i)
  TemCFEnv$site1[i]<-TemCF_locs$site[(which((TemCFEnv$s1.xCoord[i] == TemCF_locs$cellLong) & 
                                              (TemCFEnv$s1.yCoord[i] == TemCF_locs$cellLat)))]
  TemCFEnv$site2[i]<-TemCF_locs$site[(which((TemCFEnv$s2.xCoord[i] == TemCF_locs$cellLong) & 
                                              (TemCFEnv$s2.yCoord[i] == TemCF_locs$cellLat)))]
}

TemCFEnv$bio1Mean<-NA
TemCFEnv$bio12Mean<-NA
TemCFEnv$pStabMean<-NA
TemCFEnv$tStabMean<-NA
for(i in 1:nrow(TemCFEnv)){
  print(i)
  TemCFEnv$bio1Mean[i]<-mean(c(TemCFEnv$s1.bio1[i], TemCFEnv$s2.bio1[i]))
  TemCFEnv$bio12Mean[i]<-mean(c(TemCFEnv$s1.bio12[i], TemCFEnv$s2.bio12[i]))
  TemCFEnv$pStabMean[i]<-mean(c(TemCFEnv$s1.precipStability[i], TemCFEnv$s2.precipStability[i]))
  TemCFEnv$tStabMean[i]<-mean(c(TemCFEnv$s1.tempStability[i], TemCFEnv$s2.tempStability[i]))
}

TemCFEnv$geoDist<-distHaversine(matrix(c(TemCFEnv$s1.xCoord, TemCFEnv$s1.yCoord), ncol = 2),
                              matrix(c(TemCFEnv$s2.xCoord, TemCFEnv$s2.yCoord), ncol = 2))/1000

write.csv(TemCFEnv, "TemCF/TemCFEnv.csv", row.names = FALSE)










TemG_phylogdmTab<-as.data.frame(fread("TemG/TemG_phylogdmTab.csv"))
TemGEnv<-TemG_phylogdmTab[, columns]
colnames(TemGEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                     "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                     "s2.precipStability", "s1.tempStability", "s2.tempStability")
TemGEnv$biomeabbr<-"TemG"
TemG_locs<-as.data.frame(fread("TemG/TemGLocs.csv"))
TemGEnv$s1.xCoord<-round(TemGEnv$s1.xCoord, 2)
TemGEnv$s2.xCoord<-round(TemGEnv$s2.xCoord, 2)
TemGEnv$s1.yCoord<-round(TemGEnv$s1.yCoord, 2)
TemGEnv$s2.yCoord<-round(TemGEnv$s2.yCoord, 2)
TemG_locs$cellLat<-round(TemG_locs$cellLat, 2)
TemG_locs$cellLong<-round(TemG_locs$cellLong, 2)

TemGEnv$site1<-NA
TemGEnv$site2<-NA
for(i in 1:nrow(TemGEnv)){
  print(i)
  TemGEnv$site1[i]<-TemG_locs$site[(which((TemGEnv$s1.xCoord[i] == TemG_locs$cellLong) & 
                                            (TemGEnv$s1.yCoord[i] == TemG_locs$cellLat)))]
  TemGEnv$site2[i]<-TemG_locs$site[(which((TemGEnv$s2.xCoord[i] == TemG_locs$cellLong) & 
                                            (TemGEnv$s2.yCoord[i] == TemG_locs$cellLat)))]
}

TemGEnv$bio1Mean<-NA
TemGEnv$bio12Mean<-NA
TemGEnv$pStabMean<-NA
TemGEnv$tStabMean<-NA
for(i in 1:nrow(TemGEnv)){
  print(i)
  TemGEnv$bio1Mean[i]<-mean(c(TemGEnv$s1.bio1[i], TemGEnv$s2.bio1[i]))
  TemGEnv$bio12Mean[i]<-mean(c(TemGEnv$s1.bio12[i], TemGEnv$s2.bio12[i]))
  TemGEnv$pStabMean[i]<-mean(c(TemGEnv$s1.precipStability[i], TemGEnv$s2.precipStability[i]))
  TemGEnv$tStabMean[i]<-mean(c(TemGEnv$s1.tempStability[i], TemGEnv$s2.tempStability[i]))
}

TemGEnv$geoDist<-distHaversine(matrix(c(TemGEnv$s1.xCoord, TemGEnv$s1.yCoord), ncol = 2),
                              matrix(c(TemGEnv$s2.xCoord, TemGEnv$s2.yCoord), ncol = 2))/1000

write.csv(TemGEnv, "TemG/TemGEnv.csv", row.names = FALSE)









TroCF_phylogdmTab<-as.data.frame(fread("TroCF/TroCF_phylogdmTab.csv"))
TroCFEnv<-TroCF_phylogdmTab[, columns]
colnames(TroCFEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                      "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                      "s2.precipStability", "s1.tempStability", "s2.tempStability")
TroCFEnv$biomeabbr<-"TroCF"
TroCF_locs<-as.data.frame(fread("TroCF/TroCFLocs.csv"))
TroCFEnv$s1.xCoord<-round(TroCFEnv$s1.xCoord, 2)
TroCFEnv$s2.xCoord<-round(TroCFEnv$s2.xCoord, 2)
TroCFEnv$s1.yCoord<-round(TroCFEnv$s1.yCoord, 2)
TroCFEnv$s2.yCoord<-round(TroCFEnv$s2.yCoord, 2)
TroCF_locs$cellLat<-round(TroCF_locs$cellLat, 2)
TroCF_locs$cellLong<-round(TroCF_locs$cellLong, 2)

TroCFEnv$site1<-NA
TroCFEnv$site2<-NA
for(i in 1:nrow(TroCFEnv)){
  print(i)
  TroCFEnv$site1[i]<-TroCF_locs$site[(which((TroCFEnv$s1.xCoord[i] == TroCF_locs$cellLong) & 
                                              (TroCFEnv$s1.yCoord[i] == TroCF_locs$cellLat)))]
  TroCFEnv$site2[i]<-TroCF_locs$site[(which((TroCFEnv$s2.xCoord[i] == TroCF_locs$cellLong) & 
                                              (TroCFEnv$s2.yCoord[i] == TroCF_locs$cellLat)))]
}

TroCFEnv$bio1Mean<-NA
TroCFEnv$bio12Mean<-NA
TroCFEnv$pStabMean<-NA
TroCFEnv$tStabMean<-NA
for(i in 1:nrow(TroCFEnv)){
  print(i)
  TroCFEnv$bio1Mean[i]<-mean(c(TroCFEnv$s1.bio1[i], TroCFEnv$s2.bio1[i]))
  TroCFEnv$bio12Mean[i]<-mean(c(TroCFEnv$s1.bio12[i], TroCFEnv$s2.bio12[i]))
  TroCFEnv$pStabMean[i]<-mean(c(TroCFEnv$s1.precipStability[i], TroCFEnv$s2.precipStability[i]))
  TroCFEnv$tStabMean[i]<-mean(c(TroCFEnv$s1.tempStability[i], TroCFEnv$s2.tempStability[i]))
}

TroCFEnv$geoDist<-distHaversine(matrix(c(TroCFEnv$s1.xCoord, TroCFEnv$s1.yCoord), ncol = 2),
                              matrix(c(TroCFEnv$s2.xCoord, TroCFEnv$s2.yCoord), ncol = 2))/1000


write.csv(TroCFEnv, "TroCF/TroCFEnv.csv", row.names = FALSE)










TroG_phylogdmTab<-as.data.frame(fread("TroG/TroG_phylogdmTab.csv"))
TroGEnv<-TroG_phylogdmTab[, columns]
colnames(TroGEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                     "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                     "s2.precipStability", "s1.tempStability", "s2.tempStability")
TroGEnv$biomeabbr<-"TroG"
TroG_locs<-as.data.frame(fread("TroG/TroGLocs.csv"))
TroGEnv$s1.xCoord<-round(TroGEnv$s1.xCoord, 2)
TroGEnv$s2.xCoord<-round(TroGEnv$s2.xCoord, 2)
TroGEnv$s1.yCoord<-round(TroGEnv$s1.yCoord, 2)
TroGEnv$s2.yCoord<-round(TroGEnv$s2.yCoord, 2)
TroG_locs$cellLat<-round(TroG_locs$cellLat, 2)
TroG_locs$cellLong<-round(TroG_locs$cellLong, 2)

TroGEnv$site1<-NA
TroGEnv$site2<-NA
for(i in 1:nrow(TroGEnv)){
  print(i)
  TroGEnv$site1[i]<-TroG_locs$site[(which((TroGEnv$s1.xCoord[i] == TroG_locs$cellLong) & 
                                            (TroGEnv$s1.yCoord[i] == TroG_locs$cellLat)))]
  TroGEnv$site2[i]<-TroG_locs$site[(which((TroGEnv$s2.xCoord[i] == TroG_locs$cellLong) & 
                                            (TroGEnv$s2.yCoord[i] == TroG_locs$cellLat)))]
}

TroGEnv$bio1Mean<-NA
TroGEnv$bio12Mean<-NA
TroGEnv$pStabMean<-NA
TroGEnv$tStabMean<-NA
for(i in 1:nrow(TroGEnv)){
  print(i)
  TroGEnv$bio1Mean[i]<-mean(c(TroGEnv$s1.bio1[i], TroGEnv$s2.bio1[i]))
  TroGEnv$bio12Mean[i]<-mean(c(TroGEnv$s1.bio12[i], TroGEnv$s2.bio12[i]))
  TroGEnv$pStabMean[i]<-mean(c(TroGEnv$s1.precipStability[i], TroGEnv$s2.precipStability[i]))
  TroGEnv$tStabMean[i]<-mean(c(TroGEnv$s1.tempStability[i], TroGEnv$s2.tempStability[i]))
}

TroGEnv$geoDist<-distHaversine(matrix(c(TroGEnv$s1.xCoord, TroGEnv$s1.yCoord), ncol = 2),
                              matrix(c(TroGEnv$s2.xCoord, TroGEnv$s2.yCoord), ncol = 2))/1000


write.csv(TroGEnv, "TroG/TroGEnv.csv", row.names = FALSE)








Tun_phylogdmTab<-as.data.frame(fread("Tun/Tun_phylogdmTab.csv"))
TunEnv<-Tun_phylogdmTab[, columns]
colnames(TunEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
TunEnv$biomeabbr<-"Tun"
Tun_locs<-as.data.frame(fread("Tun/TunLocs.csv"))
TunEnv$s1.xCoord<-round(TunEnv$s1.xCoord, 2)
TunEnv$s2.xCoord<-round(TunEnv$s2.xCoord, 2)
TunEnv$s1.yCoord<-round(TunEnv$s1.yCoord, 2)
TunEnv$s2.yCoord<-round(TunEnv$s2.yCoord, 2)
Tun_locs$cellLat<-round(Tun_locs$cellLat, 2)
Tun_locs$cellLong<-round(Tun_locs$cellLong, 2)

TunEnv$site1<-NA
TunEnv$site2<-NA
for(i in 1:nrow(TunEnv)){
  print(i)
  TunEnv$site1[i]<-Tun_locs$site[(which((TunEnv$s1.xCoord[i] == Tun_locs$cellLong) & 
                                          (TunEnv$s1.yCoord[i] == Tun_locs$cellLat)))]
  TunEnv$site2[i]<-Tun_locs$site[(which((TunEnv$s2.xCoord[i] == Tun_locs$cellLong) & 
                                          (TunEnv$s2.yCoord[i] == Tun_locs$cellLat)))]
}

TunEnv$bio1Mean<-NA
TunEnv$bio12Mean<-NA
TunEnv$pStabMean<-NA
TunEnv$tStabMean<-NA
for(i in 1:nrow(TunEnv)){
  print(i)
  TunEnv$bio1Mean[i]<-mean(c(TunEnv$s1.bio1[i], TunEnv$s2.bio1[i]))
  TunEnv$bio12Mean[i]<-mean(c(TunEnv$s1.bio12[i], TunEnv$s2.bio12[i]))
  TunEnv$pStabMean[i]<-mean(c(TunEnv$s1.precipStability[i], TunEnv$s2.precipStability[i]))
  TunEnv$tStabMean[i]<-mean(c(TunEnv$s1.tempStability[i], TunEnv$s2.tempStability[i]))
}

TunEnv$geoDist<-distHaversine(matrix(c(TunEnv$s1.xCoord, TunEnv$s1.yCoord), ncol = 2),
                              matrix(c(TunEnv$s2.xCoord, TunEnv$s2.yCoord), ncol = 2))/1000

write.csv(TunEnv, "Tun/TunEnv.csv", row.names = FALSE)







All_phylogdmTab<-as.data.frame(fread("All/All_phylogdmTab.csv"))
AllEnv<-All_phylogdmTab[, columns]
colnames(AllEnv)<-c("distance", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.bio1", "s2.bio1", "s1.bio12", "s2.bio12", "s1.precipStability",
                    "s2.precipStability", "s1.tempStability", "s2.tempStability")
AllEnv$biomeabbr<-"All"
All_locs<-as.data.frame(fread("All/AllLocs.csv"))
remove(All_phylogdmTab)
# AllEnv$s1.xCoord<-round(AllEnv$s1.xCoord, 2)
# AllEnv$s2.xCoord<-round(AllEnv$s2.xCoord, 2)
# AllEnv$s1.yCoord<-round(AllEnv$s1.yCoord, 2)
# AllEnv$s2.yCoord<-round(AllEnv$s2.yCoord, 2)
# All_locs$cellLat<-round(All_locs$cellLat, 2)
# All_locs$cellLong<-round(All_locs$cellLong, 2)
# 
# AllEnv$site1<-NA
# AllEnv$site2<-NA
# for(i in 1:nrow(AllEnv)){
#   print(i)
#   AllEnv$site1[i]<-All_locs$site[(which((AllEnv$s1.xCoord[i] == All_locs$cellLong) & 
#                                           (AllEnv$s1.yCoord[i] == All_locs$cellLat)))]
#   AllEnv$site2[i]<-All_locs$site[(which((AllEnv$s2.xCoord[i] == All_locs$cellLong) & 
#                                           (AllEnv$s2.yCoord[i] == All_locs$cellLat)))]
# }
AllEnv$bio1Mean<-rowMeans(AllEnv[,c("s1.bio1", "s2.bio1")])
AllEnv$bio12Mean<-rowMeans(AllEnv[,c("s1.bio12","s2.bio12")])
AllEnv$pStabMean<-rowMeans(AllEnv[,c("s1.precipStability","s2.precipStability")])
AllEnv$tStabMean<-rowMeans(AllEnv[,c("s1.tempStability","s2.tempStability")])


# for(i in 1:nrow(AllEnv)){
#   print(i)
#   AllEnv$bio1Mean[i]<-mean(c(AllEnv$s1.bio1[i], AllEnv$s2.bio1[i]))
#   AllEnv$bio12Mean[i]<-mean(c(AllEnv$s1.bio12[i], AllEnv$s2.bio12[i]))
#   AllEnv$pStabMean[i]<-mean(c(AllEnv$s1.precipStability[i], AllEnv$s2.precipStability[i]))
#   AllEnv$tStabMean[i]<-mean(c(AllEnv$s1.tempStability[i], AllEnv$s2.tempStability[i]))
# }

AllEnv$geoDist<-distHaversine(matrix(c(AllEnv$s1.xCoord, AllEnv$s1.yCoord), ncol = 2),
                              matrix(c(AllEnv$s2.xCoord, AllEnv$s2.yCoord), ncol = 2))/1000

write.csv(AllEnv, "All/AllEnv.csv", row.names = FALSE)






## Plot logistic glms and slowps----

# Scale the predictor variables for graphical interpretability
DBFEnv$bio1Mean_C<-DBFEnv$bio1Mean/10
DBFEnv$bio12Mean_sqrt<-sqrt(DBFEnv$bio12Mean)
DBFEnv$geoDist_log<-log(DBFEnv$geoDist)
DBFMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = DBFEnv, family = "binomial")
summary(DBFMod)


# Make a df of the slopes of the each predictor variable
DBF_slowps<-as.data.frame(DBFMod[["coefficients"]][-1])
DBF_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
DBF_slowps$meanpbd<-mean(DBFEnv$distance)
DBF_slowps$mean_bio1<-mean(cbind((DBFEnv$s1.bio1/10),(DBFEnv$s2.bio1/10)), na.rm=TRUE)
DBF_slowps$mean_bio12<-mean(cbind(sqrt(DBFEnv$s1.bio12),sqrt(DBFEnv$s2.bio12)), na.rm=TRUE)
DBF_slowps$mean_pStab<-mean(cbind(DBFEnv$s1.precipStability, DBFEnv$s2.precipStability), na.rm =TRUE)
DBF_slowps$mean_tStab<-mean(cbind(DBFEnv$s1.tempStability, DBFEnv$s2.tempStability), na.rm = TRUE)
DBF_slowps$mean_geodist<-mean(DBFEnv$geoDist_log)
DBF_slowps$biomeabbr<-"DBF"
colnames(DBF_slowps)[1]<-"coefficients"
DBF_slowps$pval<-summary(DBFMod)$coefficients[2:6,4]
write.csv(DBF_slowps, "DBF/DBF_slowps.csv", row.names=FALSE )

# Check significance of individual predictor terms using likelihood ratio tests
anova(DBFMod,
      update(DBFMod, .~.-bio1Mean_C),
      test = "LRT")
anova(DBFMod,
      update(DBFMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(DBFMod,
      update(DBFMod, .~.-pStabMean),
      test = "LRT")
anova(DBFMod,
      update(DBFMod, .~.-tStabMean),
      test = "LRT")
anova(DBFMod,
      update(DBFMod, .~.-geoDist_log),
      test = "LRT")

# Check model residuals
# plot(DBFMod)

# Make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
DBFMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(DBFEnv$bio1Mean_C, na.rm =T),
                                                 max(DBFEnv$bio1Mean_C, na.rm =T),
                                                 length.out = 1000),
                                "bio12Mean_sqrt" = median(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(DBFEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(DBFEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(DBFEnv$geoDist_log, na.rm = T))

# Apply predict() function to dummy data
DBFMod_pred_bio1Mean<-cbind(DBFMod_dat_bio1Mean,
                            predict(DBFMod,DBFMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(DBFMod_pred_bio1Mean)
DBFMod_pred_bio1Mean$biomeabbr<-"DBF"
# How does the model look? plot it
#ggplot(DBFMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = DBFEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

DBFMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                                   max(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(DBFEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(DBFEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(DBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(DBFEnv$geoDist_log, na.rm = T))

# Repeat this process across other predictor variables

# apply predict() function to dummy data
DBFMod_pred_bio12Mean<-cbind(DBFMod_dat_bio12Mean,
                            predict(DBFMod,DBFMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(DBFMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(DBFMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


DBFMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(DBFEnv$pStabMean, na.rm =T),
                                                   max(DBFEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(DBFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(DBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(DBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
DBFMod_pred_pStabMean<-cbind(DBFMod_dat_pStabMean, 
                             predict(DBFMod,DBFMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(DBFMod_pred_pStabMean)
DBFMod_pred_pStabMean$biomeabbr<-"DBF"
# how does the model look? plot it
#ggplot(DBFMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



DBFMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(DBFEnv$tStabMean, na.rm =T),
                                                   max(DBFEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(DBFEnv$bio1Mean_C, na.rm =T),
                                "bio12Mean_sqrt" = median(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(DBFEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(DBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
DBFMod_pred_tStabMean<-cbind(DBFMod_dat_tStabMean,
                             predict(DBFMod,DBFMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(DBFMod_pred_tStabMean)
DBFMod_pred_tStabMean$biomeabbr<-"DBF"
# how does the model look? plot it
#ggplot(DBFMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

DBFMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(DBFEnv$geoDist_log, na.rm =T),
                                                   max(DBFEnv$geoDist_log, na.rm =T),
                                                   length.out = 1000),
                                   "bio1Mean_C" = median(DBFEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(DBFEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(DBFEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(DBFEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
DBFMod_pred_geoDistMean<-cbind(DBFMod_dat_geoDistMean,
                               predict(DBFMod,DBFMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(DBFMod_pred_geoDistMean)
DBFMod_pred_geoDistMean$biomeabbr<-"DBF"
# how does the model look? plot it
#ggplot(DBFMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(DBFEnv)




# Repeat this chunk of code 14 more times



DesEnv<-as.data.frame(fread("Des/DesEnv.csv"))
DesEnv$bio1Mean_C<-DesEnv$bio1Mean/10
DesEnv$bio12Mean_sqrt<-sqrt(DesEnv$bio12Mean)
DesEnv$geoDist_log<-log(DesEnv$geoDist)
DesMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = DesEnv, family = "binomial")
summary(DesMod)

Des_slowps<-as.data.frame(DesMod[["coefficients"]][-1])
Des_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
Des_slowps$meanpbd<-mean(DesEnv$distance)
Des_slowps$mean_bio1<-mean(cbind((DesEnv$s1.bio1/10),(DesEnv$s2.bio1/10)), na.rm=TRUE)
Des_slowps$mean_bio12<-mean(cbind(sqrt(DesEnv$s1.bio12),sqrt(DesEnv$s2.bio12)), na.rm=TRUE)
Des_slowps$mean_pStab<-mean(cbind(DesEnv$s1.precipStability, DesEnv$s2.precipStability), na.rm =TRUE)
Des_slowps$mean_tStab<-mean(cbind(DesEnv$s1.tempStability, DesEnv$s2.tempStability), na.rm = TRUE)
Des_slowps$mean_geodist<-mean(DesEnv$geoDist_log)
Des_slowps$biomeabbr<-"Des"
colnames(Des_slowps)[1]<-"coefficients"
Des_slowps$pval<-summary(DesMod)$coefficients[2:6,4]
write.csv(Des_slowps, "Des/Des_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(DesMod,
      update(DesMod, .~.-bio1Mean_C),
      test = "LRT")
anova(DesMod,
      update(DesMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(DesMod,
      update(DesMod, .~.-pStabMean),
      test = "LRT")
anova(DesMod,
      update(DesMod, .~.-tStabMean),
      test = "LRT")
anova(DesMod,
      update(DesMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(DesMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
DesMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(DesEnv$bio1Mean_C, na.rm =T),
                                                   max(DesEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(DesEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(DesEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(DesEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(DesEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
DesMod_pred_bio1Mean<-cbind(DesMod_dat_bio1Mean,
                            predict(DesMod,DesMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(DesMod_pred_bio1Mean)
DesMod_pred_bio1Mean$biomeabbr<-"Des"
# how does the model look? plot it
#ggplot(DesMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = DesEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

DesMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(DesEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(DesEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(DesEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(DesEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(DesEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(DesEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
DesMod_pred_bio12Mean<-cbind(DesMod_dat_bio12Mean,
                             predict(DesMod,DesMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(DesMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(DesMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


DesMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(DesEnv$pStabMean, na.rm =T),
                                                   max(DesEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(DesEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(DesEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(DesEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(DesEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
DesMod_pred_pStabMean<-cbind(DesMod_dat_pStabMean, 
                             predict(DesMod,DesMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(DesMod_pred_pStabMean)
DesMod_pred_pStabMean$biomeabbr<-"Des"
# how does the model look? plot it
#ggplot(DesMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



DesMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(DesEnv$tStabMean, na.rm =T),
                                                   max(DesEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(DesEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(DesEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(DesEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(DesEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
DesMod_pred_tStabMean<-cbind(DesMod_dat_tStabMean,
                             predict(DesMod,DesMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(DesMod_pred_tStabMean)
DesMod_pred_tStabMean$biomeabbr<-"Des"
# how does the model look? plot it
#ggplot(DesMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

DesMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(DesEnv$geoDist_log, na.rm =T),
                                                       max(DesEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(DesEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(DesEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(DesEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(DesEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
DesMod_pred_geoDistMean<-cbind(DesMod_dat_geoDistMean,
                               predict(DesMod,DesMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(DesMod_pred_geoDistMean)
DesMod_pred_geoDistMean$biomeabbr<-"Des"
# how does the model look? plot it
#ggplot(DesMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(DesEnv)








FloGEnv<-as.data.frame(fread("FloG/FloGEnv.csv"))
FloGEnv$bio1Mean_C<-FloGEnv$bio1Mean/10
FloGEnv$bio12Mean_sqrt<-sqrt(FloGEnv$bio12Mean)
FloGEnv$geoDist_log<-log(FloGEnv$geoDist)
FloGMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = FloGEnv, family = "binomial")
summary(FloGMod)

FloG_slowps<-as.data.frame(FloGMod[["coefficients"]][-1])
FloG_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
FloG_slowps$meanpbd<-mean(FloGEnv$distance)
FloG_slowps$mean_bio1<-mean(cbind((FloGEnv$s1.bio1/10),(FloGEnv$s2.bio1/10)), na.rm=TRUE)
FloG_slowps$mean_bio12<-mean(cbind(sqrt(FloGEnv$s1.bio12),sqrt(FloGEnv$s2.bio12)), na.rm=TRUE)
FloG_slowps$mean_pStab<-mean(cbind(FloGEnv$s1.precipStability, FloGEnv$s2.precipStability), na.rm =TRUE)
FloG_slowps$mean_tStab<-mean(cbind(FloGEnv$s1.tempStability, FloGEnv$s2.tempStability), na.rm = TRUE)
FloG_slowps$mean_geodist<-mean(FloGEnv$geoDist_log)
FloG_slowps$biomeabbr<-"FloG"
colnames(FloG_slowps)[1]<-"coefficients"
FloG_slowps$pval<-summary(FloGMod)$coefficients[2:6,4]
write.csv(FloG_slowps, "FloG/FloG_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(FloGMod,
      update(FloGMod, .~.-bio1Mean_C),
      test = "LRT")
anova(FloGMod,
      update(FloGMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(FloGMod,
      update(FloGMod, .~.-pStabMean),
      test = "LRT")
anova(FloGMod,
      update(FloGMod, .~.-tStabMean),
      test = "LRT")
anova(FloGMod,
      update(FloGMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(FloGMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
FloGMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(FloGEnv$bio1Mean_C, na.rm =T),
                                                   max(FloGEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(FloGEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(FloGEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(FloGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
FloGMod_pred_bio1Mean<-cbind(FloGMod_dat_bio1Mean,
                            predict(FloGMod,FloGMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(FloGMod_pred_bio1Mean)
FloGMod_pred_bio1Mean$biomeabbr<-"FloG"
# how does the model look? plot it
#ggplot(FloGMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = FloGEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

FloGMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(FloGEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(FloGEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(FloGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(FloGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
FloGMod_pred_bio12Mean<-cbind(FloGMod_dat_bio12Mean,
                             predict(FloGMod,FloGMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(FloGMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(FloGMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


FloGMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(FloGEnv$pStabMean, na.rm =T),
                                                   max(FloGEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(FloGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(FloGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(FloGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
FloGMod_pred_pStabMean<-cbind(FloGMod_dat_pStabMean, 
                             predict(FloGMod,FloGMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(FloGMod_pred_pStabMean)
FloGMod_pred_pStabMean$biomeabbr<-"FloG"
# how does the model look? plot it
#ggplot(FloGMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



FloGMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(FloGEnv$tStabMean, na.rm =T),
                                                   max(FloGEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(FloGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(FloGEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(FloGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
FloGMod_pred_tStabMean<-cbind(FloGMod_dat_tStabMean,
                             predict(FloGMod,FloGMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(FloGMod_pred_tStabMean)
FloGMod_pred_tStabMean$biomeabbr<-"FloG"
# how does the model look? plot it
#ggplot(FloGMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

FloGMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(FloGEnv$geoDist_log, na.rm =T),
                                                       max(FloGEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(FloGEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(FloGEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(FloGEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(FloGEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
FloGMod_pred_geoDistMean<-cbind(FloGMod_dat_geoDistMean,
                               predict(FloGMod,FloGMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(FloGMod_pred_geoDistMean)
FloGMod_pred_geoDistMean$biomeabbr<-"FloG"
# how does the model look? plot it
#ggplot(FloGMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(FloGEnv)









MangEnv<-as.data.frame(fread("Mang/MangEnv.csv"))
MangEnv$bio1Mean_C<-MangEnv$bio1Mean/10
MangEnv$bio12Mean_sqrt<-sqrt(MangEnv$bio12Mean)
MangEnv$geoDist_log<-log(MangEnv$geoDist)
MangMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = MangEnv, family = "binomial")
summary(MangMod)

Mang_slowps<-as.data.frame(MangMod[["coefficients"]][-1])
Mang_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
Mang_slowps$meanpbd<-mean(MangEnv$distance)
Mang_slowps$mean_bio1<-mean(cbind((MangEnv$s1.bio1/10),(MangEnv$s2.bio1/10)), na.rm=TRUE)
Mang_slowps$mean_bio12<-mean(cbind(sqrt(MangEnv$s1.bio12),sqrt(MangEnv$s2.bio12)), na.rm=TRUE)
Mang_slowps$mean_pStab<-mean(cbind(MangEnv$s1.precipStability, MangEnv$s2.precipStability), na.rm =TRUE)
Mang_slowps$mean_tStab<-mean(cbind(MangEnv$s1.tempStability, MangEnv$s2.tempStability), na.rm = TRUE)
Mang_slowps$mean_geodist<-mean(MangEnv$geoDist_log)
Mang_slowps$biomeabbr<-"Mang"
colnames(Mang_slowps)[1]<-"coefficients"
Mang_slowps$pval<-summary(MangMod)$coefficients[2:6,4]
write.csv(Mang_slowps, "Mang/Mang_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(MangMod,
      update(MangMod, .~.-bio1Mean_C),
      test = "LRT")
anova(MangMod,
      update(MangMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(MangMod,
      update(MangMod, .~.-pStabMean),
      test = "LRT")
anova(MangMod,
      update(MangMod, .~.-tStabMean),
      test = "LRT")
anova(MangMod,
      update(MangMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(MangMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
MangMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(MangEnv$bio1Mean_C, na.rm =T),
                                                   max(MangEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(MangEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(MangEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(MangEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(MangEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MangMod_pred_bio1Mean<-cbind(MangMod_dat_bio1Mean,
                            predict(MangMod,MangMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(MangMod_pred_bio1Mean)
MangMod_pred_bio1Mean$biomeabbr<-"Mang"
# how does the model look? plot it
#ggplot(MangMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = MangEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

MangMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(MangEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(MangEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(MangEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(MangEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(MangEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MangEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MangMod_pred_bio12Mean<-cbind(MangMod_dat_bio12Mean,
                             predict(MangMod,MangMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(MangMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(MangMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


MangMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(MangEnv$pStabMean, na.rm =T),
                                                   max(MangEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MangEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MangEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(MangEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MangEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MangMod_pred_pStabMean<-cbind(MangMod_dat_pStabMean, 
                             predict(MangMod,MangMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(MangMod_pred_pStabMean)
MangMod_pred_pStabMean$biomeabbr<-"Mang"
# how does the model look? plot it
#ggplot(MangMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



MangMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(MangEnv$tStabMean, na.rm =T),
                                                   max(MangEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MangEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MangEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(MangEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(MangEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MangMod_pred_tStabMean<-cbind(MangMod_dat_tStabMean,
                             predict(MangMod,MangMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(MangMod_pred_tStabMean)
MangMod_pred_tStabMean$biomeabbr<-"Mang"
# how does the model look? plot it
#ggplot(MangMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

MangMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(MangEnv$geoDist_log, na.rm =T),
                                                       max(MangEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(MangEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(MangEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(MangEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(MangEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
MangMod_pred_geoDistMean<-cbind(MangMod_dat_geoDistMean,
                               predict(MangMod,MangMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(MangMod_pred_geoDistMean)
MangMod_pred_geoDistMean$biomeabbr<-"Mang"
# how does the model look? plot it
#ggplot(MangMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(MangEnv)









MBFEnv<-as.data.frame(fread("MBF/MBFEnv.csv"))
MBFEnv$bio1Mean_C<-MBFEnv$bio1Mean/10
MBFEnv$bio12Mean_sqrt<-sqrt(MBFEnv$bio12Mean)
MBFEnv$geoDist_log<-log(MBFEnv$geoDist)
MBFMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = MBFEnv, family = "binomial")
summary(MBFMod)

MBF_slowps<-as.data.frame(MBFMod[["coefficients"]][-1])
MBF_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
MBF_slowps$meanpbd<-mean(MBFEnv$distance)
MBF_slowps$mean_bio1<-mean(cbind((MBFEnv$s1.bio1/10),(MBFEnv$s2.bio1/10)), na.rm=TRUE)
MBF_slowps$mean_bio12<-mean(cbind(sqrt(MBFEnv$s1.bio12),sqrt(MBFEnv$s2.bio12)), na.rm=TRUE)
MBF_slowps$mean_pStab<-mean(cbind(MBFEnv$s1.precipStability, MBFEnv$s2.precipStability), na.rm =TRUE)
MBF_slowps$mean_tStab<-mean(cbind(MBFEnv$s1.tempStability, MBFEnv$s2.tempStability), na.rm = TRUE)
MBF_slowps$mean_geodist<-mean(MBFEnv$geoDist_log)
MBF_slowps$biomeabbr<-"MBF"
colnames(MBF_slowps)[1]<-"coefficients"
MBF_slowps$pval<-summary(MBFMod)$coefficients[2:6,4]
write.csv(MBF_slowps, "MBF/MBF_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(MBFMod,
      update(MBFMod, .~.-bio1Mean_C),
      test = "LRT")
anova(MBFMod,
      update(MBFMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(MBFMod,
      update(MBFMod, .~.-pStabMean),
      test = "LRT")
anova(MBFMod,
      update(MBFMod, .~.-tStabMean),
      test = "LRT")
anova(MBFMod,
      update(MBFMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(MBFMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
MBFMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(MBFEnv$bio1Mean_C, na.rm =T),
                                                   max(MBFEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(MBFEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(MBFEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(MBFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MBFMod_pred_bio1Mean<-cbind(MBFMod_dat_bio1Mean,
                            predict(MBFMod,MBFMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(MBFMod_pred_bio1Mean)
MBFMod_pred_bio1Mean$biomeabbr<-"MBF"
# how does the model look? plot it
#ggplot(MBFMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = MBFEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

MBFMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(MBFEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(MBFEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(MBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MBFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MBFMod_pred_bio12Mean<-cbind(MBFMod_dat_bio12Mean,
                             predict(MBFMod,MBFMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(MBFMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(MBFMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


MBFMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(MBFEnv$pStabMean, na.rm =T),
                                                   max(MBFEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MBFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(MBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MBFMod_pred_pStabMean<-cbind(MBFMod_dat_pStabMean, 
                             predict(MBFMod,MBFMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(MBFMod_pred_pStabMean)
MBFMod_pred_pStabMean$biomeabbr<-"MBF"
# how does the model look? plot it
#ggplot(MBFMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



MBFMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(MBFEnv$tStabMean, na.rm =T),
                                                   max(MBFEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MBFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(MBFEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(MBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MBFMod_pred_tStabMean<-cbind(MBFMod_dat_tStabMean,
                             predict(MBFMod,MBFMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(MBFMod_pred_tStabMean)
MBFMod_pred_tStabMean$biomeabbr<-"MBF"
# how does the model look? plot it
#ggplot(MBFMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

MBFMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(MBFEnv$geoDist_log, na.rm =T),
                                                       max(MBFEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(MBFEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(MBFEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(MBFEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(MBFEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
MBFMod_pred_geoDistMean<-cbind(MBFMod_dat_geoDistMean,
                               predict(MBFMod,MBFMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(MBFMod_pred_geoDistMean)
MBFMod_pred_geoDistMean$biomeabbr<-"MBF"
# how does the model look? plot it
#ggplot(MBFMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(MBFEnv)









MeditEnv<-as.data.frame(fread("Medit/MeditEnv.csv"))
MeditEnv$bio1Mean_C<-MeditEnv$bio1Mean/10
MeditEnv$bio12Mean_sqrt<-sqrt(MeditEnv$bio12Mean)
MeditEnv$geoDist_log<-log(MeditEnv$geoDist)
MeditMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = MeditEnv, family = "binomial")
summary(MeditMod)

Medit_slowps<-as.data.frame(MeditMod[["coefficients"]][-1])
Medit_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
Medit_slowps$meanpbd<-mean(MeditEnv$distance)
Medit_slowps$mean_bio1<-mean(cbind((MeditEnv$s1.bio1/10),(MeditEnv$s2.bio1/10)), na.rm=TRUE)
Medit_slowps$mean_bio12<-mean(cbind(sqrt(MeditEnv$s1.bio12),sqrt(MeditEnv$s2.bio12)), na.rm=TRUE)
Medit_slowps$mean_pStab<-mean(cbind(MeditEnv$s1.precipStability, MeditEnv$s2.precipStability), na.rm =TRUE)
Medit_slowps$mean_tStab<-mean(cbind(MeditEnv$s1.tempStability, MeditEnv$s2.tempStability), na.rm = TRUE)
Medit_slowps$mean_geodist<-mean(MeditEnv$geoDist_log)
Medit_slowps$biomeabbr<-"Medit"
colnames(Medit_slowps)[1]<-"coefficients"
Medit_slowps$pval<-summary(MeditMod)$coefficients[2:6,4]
write.csv(Medit_slowps, "Medit/Medit_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(MeditMod,
      update(MeditMod, .~.-bio1Mean_C),
      test = "LRT")
anova(MeditMod,
      update(MeditMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(MeditMod,
      update(MeditMod, .~.-pStabMean),
      test = "LRT")
anova(MeditMod,
      update(MeditMod, .~.-tStabMean),
      test = "LRT")
anova(MeditMod,
      update(MeditMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(MeditMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
MeditMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(MeditEnv$bio1Mean_C, na.rm =T),
                                                   max(MeditEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(MeditEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(MeditEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(MeditEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MeditMod_pred_bio1Mean<-cbind(MeditMod_dat_bio1Mean,
                            predict(MeditMod,MeditMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(MeditMod_pred_bio1Mean)
MeditMod_pred_bio1Mean$biomeabbr<-"Medit"
# how does the model look? plot it
#ggplot(MeditMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = MeditEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

MeditMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(MeditEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(MeditEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(MeditEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MeditEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MeditMod_pred_bio12Mean<-cbind(MeditMod_dat_bio12Mean,
                             predict(MeditMod,MeditMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(MeditMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(MeditMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


MeditMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(MeditEnv$pStabMean, na.rm =T),
                                                   max(MeditEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MeditEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(MeditEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MeditEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MeditMod_pred_pStabMean<-cbind(MeditMod_dat_pStabMean, 
                             predict(MeditMod,MeditMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(MeditMod_pred_pStabMean)
MeditMod_pred_pStabMean$biomeabbr<-"Medit"
# how does the model look? plot it
#ggplot(MeditMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



MeditMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(MeditEnv$tStabMean, na.rm =T),
                                                   max(MeditEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MeditEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(MeditEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(MeditEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MeditMod_pred_tStabMean<-cbind(MeditMod_dat_tStabMean,
                             predict(MeditMod,MeditMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(MeditMod_pred_tStabMean)
MeditMod_pred_tStabMean$biomeabbr<-"Medit"
# how does the model look? plot it
#ggplot(MeditMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

MeditMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(MeditEnv$geoDist_log, na.rm =T),
                                                       max(MeditEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(MeditEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(MeditEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(MeditEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(MeditEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
MeditMod_pred_geoDistMean<-cbind(MeditMod_dat_geoDistMean,
                               predict(MeditMod,MeditMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(MeditMod_pred_geoDistMean)
MeditMod_pred_geoDistMean$biomeabbr<-"Medit"
# how does the model look? plot it
#ggplot(MeditMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(MeditEnv)









MonGEnv<-as.data.frame(fread("MonG/MonGEnv.csv"))
MonGEnv$bio1Mean_C<-MonGEnv$bio1Mean/10
MonGEnv$bio12Mean_sqrt<-sqrt(MonGEnv$bio12Mean)
MonGEnv$geoDist_log<-log(MonGEnv$geoDist)
MonGMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = MonGEnv, family = "binomial")
summary(MonGMod)

MonG_slowps<-as.data.frame(MonGMod[["coefficients"]][-1])
MonG_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
MonG_slowps$meanpbd<-mean(MonGEnv$distance)
MonG_slowps$mean_bio1<-mean(cbind((MonGEnv$s1.bio1/10),(MonGEnv$s2.bio1/10)), na.rm=TRUE)
MonG_slowps$mean_bio12<-mean(cbind(sqrt(MonGEnv$s1.bio12),sqrt(MonGEnv$s2.bio12)), na.rm=TRUE)
MonG_slowps$mean_pStab<-mean(cbind(MonGEnv$s1.precipStability, MonGEnv$s2.precipStability), na.rm =TRUE)
MonG_slowps$mean_tStab<-mean(cbind(MonGEnv$s1.tempStability, MonGEnv$s2.tempStability), na.rm = TRUE)
MonG_slowps$mean_geodist<-mean(MonGEnv$geoDist_log)
MonG_slowps$biomeabbr<-"MonG"
colnames(MonG_slowps)[1]<-"coefficients"
MonG_slowps$pval<-summary(MonGMod)$coefficients[2:6,4]
write.csv(MonG_slowps, "MonG/MonG_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(MonGMod,
      update(MonGMod, .~.-bio1Mean_C),
      test = "LRT")
anova(MonGMod,
      update(MonGMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(MonGMod,
      update(MonGMod, .~.-pStabMean),
      test = "LRT")
anova(MonGMod,
      update(MonGMod, .~.-tStabMean),
      test = "LRT")
anova(MonGMod,
      update(MonGMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(MonGMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
MonGMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(MonGEnv$bio1Mean_C, na.rm =T),
                                                   max(MonGEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(MonGEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(MonGEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(MonGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MonGMod_pred_bio1Mean<-cbind(MonGMod_dat_bio1Mean,
                            predict(MonGMod,MonGMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(MonGMod_pred_bio1Mean)
MonGMod_pred_bio1Mean$biomeabbr<-"MonG"
# how does the model look? plot it
#ggplot(MonGMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = MonGEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

MonGMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(MonGEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(MonGEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(MonGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MonGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
MonGMod_pred_bio12Mean<-cbind(MonGMod_dat_bio12Mean,
                             predict(MonGMod,MonGMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(MonGMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(MonGMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


MonGMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(MonGEnv$pStabMean, na.rm =T),
                                                   max(MonGEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MonGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(MonGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(MonGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MonGMod_pred_pStabMean<-cbind(MonGMod_dat_pStabMean, 
                             predict(MonGMod,MonGMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(MonGMod_pred_pStabMean)
MonGMod_pred_pStabMean$biomeabbr<-"MonG"
# how does the model look? plot it
#ggplot(MonGMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



MonGMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(MonGEnv$tStabMean, na.rm =T),
                                                   max(MonGEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(MonGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(MonGEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(MonGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
MonGMod_pred_tStabMean<-cbind(MonGMod_dat_tStabMean,
                             predict(MonGMod,MonGMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(MonGMod_pred_tStabMean)
MonGMod_pred_tStabMean$biomeabbr<-"MonG"
# how does the model look? plot it
#ggplot(MonGMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

MonGMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(MonGEnv$geoDist_log, na.rm =T),
                                                       max(MonGEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(MonGEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(MonGEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(MonGEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(MonGEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
MonGMod_pred_geoDistMean<-cbind(MonGMod_dat_geoDistMean,
                               predict(MonGMod,MonGMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(MonGMod_pred_geoDistMean)
MonGMod_pred_geoDistMean$biomeabbr<-"MonG"
# how does the model look? plot it
#ggplot(MonGMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(MonGEnv)








TaiEnv<-as.data.frame(fread("Tai/TaiEnv.csv"))
TaiEnv$bio1Mean_C<-TaiEnv$bio1Mean/10
TaiEnv$bio12Mean_sqrt<-sqrt(TaiEnv$bio12Mean)
TaiEnv$geoDist_log<-log(TaiEnv$geoDist)
TaiMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TaiEnv, family = "binomial")
summary(TaiMod)

Tai_slowps<-as.data.frame(TaiMod[["coefficients"]][-1])
Tai_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
Tai_slowps$meanpbd<-mean(TaiEnv$distance)
Tai_slowps$mean_bio1<-mean(cbind((TaiEnv$s1.bio1/10),(TaiEnv$s2.bio1/10)), na.rm=TRUE)
Tai_slowps$mean_bio12<-mean(cbind(sqrt(TaiEnv$s1.bio12),sqrt(TaiEnv$s2.bio12)), na.rm=TRUE)
Tai_slowps$mean_pStab<-mean(cbind(TaiEnv$s1.precipStability, TaiEnv$s2.precipStability), na.rm =TRUE)
Tai_slowps$mean_tStab<-mean(cbind(TaiEnv$s1.tempStability, TaiEnv$s2.tempStability), na.rm = TRUE)
Tai_slowps$mean_geodist<-mean(TaiEnv$geoDist_log)
Tai_slowps$biomeabbr<-"Tai"
colnames(Tai_slowps)[1]<-"coefficients"
Tai_slowps$pval<-summary(TaiMod)$coefficients[2:6,4]
write.csv(Tai_slowps, "Tai/Tai_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TaiMod,
      update(TaiMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TaiMod,
      update(TaiMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TaiMod,
      update(TaiMod, .~.-pStabMean),
      test = "LRT")
anova(TaiMod,
      update(TaiMod, .~.-tStabMean),
      test = "LRT")
anova(TaiMod,
      update(TaiMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TaiMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TaiMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TaiEnv$bio1Mean_C, na.rm =T),
                                                   max(TaiEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TaiEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TaiEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TaiEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TaiMod_pred_bio1Mean<-cbind(TaiMod_dat_bio1Mean,
                            predict(TaiMod,TaiMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TaiMod_pred_bio1Mean)
TaiMod_pred_bio1Mean$biomeabbr<-"Tai"
# how does the model look? plot it
#ggplot(TaiMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TaiEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TaiMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TaiEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TaiEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TaiEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TaiEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TaiMod_pred_bio12Mean<-cbind(TaiMod_dat_bio12Mean,
                             predict(TaiMod,TaiMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TaiMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TaiMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TaiMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TaiEnv$pStabMean, na.rm =T),
                                                   max(TaiEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TaiEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TaiEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TaiEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TaiMod_pred_pStabMean<-cbind(TaiMod_dat_pStabMean, 
                             predict(TaiMod,TaiMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TaiMod_pred_pStabMean)
TaiMod_pred_pStabMean$biomeabbr<-"Tai"
# how does the model look? plot it
#ggplot(TaiMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TaiMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TaiEnv$tStabMean, na.rm =T),
                                                   max(TaiEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TaiEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TaiEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TaiEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TaiMod_pred_tStabMean<-cbind(TaiMod_dat_tStabMean,
                             predict(TaiMod,TaiMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TaiMod_pred_tStabMean)
TaiMod_pred_tStabMean$biomeabbr<-"Tai"
# how does the model look? plot it
#ggplot(TaiMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TaiMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TaiEnv$geoDist_log, na.rm =T),
                                                       max(TaiEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TaiEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TaiEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TaiEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TaiEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TaiMod_pred_geoDistMean<-cbind(TaiMod_dat_geoDistMean,
                               predict(TaiMod,TaiMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TaiMod_pred_geoDistMean)
TaiMod_pred_geoDistMean$biomeabbr<-"Tai"
# how does the model look? plot it
#ggplot(TaiMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TaiEnv)









TemBFEnv<-as.data.frame(fread("TemBF/TemBFEnv.csv"))
TemBFEnv$bio1Mean_C<-TemBFEnv$bio1Mean/10
TemBFEnv$bio12Mean_sqrt<-sqrt(TemBFEnv$bio12Mean)
TemBFEnv$geoDist_log<-log(TemBFEnv$geoDist)
TemBFMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TemBFEnv, family = "binomial")
summary(TemBFMod)

TemBF_slowps<-as.data.frame(TemBFMod[["coefficients"]][-1])
TemBF_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
TemBF_slowps$meanpbd<-mean(TemBFEnv$distance)
TemBF_slowps$mean_bio1<-mean(cbind((TemBFEnv$s1.bio1/10),(TemBFEnv$s2.bio1/10)), na.rm=TRUE)
TemBF_slowps$mean_bio12<-mean(cbind(sqrt(TemBFEnv$s1.bio12),sqrt(TemBFEnv$s2.bio12)), na.rm=TRUE)
TemBF_slowps$mean_pStab<-mean(cbind(TemBFEnv$s1.precipStability, TemBFEnv$s2.precipStability), na.rm =TRUE)
TemBF_slowps$mean_tStab<-mean(cbind(TemBFEnv$s1.tempStability, TemBFEnv$s2.tempStability), na.rm = TRUE)
TemBF_slowps$mean_geodist<-mean(TemBFEnv$geoDist_log)
TemBF_slowps$biomeabbr<-"TemBF"
colnames(TemBF_slowps)[1]<-"coefficients"
TemBF_slowps$pval<-summary(TemBFMod)$coefficients[2:6,4]
write.csv(TemBF_slowps, "TemBF/TemBF_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TemBFMod,
      update(TemBFMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TemBFMod,
      update(TemBFMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TemBFMod,
      update(TemBFMod, .~.-pStabMean),
      test = "LRT")
anova(TemBFMod,
      update(TemBFMod, .~.-tStabMean),
      test = "LRT")
anova(TemBFMod,
      update(TemBFMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TemBFMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TemBFMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TemBFEnv$bio1Mean_C, na.rm =T),
                                                   max(TemBFEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TemBFEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TemBFEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TemBFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemBFMod_pred_bio1Mean<-cbind(TemBFMod_dat_bio1Mean,
                            predict(TemBFMod,TemBFMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TemBFMod_pred_bio1Mean)
TemBFMod_pred_bio1Mean$biomeabbr<-"TemBF"
# how does the model look? plot it
#ggplot(TemBFMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TemBFEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TemBFMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TemBFEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TemBFEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TemBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemBFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemBFMod_pred_bio12Mean<-cbind(TemBFMod_dat_bio12Mean,
                             predict(TemBFMod,TemBFMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TemBFMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TemBFMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TemBFMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TemBFEnv$pStabMean, na.rm =T),
                                                   max(TemBFEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemBFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TemBFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemBFMod_pred_pStabMean<-cbind(TemBFMod_dat_pStabMean, 
                             predict(TemBFMod,TemBFMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TemBFMod_pred_pStabMean)
TemBFMod_pred_pStabMean$biomeabbr<-"TemBF"
# how does the model look? plot it
#ggplot(TemBFMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TemBFMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TemBFEnv$tStabMean, na.rm =T),
                                                   max(TemBFEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemBFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TemBFEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TemBFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemBFMod_pred_tStabMean<-cbind(TemBFMod_dat_tStabMean,
                             predict(TemBFMod,TemBFMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TemBFMod_pred_tStabMean)
TemBFMod_pred_tStabMean$biomeabbr<-"TemBF"
# how does the model look? plot it
#ggplot(TemBFMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TemBFMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TemBFEnv$geoDist_log, na.rm =T),
                                                       max(TemBFEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TemBFEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TemBFEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TemBFEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TemBFEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TemBFMod_pred_geoDistMean<-cbind(TemBFMod_dat_geoDistMean,
                               predict(TemBFMod,TemBFMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TemBFMod_pred_geoDistMean)
TemBFMod_pred_geoDistMean$biomeabbr<-"TemBF"
# how does the model look? plot it
#ggplot(TemBFMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TemBFEnv)







TemCFEnv<-as.data.frame(fread("TemCF/TemCFEnv.csv"))
TemCFEnv$bio1Mean_C<-TemCFEnv$bio1Mean/10
TemCFEnv$bio12Mean_sqrt<-sqrt(TemCFEnv$bio12Mean)
TemCFEnv$geoDist_log<-log(TemCFEnv$geoDist)
TemCFMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TemCFEnv, family = "binomial")
summary(TemCFMod)

TemCF_slowps<-as.data.frame(TemCFMod[["coefficients"]][-1])
TemCF_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
TemCF_slowps$meanpbd<-mean(TemCFEnv$distance)
TemCF_slowps$mean_bio1<-mean(cbind((TemCFEnv$s1.bio1/10),(TemCFEnv$s2.bio1/10)), na.rm=TRUE)
TemCF_slowps$mean_bio12<-mean(cbind(sqrt(TemCFEnv$s1.bio12),sqrt(TemCFEnv$s2.bio12)), na.rm=TRUE)
TemCF_slowps$mean_pStab<-mean(cbind(TemCFEnv$s1.precipStability, TemCFEnv$s2.precipStability), na.rm =TRUE)
TemCF_slowps$mean_tStab<-mean(cbind(TemCFEnv$s1.tempStability, TemCFEnv$s2.tempStability), na.rm = TRUE)
TemCF_slowps$mean_geodist<-mean(TemCFEnv$geoDist_log)
TemCF_slowps$biomeabbr<-"TemCF"
colnames(TemCF_slowps)[1]<-"coefficients"
TemCF_slowps$pval<-summary(TemCFMod)$coefficients[2:6,4]
write.csv(TemCF_slowps, "TemCF/TemCF_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TemCFMod,
      update(TemCFMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TemCFMod,
      update(TemCFMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TemCFMod,
      update(TemCFMod, .~.-pStabMean),
      test = "LRT")
anova(TemCFMod,
      update(TemCFMod, .~.-tStabMean),
      test = "LRT")
anova(TemCFMod,
      update(TemCFMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TemCFMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TemCFMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TemCFEnv$bio1Mean_C, na.rm =T),
                                                   max(TemCFEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TemCFEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TemCFEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TemCFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemCFMod_pred_bio1Mean<-cbind(TemCFMod_dat_bio1Mean,
                            predict(TemCFMod,TemCFMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TemCFMod_pred_bio1Mean)
TemCFMod_pred_bio1Mean$biomeabbr<-"TemCF"
# how does the model look? plot it
#ggplot(TemCFMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TemCFEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TemCFMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TemCFEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TemCFEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TemCFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemCFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemCFMod_pred_bio12Mean<-cbind(TemCFMod_dat_bio12Mean,
                             predict(TemCFMod,TemCFMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TemCFMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TemCFMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TemCFMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TemCFEnv$pStabMean, na.rm =T),
                                                   max(TemCFEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemCFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TemCFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemCFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemCFMod_pred_pStabMean<-cbind(TemCFMod_dat_pStabMean, 
                             predict(TemCFMod,TemCFMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TemCFMod_pred_pStabMean)
TemCFMod_pred_pStabMean$biomeabbr<-"TemCF"
# how does the model look? plot it
#ggplot(TemCFMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TemCFMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TemCFEnv$tStabMean, na.rm =T),
                                                   max(TemCFEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemCFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TemCFEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TemCFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemCFMod_pred_tStabMean<-cbind(TemCFMod_dat_tStabMean,
                             predict(TemCFMod,TemCFMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TemCFMod_pred_tStabMean)
TemCFMod_pred_tStabMean$biomeabbr<-"TemCF"
# how does the model look? plot it
#ggplot(TemCFMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TemCFMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TemCFEnv$geoDist_log, na.rm =T),
                                                       max(TemCFEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TemCFEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TemCFEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TemCFEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TemCFEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TemCFMod_pred_geoDistMean<-cbind(TemCFMod_dat_geoDistMean,
                               predict(TemCFMod,TemCFMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TemCFMod_pred_geoDistMean)
TemCFMod_pred_geoDistMean$biomeabbr<-"TemCF"
# how does the model look? plot it
#ggplot(TemCFMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TemCFEnv)






TemGEnv<-as.data.frame(fread("TemG/TemGEnv.csv"))
TemGEnv$bio1Mean_C<-TemGEnv$bio1Mean/10
TemGEnv$bio12Mean_sqrt<-sqrt(TemGEnv$bio12Mean)
TemGEnv$geoDist_log<-log(TemGEnv$geoDist)
TemGMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TemGEnv, family = "binomial")
summary(TemGMod)

TemG_slowps<-as.data.frame(TemGMod[["coefficients"]][-1])
TemG_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
TemG_slowps$meanpbd<-mean(TemGEnv$distance)
TemG_slowps$mean_bio1<-mean(cbind((TemGEnv$s1.bio1/10),(TemGEnv$s2.bio1/10)), na.rm=TRUE)
TemG_slowps$mean_bio12<-mean(cbind(sqrt(TemGEnv$s1.bio12),sqrt(TemGEnv$s2.bio12)), na.rm=TRUE)
TemG_slowps$mean_pStab<-mean(cbind(TemGEnv$s1.precipStability, TemGEnv$s2.precipStability), na.rm =TRUE)
TemG_slowps$mean_tStab<-mean(cbind(TemGEnv$s1.tempStability, TemGEnv$s2.tempStability), na.rm = TRUE)
TemG_slowps$mean_geodist<-mean(TemGEnv$geoDist_log)
TemG_slowps$biomeabbr<-"TemG"
colnames(TemG_slowps)[1]<-"coefficients"
TemG_slowps$pval<-summary(TemGMod)$coefficients[2:6,4]
write.csv(TemG_slowps, "TemG/TemG_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TemGMod,
      update(TemGMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TemGMod,
      update(TemGMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TemGMod,
      update(TemGMod, .~.-pStabMean),
      test = "LRT")
anova(TemGMod,
      update(TemGMod, .~.-tStabMean),
      test = "LRT")
anova(TemGMod,
      update(TemGMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TemGMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TemGMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TemGEnv$bio1Mean_C, na.rm =T),
                                                   max(TemGEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TemGEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TemGEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TemGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemGMod_pred_bio1Mean<-cbind(TemGMod_dat_bio1Mean,
                            predict(TemGMod,TemGMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TemGMod_pred_bio1Mean)
TemGMod_pred_bio1Mean$biomeabbr<-"TemG"
# how does the model look? plot it
#ggplot(TemGMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TemGEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TemGMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TemGEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TemGEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TemGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TemGMod_pred_bio12Mean<-cbind(TemGMod_dat_bio12Mean,
                             predict(TemGMod,TemGMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TemGMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TemGMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TemGMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TemGEnv$pStabMean, na.rm =T),
                                                   max(TemGEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TemGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TemGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemGMod_pred_pStabMean<-cbind(TemGMod_dat_pStabMean, 
                             predict(TemGMod,TemGMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TemGMod_pred_pStabMean)
TemGMod_pred_pStabMean$biomeabbr<-"TemG"
# how does the model look? plot it
#ggplot(TemGMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TemGMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TemGEnv$tStabMean, na.rm =T),
                                                   max(TemGEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TemGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TemGEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TemGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TemGMod_pred_tStabMean<-cbind(TemGMod_dat_tStabMean,
                             predict(TemGMod,TemGMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TemGMod_pred_tStabMean)
TemGMod_pred_tStabMean$biomeabbr<-"TemG"
# how does the model look? plot it
#ggplot(TemGMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TemGMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TemGEnv$geoDist_log, na.rm =T),
                                                       max(TemGEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TemGEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TemGEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TemGEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TemGEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TemGMod_pred_geoDistMean<-cbind(TemGMod_dat_geoDistMean,
                               predict(TemGMod,TemGMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TemGMod_pred_geoDistMean)
TemGMod_pred_geoDistMean$biomeabbr<-"TemG"
# how does the model look? plot it
#ggplot(TemGMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TemGEnv)








TroCFEnv<-as.data.frame(fread("TroCF/TroCFEnv.csv"))
TroCFEnv$bio1Mean_C<-TroCFEnv$bio1Mean/10
TroCFEnv$bio12Mean_sqrt<-sqrt(TroCFEnv$bio12Mean)
TroCFEnv$geoDist_log<-log(TroCFEnv$geoDist)
TroCFMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TroCFEnv, family = "binomial")
summary(TroCFMod)

TroCF_slowps<-as.data.frame(TroCFMod[["coefficients"]][-1])
TroCF_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
TroCF_slowps$meanpbd<-mean(TroCFEnv$distance)
TroCF_slowps$mean_bio1<-mean(cbind((TroCFEnv$s1.bio1/10),(TroCFEnv$s2.bio1/10)), na.rm=TRUE)
TroCF_slowps$mean_bio12<-mean(cbind(sqrt(TroCFEnv$s1.bio12),sqrt(TroCFEnv$s2.bio12)), na.rm=TRUE)
TroCF_slowps$mean_pStab<-mean(cbind(TroCFEnv$s1.precipStability, TroCFEnv$s2.precipStability), na.rm =TRUE)
TroCF_slowps$mean_tStab<-mean(cbind(TroCFEnv$s1.tempStability, TroCFEnv$s2.tempStability), na.rm = TRUE)
TroCF_slowps$mean_geodist<-mean(TroCFEnv$geoDist_log)
TroCF_slowps$biomeabbr<-"TroCF"
colnames(TroCF_slowps)[1]<-"coefficients"
TroCF_slowps$pval<-summary(TroCFMod)$coefficients[2:6,4]
write.csv(TroCF_slowps, "TroCF/TroCF_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TroCFMod,
      update(TroCFMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TroCFMod,
      update(TroCFMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TroCFMod,
      update(TroCFMod, .~.-pStabMean),
      test = "LRT")
anova(TroCFMod,
      update(TroCFMod, .~.-tStabMean),
      test = "LRT")
anova(TroCFMod,
      update(TroCFMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TroCFMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TroCFMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TroCFEnv$bio1Mean_C, na.rm =T),
                                                   max(TroCFEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TroCFEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TroCFEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TroCFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TroCFMod_pred_bio1Mean<-cbind(TroCFMod_dat_bio1Mean,
                            predict(TroCFMod,TroCFMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TroCFMod_pred_bio1Mean)
TroCFMod_pred_bio1Mean$biomeabbr<-"TroCF"
# how does the model look? plot it
#ggplot(TroCFMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TroCFEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TroCFMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TroCFEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TroCFEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TroCFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TroCFEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TroCFMod_pred_bio12Mean<-cbind(TroCFMod_dat_bio12Mean,
                             predict(TroCFMod,TroCFMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TroCFMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TroCFMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TroCFMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TroCFEnv$pStabMean, na.rm =T),
                                                   max(TroCFEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TroCFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TroCFEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TroCFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TroCFMod_pred_pStabMean<-cbind(TroCFMod_dat_pStabMean, 
                             predict(TroCFMod,TroCFMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TroCFMod_pred_pStabMean)
TroCFMod_pred_pStabMean$biomeabbr<-"TroCF"
# how does the model look? plot it
#ggplot(TroCFMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TroCFMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TroCFEnv$tStabMean, na.rm =T),
                                                   max(TroCFEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TroCFEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TroCFEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TroCFEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TroCFMod_pred_tStabMean<-cbind(TroCFMod_dat_tStabMean,
                             predict(TroCFMod,TroCFMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TroCFMod_pred_tStabMean)
TroCFMod_pred_tStabMean$biomeabbr<-"TroCF"
# how does the model look? plot it
#ggplot(TroCFMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TroCFMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TroCFEnv$geoDist_log, na.rm =T),
                                                       max(TroCFEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TroCFEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TroCFEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TroCFEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TroCFEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TroCFMod_pred_geoDistMean<-cbind(TroCFMod_dat_geoDistMean,
                               predict(TroCFMod,TroCFMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TroCFMod_pred_geoDistMean)
TroCFMod_pred_geoDistMean$biomeabbr<-"TroCF"
# how does the model look? plot it
#ggplot(TroCFMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TroCFEnv)







TroGEnv<-as.data.frame(fread("TroG/TroGEnv.csv"))
TroGEnv$bio1Mean_C<-TroGEnv$bio1Mean/10
TroGEnv$bio12Mean_sqrt<-sqrt(TroGEnv$bio12Mean)
TroGEnv$geoDist_log<-log(TroGEnv$geoDist)
TroGMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TroGEnv, family = "binomial")
summary(TroGMod)

TroG_slowps<-as.data.frame(TroGMod[["coefficients"]][-1])
TroG_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
TroG_slowps$meanpbd<-mean(TroGEnv$distance)
TroG_slowps$mean_bio1<-mean(cbind((TroGEnv$s1.bio1/10),(TroGEnv$s2.bio1/10)), na.rm=TRUE)
TroG_slowps$mean_bio12<-mean(cbind(sqrt(TroGEnv$s1.bio12),sqrt(TroGEnv$s2.bio12)), na.rm=TRUE)
TroG_slowps$mean_pStab<-mean(cbind(TroGEnv$s1.precipStability, TroGEnv$s2.precipStability), na.rm =TRUE)
TroG_slowps$mean_tStab<-mean(cbind(TroGEnv$s1.tempStability, TroGEnv$s2.tempStability), na.rm = TRUE)
TroG_slowps$mean_geodist<-mean(TroGEnv$geoDist_log)
TroG_slowps$biomeabbr<-"TroG"
colnames(TroG_slowps)[1]<-"coefficients"
TroG_slowps$pval<-summary(TroGMod)$coefficients[2:6,4]
write.csv(TroG_slowps, "TroG/TroG_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TroGMod,
      update(TroGMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TroGMod,
      update(TroGMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TroGMod,
      update(TroGMod, .~.-pStabMean),
      test = "LRT")
anova(TroGMod,
      update(TroGMod, .~.-tStabMean),
      test = "LRT")
anova(TroGMod,
      update(TroGMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TroGMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TroGMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TroGEnv$bio1Mean_C, na.rm =T),
                                                   max(TroGEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TroGEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TroGEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TroGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TroGMod_pred_bio1Mean<-cbind(TroGMod_dat_bio1Mean,
                            predict(TroGMod,TroGMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TroGMod_pred_bio1Mean)
TroGMod_pred_bio1Mean$biomeabbr<-"TroG"
# how does the model look? plot it
#ggplot(TroGMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TroGEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TroGMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TroGEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TroGEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TroGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TroGEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TroGMod_pred_bio12Mean<-cbind(TroGMod_dat_bio12Mean,
                             predict(TroGMod,TroGMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TroGMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TroGMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TroGMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TroGEnv$pStabMean, na.rm =T),
                                                   max(TroGEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TroGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TroGEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TroGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TroGMod_pred_pStabMean<-cbind(TroGMod_dat_pStabMean, 
                             predict(TroGMod,TroGMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TroGMod_pred_pStabMean)
TroGMod_pred_pStabMean$biomeabbr<-"TroG"
# how does the model look? plot it
#ggplot(TroGMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TroGMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TroGEnv$tStabMean, na.rm =T),
                                                   max(TroGEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TroGEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TroGEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TroGEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TroGMod_pred_tStabMean<-cbind(TroGMod_dat_tStabMean,
                             predict(TroGMod,TroGMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TroGMod_pred_tStabMean)
TroGMod_pred_tStabMean$biomeabbr<-"TroG"
# how does the model look? plot it
#ggplot(TroGMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TroGMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TroGEnv$geoDist_log, na.rm =T),
                                                       max(TroGEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TroGEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TroGEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TroGEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TroGEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TroGMod_pred_geoDistMean<-cbind(TroGMod_dat_geoDistMean,
                               predict(TroGMod,TroGMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TroGMod_pred_geoDistMean)
TroGMod_pred_geoDistMean$biomeabbr<-"TroG"
# how does the model look? plot it
#ggplot(TroGMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TroGEnv)






TunEnv<-as.data.frame(fread("Tun/TunEnv.csv"))
TunEnv$bio1Mean_C<-TunEnv$bio1Mean/10
TunEnv$bio12Mean_sqrt<-sqrt(TunEnv$bio12Mean)
TunEnv$geoDist_log<-log(TunEnv$geoDist)
TunMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = TunEnv, family = "binomial")
summary(TunMod)

Tun_slowps<-as.data.frame(TunMod[["coefficients"]][-1])
Tun_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
Tun_slowps$meanpbd<-mean(TunEnv$distance)
Tun_slowps$mean_bio1<-mean(cbind((TunEnv$s1.bio1/10),(TunEnv$s2.bio1/10)), na.rm=TRUE)
Tun_slowps$mean_bio12<-mean(cbind(sqrt(TunEnv$s1.bio12),sqrt(TunEnv$s2.bio12)), na.rm=TRUE)
Tun_slowps$mean_pStab<-mean(cbind(TunEnv$s1.precipStability, TunEnv$s2.precipStability), na.rm =TRUE)
Tun_slowps$mean_tStab<-mean(cbind(TunEnv$s1.tempStability, TunEnv$s2.tempStability), na.rm = TRUE)
Tun_slowps$mean_geodist<-mean(TunEnv$geoDist_log)
Tun_slowps$biomeabbr<-"Tun"
colnames(Tun_slowps)[1]<-"coefficients"
Tun_slowps$pval<-summary(TunMod)$coefficients[2:6,4]
write.csv(Tun_slowps, "Tun/Tun_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(TunMod,
      update(TunMod, .~.-bio1Mean_C),
      test = "LRT")
anova(TunMod,
      update(TunMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(TunMod,
      update(TunMod, .~.-pStabMean),
      test = "LRT")
anova(TunMod,
      update(TunMod, .~.-tStabMean),
      test = "LRT")
anova(TunMod,
      update(TunMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(TunMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
TunMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(TunEnv$bio1Mean_C, na.rm =T),
                                                   max(TunEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(TunEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(TunEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(TunEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(TunEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TunMod_pred_bio1Mean<-cbind(TunMod_dat_bio1Mean,
                            predict(TunMod,TunMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(TunMod_pred_bio1Mean)
TunMod_pred_bio1Mean$biomeabbr<-"Tun"
# how does the model look? plot it
#ggplot(TunMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = TunEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

TunMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(TunEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(TunEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(TunEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(TunEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(TunEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TunEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
TunMod_pred_bio12Mean<-cbind(TunMod_dat_bio12Mean,
                             predict(TunMod,TunMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(TunMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(TunMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


TunMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(TunEnv$pStabMean, na.rm =T),
                                                   max(TunEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TunEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TunEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(TunEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(TunEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TunMod_pred_pStabMean<-cbind(TunMod_dat_pStabMean, 
                             predict(TunMod,TunMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(TunMod_pred_pStabMean)
TunMod_pred_pStabMean$biomeabbr<-"Tun"
# how does the model look? plot it
#ggplot(TunMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



TunMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(TunEnv$tStabMean, na.rm =T),
                                                   max(TunEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(TunEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(TunEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(TunEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(TunEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
TunMod_pred_tStabMean<-cbind(TunMod_dat_tStabMean,
                             predict(TunMod,TunMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(TunMod_pred_tStabMean)
TunMod_pred_tStabMean$biomeabbr<-"Tun"
# how does the model look? plot it
#ggplot(TunMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

TunMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(TunEnv$geoDist_log, na.rm =T),
                                                       max(TunEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(TunEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(TunEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(TunEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(TunEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
TunMod_pred_geoDistMean<-cbind(TunMod_dat_geoDistMean,
                               predict(TunMod,TunMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(TunMod_pred_geoDistMean)
TunMod_pred_geoDistMean$biomeabbr<-"Tun"
# how does the model look? plot it
#ggplot(TunMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(TunEnv)







AllEnv<-as.data.frame(fread("All/AllEnv.csv"))
AllEnv$bio1Mean_C<-AllEnv$bio1Mean/10
AllEnv$bio12Mean_sqrt<-sqrt(AllEnv$bio12Mean)
AllEnv$geoDist_log<-log(AllEnv$geoDist)
AllMod<-glm(formula = distance~bio1Mean_C+ bio12Mean_sqrt + pStabMean + tStabMean + geoDist_log, 
            data = AllEnv, family = "binomial")
summary(AllMod)

All_slowps<-as.data.frame(AllMod[["coefficients"]][-1])
All_slowps$predictor<-c("bio1Mean", "bio12Mean", "pStabMean", "tStabMean", "geoDist")
All_slowps$meanpbd<-mean(AllEnv$distance)
All_slowps$mean_bio1<-mean(cbind((AllEnv$s1.bio1/10),(AllEnv$s2.bio1/10)), na.rm=TRUE)
All_slowps$mean_bio12<-mean(cbind(sqrt(AllEnv$s1.bio12),sqrt(AllEnv$s2.bio12)), na.rm=TRUE)
All_slowps$mean_pStab<-mean(cbind(AllEnv$s1.precipStability, AllEnv$s2.precipStability), na.rm =TRUE)
All_slowps$mean_tStab<-mean(cbind(AllEnv$s1.tempStability, AllEnv$s2.tempStability), na.rm = TRUE)
All_slowps$mean_geodist<-mean(AllEnv$geoDist_log)
All_slowps$biomeabbr<-"All"
colnames(All_slowps)[1]<-"coefficients"
All_slowps$pval<-summary(AllMod)$coefficients[2:6,4]
write.csv(All_slowps, "All/All_slowps.csv", row.names=FALSE )


# check significance of individual predictor terms using likelihood ratio tests
anova(AllMod,
      update(AllMod, .~.-bio1Mean_C),
      test = "LRT")
anova(AllMod,
      update(AllMod, .~.-bio12Mean_sqrt),
      test = "LRT")
anova(AllMod,
      update(AllMod, .~.-pStabMean),
      test = "LRT")
anova(AllMod,
      update(AllMod, .~.-tStabMean),
      test = "LRT")
anova(AllMod,
      update(AllMod, .~.-geoDist_log),
      test = "LRT")

#check model residuals
# plot(AllMod)

# make one dummy data frame for each predictor variable,
# where only that one variable is allowed to vary,
# and the others are kept at their median values
AllMod_dat_bio1Mean<-data.frame("bio1Mean_C" = seq(min(AllEnv$bio1Mean_C, na.rm =T),
                                                   max(AllEnv$bio1Mean_C, na.rm =T),
                                                   length.out = 1000),
                                "bio12Mean_sqrt" = median(AllEnv$bio12Mean_sqrt, na.rm =T),
                                "pStabMean" = median(AllEnv$pStabMean, na.rm =T),
                                "tStabMean" = median(AllEnv$tStabMean, na.rm =T),
                                "geoDist_log" = median(AllEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
AllMod_pred_bio1Mean<-cbind(AllMod_dat_bio1Mean,
                            predict(AllMod,AllMod_dat_bio1Mean, type = "response", se.fit = TRUE))
names(AllMod_pred_bio1Mean)
AllMod_pred_bio1Mean$biomeabbr<-"All"
# how does the model look? plot it
#ggplot(AllMod_pred_bio1Mean, aes(x = bio1Mean_C, y = fit)) +geom_line() +
  #ylim0,1) #+
# geom_point(inherit.aes = F,
#            data = AllEnv,
#            mapping = aes(x = bio1Mean,
#                          y = distance))

AllMod_dat_bio12Mean<-data.frame("bio12Mean_sqrt" = seq(min(AllEnv$bio12Mean_sqrt, na.rm =T),
                                                       max(AllEnv$bio12Mean_sqrt, na.rm =T),
                                                       length.out = 1000),
                                 "bio1Mean_C" = median(AllEnv$bio1Mean_C, na.rm =T),
                                 "pStabMean" = median(AllEnv$pStabMean, na.rm =T),
                                 "tStabMean" = median(AllEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(AllEnv$geoDist_log, na.rm = T))

# apply predict() function to dummy data
AllMod_pred_bio12Mean<-cbind(AllMod_dat_bio12Mean,
                             predict(AllMod,AllMod_dat_bio12Mean, type = "response", se.fit = TRUE))
names(AllMod_pred_bio12Mean)
# how does the model look? plot it
#ggplot(AllMod_pred_bio12Mean, aes(x = bio12Mean_sqrt, y = fit)) +geom_line() +
  #ylim0,1) #+


AllMod_dat_pStabMean<-data.frame("pStabMean" = seq(min(AllEnv$pStabMean, na.rm =T),
                                                   max(AllEnv$pStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(AllEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(AllEnv$bio12Mean_sqrt, na.rm =T),
                                 "tStabMean" = median(AllEnv$tStabMean, na.rm =T),
                                 "geoDist_log" = median(AllEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
AllMod_pred_pStabMean<-cbind(AllMod_dat_pStabMean, 
                             predict(AllMod,AllMod_dat_pStabMean, type = "response", se.fit = TRUE))
names(AllMod_pred_pStabMean)
AllMod_pred_pStabMean$biomeabbr<-"All"
# how does the model look? plot it
#ggplot(AllMod_pred_pStabMean, aes(x = pStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+



AllMod_dat_tStabMean<-data.frame("tStabMean" = seq(min(AllEnv$tStabMean, na.rm =T),
                                                   max(AllEnv$tStabMean, na.rm =T),
                                                   length.out = 1000),
                                 "bio1Mean_C" = median(AllEnv$bio1Mean_C, na.rm =T),
                                 "bio12Mean_sqrt" = median(AllEnv$bio12Mean_sqrt, na.rm =T),
                                 "pStabMean" = median(AllEnv$pStabMean, na.rm =T),
                                 "geoDist_log" = median(AllEnv$geoDist_log, na.rm = T))
# apply predict() function to dummy data
AllMod_pred_tStabMean<-cbind(AllMod_dat_tStabMean,
                             predict(AllMod,AllMod_dat_tStabMean, type = "response", se.fit = TRUE))
names(AllMod_pred_tStabMean)
AllMod_pred_tStabMean$biomeabbr<-"All"
# how does the model look? plot it
#ggplot(AllMod_pred_tStabMean, aes(x = tStabMean, y = fit)) +geom_line() +
  #ylim0,1) #+

AllMod_dat_geoDistMean<-data.frame("geoDist_log" = seq(min(AllEnv$geoDist_log, na.rm =T),
                                                       max(AllEnv$geoDist_log, na.rm =T),
                                                       length.out = 1000),
                                   "bio1Mean_C" = median(AllEnv$bio1Mean_C, na.rm =T),
                                   "bio12Mean_sqrt" = median(AllEnv$bio12Mean_sqrt, na.rm =T),
                                   "pStabMean" = median(AllEnv$pStabMean, na.rm =T),
                                   "tStabMean" = median(AllEnv$tStabMean, na.rm = T))
# apply predict() function to dummy data
AllMod_pred_geoDistMean<-cbind(AllMod_dat_geoDistMean,
                               predict(AllMod,AllMod_dat_geoDistMean, type = "response", se.fit = TRUE))
names(AllMod_pred_geoDistMean)
AllMod_pred_geoDistMean$biomeabbr<-"All"
# how does the model look? plot it
#ggplot(AllMod_pred_geoDistMean, aes(x = geoDist_log, y = fit)) +geom_line() +
  #ylim0,1) #+
remove(AllEnv)




ModDat_bio1<-rbind(DBFMod_pred_bio1Mean, DesMod_pred_bio1Mean, FloGMod_pred_bio1Mean,
                   MangMod_pred_bio1Mean, MBFMod_pred_bio1Mean, MeditMod_pred_bio1Mean, 
                   MonGMod_pred_bio1Mean, TaiMod_pred_bio1Mean, TemBFMod_pred_bio1Mean,
                   TemCFMod_pred_bio1Mean, TemGMod_pred_bio1Mean, TroCFMod_pred_bio1Mean,
                   TroGMod_pred_bio1Mean, TunMod_pred_bio1Mean, AllMod_pred_bio1Mean)
ModDat_bio1$biomeabbr[which(ModDat_bio1$biomeabbr == "All")]<-"Global"
ModDat_bio12<-rbind(DBFMod_pred_bio12Mean, DesMod_pred_bio12Mean, FloGMod_pred_bio12Mean,
                    MangMod_pred_bio12Mean, MBFMod_pred_bio12Mean, MeditMod_pred_bio12Mean, 
                    MonGMod_pred_bio12Mean, TaiMod_pred_bio12Mean, TemBFMod_pred_bio12Mean,
                    TemCFMod_pred_bio12Mean, TemGMod_pred_bio12Mean, TroCFMod_pred_bio12Mean,
                    TroGMod_pred_bio12Mean, TunMod_pred_bio12Mean, AllMod_pred_bio12Mean)
ModDat_pStab<-rbind(DBFMod_pred_pStabMean, DesMod_pred_pStabMean, FloGMod_pred_pStabMean,
                    MangMod_pred_pStabMean, MBFMod_pred_pStabMean, MeditMod_pred_pStabMean, 
                    MonGMod_pred_pStabMean, TaiMod_pred_pStabMean, TemBFMod_pred_pStabMean,
                    TemCFMod_pred_pStabMean, TemGMod_pred_pStabMean, TroCFMod_pred_pStabMean,
                    TroGMod_pred_pStabMean, TunMod_pred_pStabMean, AllMod_pred_pStabMean)
ModDat_pStab$biomeabbr[which(ModDat_pStab$biomeabbr == "All")]<-"Global"
ModDat_tStab<-rbind(DBFMod_pred_tStabMean, DesMod_pred_tStabMean, FloGMod_pred_tStabMean,
                    MangMod_pred_tStabMean, MBFMod_pred_tStabMean, MeditMod_pred_tStabMean, 
                    MonGMod_pred_tStabMean, TaiMod_pred_tStabMean, TemBFMod_pred_tStabMean,
                    TemCFMod_pred_tStabMean, TemGMod_pred_tStabMean, TroCFMod_pred_tStabMean,
                    TroGMod_pred_tStabMean, TunMod_pred_tStabMean, AllMod_pred_tStabMean)
ModDat_tStab$biomeabbr[which(ModDat_tStab$biomeabbr == "All")]<-"Global"
ModDat_geoDist<-rbind(DBFMod_pred_geoDistMean, DesMod_pred_geoDistMean, FloGMod_pred_geoDistMean,
                    MangMod_pred_geoDistMean, MBFMod_pred_geoDistMean, MeditMod_pred_geoDistMean, 
                    MonGMod_pred_geoDistMean, TaiMod_pred_geoDistMean, TemBFMod_pred_geoDistMean,
                    TemCFMod_pred_geoDistMean, TemGMod_pred_geoDistMean, TroCFMod_pred_geoDistMean,
                    TroGMod_pred_geoDistMean, TunMod_pred_geoDistMean, AllMod_pred_geoDistMean)
ModDat_geoDist$biomeabbr[which(ModDat_geoDist$biomeabbr == "All")]<-"Global"

ModDat_bio1$biomeabbr<-as.factor(ModDat_bio1$biomeabbr)
ModDat_bio12$biomeabbr<-as.factor(ModDat_bio1$biomeabbr)
ModDat_pStab$biomeabbr<-as.factor(ModDat_pStab$biomeabbr)
ModDat_tStab$biomeabbr<-as.factor(ModDat_tStab$biomeabbr)
ModDat_geoDist$biomeabbr<-as.factor(ModDat_geoDist$biomeabbr)

ModDat_bio1 <- ModDat_bio1 %>%
  mutate(lwr = fit - 1.96 * se.fit,
         upr = fit + 1.96 * se.fit)

ModDat_bio12<- ModDat_bio12 %>%
  mutate(lwr = fit - 1.96 * se.fit,
         upr = fit + 1.96 * se.fit)

ModDat_pStab<- ModDat_pStab %>%
  mutate(lwr = fit - 1.96 * se.fit,
         upr = fit + 1.96 * se.fit)

ModDat_tStab<- ModDat_tStab %>%
  mutate(lwr = fit - 1.96 * se.fit,
         upr = fit + 1.96 * se.fit)

ModDat_geoDist<- ModDat_geoDist %>%
  mutate(lwr = fit - 1.96 * se.fit,
         upr = fit + 1.96 * se.fit)

All_slowps$biomeabbr<-"Global"
slowps<-rbind(DBF_slowps, Des_slowps, FloG_slowps, Mang_slowps,
              MBF_slowps, Medit_slowps, MonG_slowps, Tai_slowps,
              TemBF_slowps, TemCF_slowps, TemG_slowps, TroCF_slowps,
              TroG_slowps, Tun_slowps, All_slowps)

slowps$biomeabbr<-as.factor(slowps$biomeabbr)
slowps_mean_bio1<-slowps[which(slowps$predictor == "bio1Mean"), ]
slowps_mean_bio12<-slowps[which(slowps$predictor == "bio12Mean"), ]
slowps_mean_pStab<-slowps[which(slowps$predictor == "pStabMean"), ]
slowps_mean_tStab<-slowps[which(slowps$predictor == "tStabMean"), ]
slowps_mean_geoDist<-slowps[which(slowps$predictor == "geoDist"), ]

write.csv(ModDat_bio1, "ModDat_bio1.csv", row.names = FALSE)
write.csv(ModDat_bio12, "ModDat_bio12.csv", row.names = FALSE)
write.csv(ModDat_pStab, "ModDat_pStab.csv", row.names = FALSE)
write.csv(ModDat_tStab, "ModDat_tStab.csv", row.names = FALSE)
write.csv(ModDat_geoDist, "ModDat_geoDist.csv", row.names = FALSE)
write.csv(slowps_mean_bio1, "slowps_mean_bio1.csv", row.names = FALSE)
write.csv(slowps_mean_bio12, "slowps_mean_bio12.csv", row.names = FALSE)
write.csv(slowps_mean_pStab, "slowps_mean_pStab.csv", row.names = FALSE)
write.csv(slowps_mean_tStab, "slowps_mean_tStab.csv", row.names = FALSE)
write.csv(slowps_mean_geoDist, "slowps_mean_geoDist.csv", row.names = FALSE)


ModDat_bio1<-read.csv("ModDat_bio1.csv")
ModDat_bio12<-read.csv("ModDat_bio12.csv")
ModDat_pStab<-read.csv("ModDat_pStab.csv")
ModDat_tStab<-read.csv("ModDat_tStab.csv")
ModDat_geoDist<-read.csv("ModDat_geoDist.csv")
slowps_mean_bio1<-read.csv("slowps_mean_bio1.csv")
slowps_mean_bio12<-read.csv("slowps_mean_bio12.csv")
slowps_mean_pStab<-read.csv("slowps_mean_pStab.csv")
slowps_mean_tStab<-read.csv("slowps_mean_tStab.csv")
slowps_mean_geoDist<-read.csv("slowps_mean_geoDist.csv")


slowps_mean_bio1<-slowps_mean_bio1[-which(slowps_mean_bio1$biomeabbr == "Global"),]
slowps_mean_bio12<-slowps_mean_bio12[-which(slowps_mean_bio12$biomeabbr == "Global"),]
slowps_mean_pStab<-slowps_mean_bio1[-which(slowps_mean_pStab$biomeabbr == "Global"),]
slowps_mean_tStab<-slowps_mean_bio1[-which(slowps_mean_tStab$biomeabbr == "Global"),]
slowps_mean_geoDist<-slowps_mean_geoDist[-which(slowps_mean_geoDist$biomeabbr == "Global"),]


Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")

#Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "black", "#785ef0", "#dc267f", "#fe6100",
#        "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")

# Plot with confidence intervals
p1 <- ggplot(ModDat_bio1, aes(x = bio1Mean_C, y = fit, color = biomeabbr)) +
  geom_line() +  # Plot the fitted line
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = biomeabbr), alpha = 0.2, color = NA) +  # Confidence interval
  scale_color_manual(values = Hexxs) +  
  scale_fill_manual(values = Hexxs) +  # To match fill color with line color
  labs(color = "Biome", fill = "Biome", x = "Annual Temperature (Celsius)", y = "PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

p1_slope<-ggplot(slowps_mean_bio1, aes(x = mean_bio1, y = coefficients, color = biomeabbr)) +
  geom_point()+
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = "Mean Annual Temperature", y = "Slope") +
  theme_minimal()+
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")


p2 <- ggplot(ModDat_bio12, aes(x = bio12Mean_sqrt, y = fit, color = biomeabbr)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = biomeabbr), alpha = 0.2, color = NA) +  # Confidence interval 
  scale_color_manual(values = Hexxs) + 
  scale_fill_manual(values = Hexxs) +  # To match fill color with line color 
  labs(color = "Biome", 
       x = bquote(sqrt("Annual Precipitation (mm)")),  # Square root symbol in x-axis label
       y = "PBD") + 
  theme_minimal() + 
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

p2_slope<-ggplot(slowps_mean_bio12, aes(x = mean_bio12, y = coefficients, color = biomeabbr)) +
  geom_point()+
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = bquote(sqrt("Mean Annual Precipitation (mm)")), y = "Slope") +
  theme_minimal()+
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")

p3<-ggplot(ModDat_pStab, aes(x = pStabMean, y = fit, color = biomeabbr)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = biomeabbr), alpha = 0.2, color = NA) +  # Confidence interval
  scale_color_manual(values = Hexxs) +  
  scale_fill_manual(values = Hexxs) +  # To match fill color with line color
  labs(color = "Biome", x = "Precipitation Stability", y = "PBD") +
  theme_minimal() +
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")

p3_slope<-ggplot(slowps_mean_pStab, aes(x = mean_pStab, y = coefficients, color = biomeabbr)) +
  geom_point() +
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = "Mean Precipitation Stability", y = "Slope") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "right")  # Use theme() for font customization

p4<-ggplot(ModDat_tStab, aes(x = tStabMean, y = fit, color = biomeabbr)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = biomeabbr), alpha = 0.2, color = NA) +  # Confidence interval
  scale_color_manual(values = Hexxs) +  
  scale_fill_manual(values = Hexxs) +  # To match fill color with line color
  labs(color = "Biome", x = "Temperature Stability", y = "PBD") +
  theme_minimal() +
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")

p4_slope<-ggplot(slowps_mean_tStab, aes(x = mean_tStab, y = coefficients, color = biomeabbr)) +
  geom_point()+
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = "Mean Temperature Stability", y = "Slope") +
  theme_minimal()+
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")

p5<-ggplot(ModDat_geoDist, aes(x = geoDist_log, y = fit, color = biomeabbr)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = biomeabbr), alpha = 0.2, color = NA) +  # Confidence interval
  scale_color_manual(values = Hexxs) +  
  scale_fill_manual(values = Hexxs) +  # To match fill color with line color
  labs(color = "Biome", x = "log(Geographic Distance (km))", y = "PBD") +
  theme_minimal() +
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")

p5_slope<-ggplot(slowps_mean_geoDist, aes(x = mean_geodist, y = coefficients, color = biomeabbr)) +
  geom_point()+
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = "Mean log(Geographic Distance (km))", y = "Slope") +
  theme_minimal()+
  theme(text = element_text(family ="Times New Roman"), legend.position = "none")
legend<-get_legend(p3_slope)
quartz(w=3.375, h = 5)
p3_slope<-ggplot(slowps_mean_pStab, aes(x = mean_pStab, y = coefficients, color = biomeabbr)) +
  geom_point() +
  scale_color_manual(values = Hexxs) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  labs(color = "Biome", x = "Mean Precipitation Stability", y = "Slope") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")  # Use theme() for font customization


quartz(w=3.375, h = 3)
p1
quartz(w=3.375, h = 3)
p1_slope
quartz(w=3.375, h = 3)
p2
quartz(w=3.375, h = 3)
p2_slope
quartz(w=3.375, h = 3)
p3
quartz(w=3.375, h = 3)
p3_slope
quartz(w=3.375, h = 3)
p4
quartz(w=3.375, h = 3)
p4_slope
quartz(w=3.375, h = 3)
p5
quartz(w=3.375, h = 3)
p5_slope
quartz(w=3.375, h = 3)
p1_slope
quartz(w=3.375, h = 3)
p1_slope

plot.new()
legend("topright",                     # Adjust position as needed
       legend = unique(ModDat_bio1$biomeabbr),  # Biome names
       col = Hexxs,                     # Colors for each biome
       pch = 19,                        # Solid circle for dots
       title = "Biome",
       cex = 0.8)    


## Biome-scale analysis of subset of predictor variables -----

# Summary stats for bd_summary by biome
bd_stats<-bd_summary %>%
  group_by(biomeabbr) %>% 
  summarize(
    pbd_mean = mean(physim, na.rm=TRUE),
    pbd_median = median(physim, na.rm = TRUE))  %>% 
  ungroup()

# Remove the global-scale
slowps_mean_bio1<-slowps_mean_bio1[-which(slowps_mean_bio1$biomeabbr == 'Global'),]

# Add the environmental variables for each biome into the new df
bd_stats$mean_bio1<-slowps_mean_bio1$mean_bio1
bd_stats$meanbio12<-slowps_mean_bio1$mean_bio12
bd_stats$meanpStab<-slowps_mean_bio1$mean_pStab
bd_stats$meantStab<-slowps_mean_bio1$mean_tStab
bd_stats$meanGeoDist<-slowps_mean_bio1$mean_geodist
bd_stats$biomeabbr<-as.factor(bd_stats$biomeabbr)

# Fit a beta regression and perform model selection
betatest<-betareg(pbd_mean ~ mean_bio1 + meanbio12 + meanpStab + meantStab + meanGeoDist,
        data = bd_stats, na.action = na.fail)
saveadredge<-dredge(betatest)

summary(betatest)


Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")



# Plot Relationships between predictors and mean PBD


# Starting with Mean Annual Temperature
quartz(w=2.8, h = 2)
bio1_plot <- ggplot(bd_stats, aes(x = mean_bio1, y = pbd_mean, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "Mean Annual Mean Temperature", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

bio1_plot


# Mean annual precip
quartz(w=2.8, h = 2)
bio12_plot <- ggplot(bd_stats, aes(x = meanbio12, y = pbd_mean, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "Mean Annual Mean Precipitation", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

bio12_plot




# Precip Stability
quartz(w=2.8, h = 2)
pStab_plot <- ggplot(bd_stats, aes(x = meanpStab, y = pbd_mean, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "Mean Precipitation Stability", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

pStab_plot



# Temperature Stability
quartz(w=2.8, h = 2)
tStab_plot <- ggplot(bd_stats, aes(x = meantStab, y = pbd_mean, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "Mean Temperature Stability", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

tStab_plot


# Geographic Distance
quartz(w=2.8, h = 2)
geoDist_plot <- ggplot(bd_stats, aes(x = meanGeoDist, y = pbd_mean, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "Mean Geographic Distance", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")

geoDist_plot


## Connectivity regressions and plotting ----

# Combine all of the biome-specific connectivity metrics into one df
connectmets<-rbind(DBF_con, Des_con, FloG_con, Mang_con, MBF_con, Medit_con,
                   MonG_con, Tai_con, TemBF_con, TemCF_con, TemG_con, TroCF_con,
                   TroG_con, Tun_con)

remove(DBF_con, Des_con, FloG_con, Mang_con, MBF_con, Medit_con,
       MonG_con, Tai_con, TemBF_con, TemCF_con, TemG_con, TroCF_con,
       TroG_con, Tun_con)

# Format our dfs and get a df of just the biomes
connectmets$pbd<-NA
biomes<-unique(bd_summary$biomeabbr)
connectmets$biomeabbr<-as.factor(connectmets$biomeabbr)
bd_summary$biomeabbr<-as.factor(bd_summary$biomeabbr)

# Calculate the mean PBD of each biome and assign to connectmets
for(biome in biomes){
  print(biome)
  mean_pbd<-mean(bd_summary$physor[bd_summary$biomeabbr == biome], na.rm = TRUE)
  connectmets$pbd[connectmets$biomeabbr == biome]<-mean_pbd
}

# get the cluster count for each biome
connectmets$clust_count<-NA
for(biome in biomes){
  num_clust<-max(connectmets$cluster[connectmets$biomeabbr == biome])
  connectmets$clust_count[connectmets$biomeabbr == biome]<-num_clust
}

# Convert CPL to km
connectmets$CPL_km<-connectmets$CPL/1000

# Plot a pairwise relationship between the variables to determine if there are
  # any v clear patterns
pairs(~ pbd+ clust_count + CPL + SLC + MSC + meanAWF, data = connectmets)
#maybe CPL, idk


# Create a summary df for the biome-scale metrics
mets_abbr<-as.data.frame(unique(connectmets$biomeabbr))
mets_abbr$clust_count<-unique(connectmets$clust_count)
mets_abbr$CPL<-unique(connectmets$CPL_km)
mets_abbr$meanAWF<-unique(connectmets$meanAWF)
mets_abbr$SLC<-unique(connectmets$SLC)
mets_abbr$MSC<-unique(connectmets$MSC)
mets_abbr$pbd<-unique(connectmets$pbd)
mets_abbr$biomeabbr<-mets_abbr$`unique(connectmets$biomeabbr)`
mets_abbr<-mets_abbr[,-1]

# Model a beta regression
model <-betareg(pbd ~ clust_count + log(CPL) + log(SLC) + log(MSC) + log(meanAWF),
                data = mets_abbr)
summary(model)
plot(model)

# Model selection using dredge
options(na.action = "na.fail")
model_select<-dredge(model)


# Plot CPL, the only significant variable
plot(pbd~log(CPL), data=mets_abbr)
model_3<-betareg(pbd~log(CPL), data = mets_abbr)
summary(model_3)

Hexxs<-c("#6dcf6f", "#1a31a1", "#648fff", "#785ef0", "#dc267f", "#fe6100",
         "#ffb000", "#e88b4c", "#c1617f", "#712895", "#563a3c", "#83792a", "#7f8a50",  "#b6cd5f")

# Plot with confidence intervals
quartz(w=3.375, h = 3)
con_plot <- ggplot(mets_abbr, aes(x = log(CPL), y = pbd, color = biomeabbr)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  scale_color_manual(values = Hexxs) +  
  labs(color = "Biome", fill = "Biome", x = "log(CPL)", y = "Mean PBD") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"), legend.position = "none")+
  ylim(c(0.5, 0.75))
  
con_plot



## Proportion Summary Coefficient Heatmaps ----

# Read in the CSV files we made from the GDM summary outputs

# Set up the phylogenetic GDM heatmap
phylogenetic <- read.csv("PhyGDMsums.csv")
taxonomic <- read.csv("TaxGDMsums.csv")
phylogenetic <- phylogenetic[,-2]
rownames(phylogenetic) <- phylogenetic[,1]
phylogenetic <- phylogenetic[,-1]
names(phylogenetic) <- c("Geographic Distance", "Elevation","Fire frequency", "Precip. stability", "Temp. stability",
                         "Mean annual temp.","Mean dirunal temp. range", "Isothermality", "Temp. seasonality", "Max temp. warmest month",
                         "Min temp. coldest month", "Temp. annual range", "Mean temp. wettest quarter","Mean temp. driest quarter","Mean temp. warmest quarter",
                         "Mean temp. coldest quarter", "Mean annual precip.", "Precip. wettest month","Precip. driest month", "Precip. seasonality",
                         "Precip. wettest quarter","Precip. driest quarter","Precip. warmest quarter","Precip. coldest quarter")
# Test it just to see
heatmap(as.matrix(phylogenetic), Rowv = NA, Colv = NA, par(family = "Times", cex = 1))


# Set up the taxonomic GDM heatmap
taxonomic <- taxonomic[,-2]
rownames(taxonomic) <- taxonomic[,1]
taxonomic <- taxonomic[,-1]
names(taxonomic) <- c("Geographic Distance", "Elevation","Fire frequency", "Precip. stability", "Temp. stability",
                      "Mean annual temp.","Mean dirunal temp. range", "Isothermality", "Temp. seasonality", "Max temp. warmest month",
                      "Min temp. coldest month", "Temp. annual range", "Mean temp. wettest quarter","Mean temp. driest quarter","Mean temp. warmest quarter",
                      "Mean temp. coldest quarter", "Mean annual precip.", "Precip. wettest month","Precip. driest month", "Precip. seasonality",
                      "Precip. wettest quarter","Precip. driest quarter","Precip. warmest quarter","Precip. coldest quarter")
# Test it just to see
heatmap(as.matrix(taxonomic),Rowv = NA, Colv = NA, colorRampPalette(c("white", "black"))(100),main = "Taxonomic GDM")


# Create the df using the proportions, not the raw values
taxon_scaled <- taxonomic
for(i in 1:nrow(taxon_scaled)){
  taxon_scaled[i,] <- taxon_scaled[i,]/max(taxon_scaled[i,])
}

phylo_scaled <- phylogenetic
for(i in 1:nrow(phylo_scaled)){
  phylo_scaled[i,] <- phylo_scaled[i,]/max(phylo_scaled[i,])
}

# And the differences
difference <- phylo_scaled-taxon_scaled
difference <- round(difference, 1)


# Plot the difference heatmap
heatmap(as.matrix(difference), Rowv = NA, Colv = NA,
        col = colorRampPalette(c("lightblue", "white","lavender","purple"))(100),scale = "none")

image.plot(
  legend.only = TRUE,  # Only plot the legend
  col = colorRampPalette(c("lightblue", "white","lavender","purple"))(100),  # Use the same colors as the heatmap
  zlim = c(min(difference),max(difference)),  # Set the range of values to match the heatmap
  axis.args = list(at = seq(min(difference),max(difference), by = 0.2), labels = seq(min(difference),max(difference), by = 0.2)),  # Define axis labels
  legend.shrink = 0.9,  # Adjust the size of the legend
  mar = c(3, 1, 1, 1)  # Adjust the margins of the legend plot (top, bottom, left, right)
)

quartz(w=6.97, h=3)
par(mfrow = c(1,3))

# Plot the heatmaps
heatmap(as.matrix(taxon_scaled),Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "black"))(100))

heatmap(as.matrix(phylogenetic),Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "black"))(100))

heatmap(as.matrix(difference), Rowv = NA, Colv = NA,
        col = colorRampPalette(c("lightblue", "white","lavender","purple"))(100),scale = "none")

# Plot the difference legend
image.plot(
  legend.only = TRUE,  # Only plot the legend
  col = colorRampPalette(c("lightblue", "white","lavender","purple"))(100),  # Use the same colors as the heatmap
  zlim = c(min(difference),max(difference)),  # Set the range of values to match the heatmap
  axis.args = list(at = seq(min(difference),max(difference), by = 0.3), labels = seq(min(difference),max(difference), by = 0.3)),  # Define axis labels
  legend.shrink = 0.9,  # Adjust the size of the legend
  mar = c(3, 1, 1, 1)  # Adjust the margins of the legend plot (top, bottom, left, right)
)

# Plot the proportion summary coeff. legend
image.plot(
  legend.only = TRUE,  # Only plot the legend
  col = colorRampPalette(c("white", "black"))(100),  # Colors for the heatmap
  zlim = c(0, 1),  # Set the range of values
  axis.args = list(at = seq(min(difference), max(difference), by = 0.2), 
                   labels = seq(min(difference), max(difference), by = 0.2), 
                   cex.axis = 0.8),  # Customize axis labels
  legend.shrink = 0.9,  # Adjust the size of the legend
  horizontal = FALSE,  # Make it vertical (default), necessary for the left-side placement
  legend.mar = 43,  # Add space for the legend on the left
  legend.lab = ""  # No text title, avoids the error
)

## T-test of Tax vs Phy GDM results----

# Is taxonomic turnover significantly more predicted by geographic distance than phylogenetic
ttezt<-t.test(x=taxon_scaled$Geographic, y=phylo_scaled$Geographic, 
              alternative = 'greater', paired = TRUE)
ttezt

phylo_scaled1<-phylo_scaled[,-1]
phylo_long<- phylo_scaled1 %>% 
  pivot_longer(cols = 1:23, values_to = "bioclimatic")
taxon_scaled1<-taxon_scaled[,-1]
taxon_long<- taxon_scaled1 %>% 
  pivot_longer(cols = 1:23, values_to = "bioclimatic")

# Is phylogenetic turnover significantly more predicted by bioclimatic distance than taxonomic is?
ttezt1<-t.test(x=phylo_long$bioclimatic, y= taxon_long$bioclimatic,
               alternative = 'greater', paired = TRUE)
ttezt1








