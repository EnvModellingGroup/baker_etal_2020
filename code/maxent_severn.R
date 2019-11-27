library(raster)
library(dismo)
library(rgeos)
library(rJava)
library(sp)
library(rgdal)
library(ncdf4)
library(kernlab)
library(grDevices)


args <- commandArgs(trailingOnly=TRUE)
temp_dir <- paste0(args[1],"/")
csv_file <- paste0("../data/",args[2])
layers <- strsplit(args[3], ",")[[1]]

# these are categorical layers that need special treatment
categorical_layers <- c(1, 3, 8)
# our mapping of numbers to layer names - note the have to be the same filenames for
# both scenarios
pd_layers <- c("always_dry_masked.tif",
               "bathy_masked.tif",
               "intertidal_masked.tif",
               "max_elev_masked.tif",
               "max_vel_masked.tif",
               "av_vel_masked.tif",
               "min_elev_masked.tif",
               "subtidal_masked.tif",
               "tidal_range_masked.tif",
               "dover_sole.tif",
               "bib.tif",
               "poor_cod.tif",
               "sand_goby.tif",
               "brown_shrimp.tif",
               "velvet_swimming_crab.tif",
               "catworm.tif",
               "whelk.tif",
               "mud_shrimp.tif",
               "periwinkle.tif",
               "lugworm.tif",
               "pink_clam.tif",
               "mud_snail.tif",
               "sand_digger_shrimp.tif"
               )

# read in our point data
locs = read.csv(csv_file, header=T, sep=",")

# build vector of file names to read in, with and without barrage
with_barrage_layers <- layers
no_barrage_layers <- layers
for (i in 1:length(no_barrage_layers)) {
    no_barrage_layers[i] <- paste0(getwd(),"/../data/unaltered/",pd_layers[as.numeric(no_barrage_layers[i])])
}
for (i in 1:length(with_barrage_layers)) {
    with_barrage_layers[i] <- paste0(getwd(),"/../data/with_barrage/",pd_layers[as.numeric(with_barrage_layers[i])])
}

# Load environmental rasters, No turbines or lagoons. The Severn as it is today.
env_layers<-stack(no_barrage_layers)
env_layers_barrage<-stack(with_barrage_layers)

# we need to set some layers as categorical
# make a vector to store this
factors <- rep(NA, sum(categorical_layers %in% layers))
if (length(factors) == 0) {
    factors <- NA
} else {
    i <- 1
    if (1 %in% layers) { 
        factors[i] <- sub('\\.tif$', '', pd_layers[1])
        i <- i+1
    }
    if (3 %in% layers) {
        factors[i] <- sub('\\.tif$', '', pd_layers[3])
        i <- i+1
    }
    if (8 %in% layers) {
        factors[i] <- sub('\\.tif$', '', pd_layers[8])
        i <- i+1
    }
}


# we then create a mask using th bathymetry layer
mask<-raster(paste0(getwd(),"/../data/mask2.tif"))

# Extract mask values to table of species co-ordinates
locs_ext=extract(mask, locs[,c("X","Y")])
# Build a data frame of species occurrence data and depth data
locs = data.frame(locs, locs_ext)
# Remove points with 0 or NA values for mask
locs <- subset(locs, !is.na(locs_ext))
#locs <- subset(locs, locs != 0) 
e = extent(mask)
# Create sequences of X and Y values to define a grid
# this same resoution as the grid
xgrid = seq(e@xmin,e@xmax,200)
ygrid = seq(e@ymin,e@ymax,200)

# Identify occurrence points within each grid cell, then draw one at random
subs = c()
for(i in 1:(length(xgrid)-1)) {
    for(j in 1:(length(ygrid)-1)) {
        gridsq = subset(locs, Y > ygrid[j] & Y < ygrid[j+1] & X > xgrid[i] & X < xgrid[i+1])
        if(dim(gridsq)[1]>0) {
            subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ])
        }
    }
}

# Assign correct co-ordinate reference system to subset
coordinates <- cbind(subs$X, subs$Y)
subs_df <- SpatialPointsDataFrame(coordinates, subs, proj4string=CRS("+proj=utm +zone=30 ellps=WGS84"))
# we create 20,000 random "background points". There are other ways to do this, but start with this.
psa <- randomPoints(mask, 20000, ext=e)

# Pull environmental data for the sumbsampled-presence points from the raster stack
presence_uk= extract(env_layers, subs_df[,c("X","Y")])
# Pull environmental data for the pseudo-absence points from the raster stack
pseudo_uk = extract(env_layers, psa)

# Build some useful dataframes, with two columns of coordinates followed by the environmental variables. For the presence points:
presence_uk = data.frame(X=subs_df$X, Y=subs_df$Y, presence_uk)
presence_uk <- na.omit(presence_uk)
presence_uk <- presence_uk[is.finite(rowSums(presence_uk)),]

# Convert psa from atomic vector matrix to data.frame
psapoints=data.frame(psa)
# Bind co-ordinates
coordinates <- cbind(psapoints$x, psapoints$y)
# Create spatial data frame of pseudo absences
psadf <- SpatialPointsDataFrame(coordinates, psapoints, proj4string=CRS("+proj=utm +zone=30 ellps=WGS84"))

# Build dataframe, with two columns of coordinates followed by the environmental variabless. For the pseudo-absences:
psadfx = psadf@coords
colnames(psadfx) = c("X","Y")
pseudo_uk = data.frame(cbind(psadfx,pseudo_uk))

# Setup for a k-fold test where k == 5
# Vector of group assignments splitting the subsampled presence points data fram with environmental data into 5 groups
group_p = kfold(presence_uk, 5)
group_a = kfold(pseudo_uk, 5)
# create output required for the loop
evaluations = list(5)
models = list(5)
# This is our k-fold test. You will want to spend a bit of time making predictions on each of the 5 sub-models
# created here to check you can make decent predictions even with missing data
for (test in 1:5) {
    # Then we use test and the kfold groupings to divide the presence and absence points:
    train_p = presence_uk[group_p!=test, c("X","Y")]
    train_a = pseudo_uk[group_a!=test, c("X","Y")]
    test_p = presence_uk[group_p==test, c("X","Y")]
    test_a = pseudo_uk[group_a==test, c("X","Y")]
    # Now, estimate a maxent model using the "training" points and the environmental data. This may take a few moments to run:
    if (is.na(factors)){
        models[test] = maxent(env_layers, p=train_p, a=train_a)
    } else {
        models[test] = maxent(env_layers, p=train_p, a=train_a, factors = factors)
    }
    # To validate the model, we use the appropriately named function.
    # Produces warning message about implicit list embedding being depreciated. May fail in future versions of R
    evaluations[test] = evaluate(test_p, test_a, models[[test]], env_layers)
}

# dump a load of output to a file
con <- file(paste0(temp_dir,"temp.log"),open="wt")
sink(con, append=TRUE)

# prints out the AUC for the k-fold tests
# ideally should be > 0.75 for all
cat("K-FOLD AUC: ")
for (test in 1:5) {
    cat(paste0(evaluations[[test]]@auc,","))
}
cat("\n")


# Assess Spatial Sorting Bias (SSB)
pres_train_me <- train_p
pres_test_me <- test_p
back_train_me <- train_a
back_test_me <- test_a
sb <- ssb(pres_test_me, back_test_me, pres_train_me)

# Adjust for SSB if present via distance based point-wise sampling
i <- pwdSample(pres_test_me, back_test_me, pres_train_me, n=1, tr=0.1)
pres_test_pwd_me <- pres_test_me[!is.na(i[,1]), ]
back_test_pwd_me <- back_test_me[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd_me, back_test_pwd_me, pres_train_me)

# Your spatial bias test = should be 1 (or close)
cat("SSB evaluation:")
cat(paste0(sb[1]/sb[2],","))
cat(paste0(sb2[1]/sb2[2],","))
cat("\n")

pres_points = presence_uk[c("X","Y")]
abs_points = pseudo_uk[c("X","Y")]
# create full maxent with all points
if (is.na(factors)) {
    model <- maxent(env_layers, p=pres_points, a=abs_points)    
} else {
    model <- maxent(env_layers, p=pres_points, a=abs_points, factors = factors)
}
# Evaluation code collected for easy recall of results back to back
evaluate_full <- evaluate(pres_points,abs_points, model, env_layers)
cat("AUC:")
cat(evaluate_full@auc)
cat("\n")

cat("Threshold:")
cat(threshold(evaluate_full)$spec_sens)
cat("\n")


# output response functions
pdf(paste0(temp_dir,"temp_response.pdf"))
response(model)
dev.off()
sink()

pred_me_uk5 = predict(model, env_layers, progress='')
writeRaster(pred_me_uk5, paste0(temp_dir,"temp_pd.tif"), options="INTERLEAVE=BAND", overwrite=TRUE)

# Now, predict a maxent model using the "trained model" and the environmental data for the Pentland Firth. This may take a few moments to run:
pred_1 = predict(model, env_layers_barrage, progress='')
writeRaster(pred_1, paste0(temp_dir,"temp_barrage.tif"), options="INTERLEAVE=BAND", overwrite=TRUE)
