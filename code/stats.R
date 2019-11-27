library("ggplot2")
library("ggsci")

sdm_data <- read.csv("severn_barrage.csv", header=T, stringsAsFactors=F)
sdm_data_web <- read.csv("severn_barrage_linked.csv", header=T, stringsAsFactors=F)

# create mean of values for each species, along with SDs for Areas and AUC
means <- aggregate(sdm_data[, 2:4], list(sdm_data$Taxon), mean)
sds <- aggregate(sdm_data[, 2:4], list(sdm_data$Taxon), sd)

means_w <- aggregate(sdm_data_web[, 2:4], list(sdm_data_web$Taxon), mean)
sds_w <- aggregate(sdm_data_web[, 2:4], list(sdm_data_web$Taxon), sd)

# get a list of species
species <- unique(sdm_data$Taxon)
species_web <- unique(sdm_data_web$Taxon)

# let's do paired Wilcoxon tests to check differences are stat sig
dif_data <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Taxon", "p.val", "Mean", "Mean.Barrage","Mean.Diff","Perc.diff")
colnames(dif_data) <- x
for (s in species) {
    data <- subset(sdm_data, sdm_data$Taxon == s)
    wt <- wilcox.test(data$Area, data$Area.Barrage, paired = TRUE, alternative = "two.sided")
    mean.area <- mean(data$Area)
    mean.area.barrage <- mean(data$Area.Barrage)
    mean.diff <- mean.area.barrage - mean.area
    df <- data.frame(Taxon=s,p.val=wt$p.value,Mean=mean.area,Mean.Barrage=mean.area.barrage,Mean.Diff=mean.diff,Perc.diff=(mean.diff/mean.area))
    dif_data <- rbind(dif_data,df)
}

dif_data_web <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Taxon", "p.val", "Mean", "Mean.Barrage","Mean.Diff","Perc.diff")
colnames(dif_data) <- x
for (s in species_web) {
    data <- subset(sdm_data_web, sdm_data_web$Taxon == s)
    wt <- wilcox.test(data$Area, data$Area.Barrage, paired = TRUE, alternative = "two.sided")
    mean.area <- mean(data$Area)
    mean.area.barrage <- mean(data$Area.Barrage)
    mean.diff <- mean.area.barrage - mean.area
    df <- data.frame(Taxon=s,p.val=wt$p.value,Mean=mean.area,Mean.Barrage=mean.area.barrage,Mean.Diff=mean.diff,Perc.diff=(mean.diff/mean.area))
    dif_data_web <- rbind(dif_data_web,df)
}



# We now construct our master data frame so we can plot stuff
all_data <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Taxon", "State", "Area", "Normalised.Area")
colnames(all_data) <- x
for (s in species) {
    data <- subset(sdm_data, sdm_data$Taxon == s)
    area <- data.frame(area=data$Area)
    area.barrage <- data.frame(area=data$Area.Barrage)
    area$State <- "Unaltered"
    area.barrage$State <- "Barrage"
    area$Taxon <- s
    area.barrage$Taxon <- s
    # normalise areas to the unaltered mean area
    area.barrage$Normalised.Area <- (area.barrage$area - mean(area$area)) / (max(area$area) - min(area$area))
    area$Normalised.Area <- (area$area - mean(area$area)) / (max(area$area) - min(area$area))
    df <- rbind(area, area.barrage)    
    all_data <- rbind(all_data, df)
}

p <- ggplot(all_data, aes(x=Normalised.Area, fill = State)) + geom_density(colour="black",alpha=0.6) + theme_classic() + theme(text=element_text(size=20)) + scale_fill_jco()
p <- p + xlab("Normalised Area") + ylab("Density")
ggsave("SDM_normalised_area.pdf", p)

# We now construct our master data frame so we can plot stuff
all_data_web <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Taxon", "State", "Area", "Normalised.Area")
colnames(all_data_web) <- x
for (s in species_web) {
    data <- subset(sdm_data_web, sdm_data_web$Taxon == s)
    area <- data.frame(area=data$Area)
    area.barrage <- data.frame(area=data$Area.Barrage)
    area$State <- "Unaltered"
    area.barrage$State <- "Barrage"
    area$Taxon <- s
    area.barrage$Taxon <- s
    # normalise areas to the unaltered mean area
    area.barrage$Normalised.Area <- (area.barrage$area - mean(area$area)) / (max(area$area) - min(area$area))
    area$Normalised.Area <- (area$area - mean(area$area)) / (max(area$area) - min(area$area))
    df <- rbind(area, area.barrage)    
    all_data_web <- rbind(all_data_web, df)
}

p <- ggplot(all_data_web, aes(x=Normalised.Area, fill = State)) + geom_density(colour="black",alpha=0.6) + theme_classic() + theme(text=element_text(size=20)) + scale_fill_jco()
p <- p + xlab("Normalised Area") + ylab("Density")
ggsave("SDM_web_normalised_area.pdf", p)


# We now construct our master data frame so we can plot stuff
all_data.all<- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Taxon", "State", "Area", "Normalised.Area", "Type")
colnames(all_data.all) <- x
for (s in species) {
    data <- subset(sdm_data, sdm_data$Taxon == s)    
    area <- data.frame(area=data$Area)
    area.barrage <- data.frame(area=data$Area.Barrage)
    area$State <- "Unaltered"
    area.barrage$State <- "Barrage"
    area$Taxon <- s
    area.barrage$Taxon <- s
    area$Type <- "No.Web"
    area.barrage$Type <- "No.Web"
    area.barrage$Normalised.Area <- (area.barrage$area - mean(area$area)) / (max(area$area) - min(area$area))
    area$Normalised.Area <- (area$area - mean(area$area)) / (max(area$area) - min(area$area))
    df <- rbind(area, area.barrage) 
    all_data.all <- rbind(all_data.all, df)

}
for (s in species_web) {
    data.web <- subset(sdm_data_web, sdm_data_web$Taxon == s)
    area.web <- data.frame(area=data.web$Area)
    area.barrage.web <- data.frame(area=data.web$Area.Barrage)
    area.web$State <- "Unaltered"
    area.barrage.web$State <- "Barrage"
    area.web$Taxon <- s
    area.barrage.web$Taxon <- s
    area.web$Type <- "Web"
    area.barrage.web$Type <- "Web"
    area.barrage.web$Normalised.Area <- (area.barrage.web$area - mean(area.web$area)) / (max(area.web$area) - min(area.web$area))
    area.web$Normalised.Area <- (area.web$area - mean(area.web$area)) / (max(area.web$area) - min(area.web$area))
    df <- rbind(area.web, area.barrage.web)
    all_data.all <- rbind(all_data.all, df)
}
p <- ggplot(all_data.all, aes(x=Normalised.Area, fill = interaction(State,Type), )) + geom_density(colour="black",alpha=0.6) + theme_classic() + theme(text=element_text(size=20)) + scale_fill_jco(name = "Simulation", labels = c("Barrage: no food web", "Unaltered: no food web", "Barrage: with food web", "Unaltered: with food web"))
p <- p + xlab("Normalised Area") + ylab("Density")
ggsave("SDM_area_dist.pdf", width=297, height=210, units = "mm", p)


# example from velvet swimming crab
s <- "velvet_swimming_crab"
data <- subset(sdm_data, sdm_data$Taxon == s)
area <- data.frame(area=data$Area)
area.barrage <- data.frame(area=data$Area.Barrage)
area$State <- "Unaltered"
area.barrage$State <- "Barrage"
area$Taxon <- s
area.barrage$Taxon <- s
area$Type <- "No.Web"
area.barrage$Type <- "No.Web"
data.web <- subset(sdm_data_web, sdm_data_web$Taxon == s)
area.web <- data.frame(area=data.web$Area)
area.barrage.web <- data.frame(area=data.web$Area.Barrage)
area.web$State <- "Unaltered"
area.barrage.web$State <- "Barrage"
area.web$Taxon <- s
area.barrage.web$Taxon <- s
area.web$Type <- "Web"
area.barrage.web$Type <- "Web"
df <- rbind(area, area.barrage) 
df <- rbind(df, area.web)
df <- rbind(df, area.barrage.web)
p <- ggplot(df, aes(x=area, fill = interaction(State,Type), )) + geom_density(colour="black",alpha=0.6) + theme_classic() + theme(text=element_text(size=20)) + scale_fill_jco(name = "Simulation", labels = c("Barrage: no food web", "Unaltered: no food web", "Barrage: with food web", "Unaltered: with food web"))
p <- p + xlab(expression(Area~(km^2))) + ylab("Density")
ggsave("n_puber_dist.pdf", width=297, height=210, units = "mm", p)
