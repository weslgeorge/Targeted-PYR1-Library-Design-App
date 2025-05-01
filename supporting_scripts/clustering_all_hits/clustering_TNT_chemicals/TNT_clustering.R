library(tidyverse) # data wrangling
library(ChemmineR) # clutering
# library(ChemmineOB) # doesnt work 
library(gplots)
library(phyclust)
# reading in data files
all_df<-read.csv(file = "FDA_Sensor_Sequences_supporting_File 1.xlsx - combined.csv")
tnt_chems_df <- read.csv(file = "tnt_chemicals.csv")
# tnt_library<-c("2,4-Dinitrophenol")
hits_df<-all_df %>% 
  filter(Hit. == "Yes") %>% 
    filter(Library_name %in% tnt_library) %>% 
  rbind(tnt_chems_df) %>% 
  filter(Library_name!="")


################# Using chemmineR for clustering ############### 
# see link for good overview https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Clustering

# when first making this, you need an sdf file for your chemical list, you can get one at the link below based on the smiles in the hit list
smiles_hits_df <- hits_df %>%
  select(canonical_smiles)

write(smiles_hits_df$canonical_smiles,file = "smiles_hits.txt", sep = "\n") # writes file with all smiles in it
# I used https://cactus.nci.nih.gov/translate/ to translate smiles to SDF
tnt_sdf <- read.SDFset("tclcactvs000CvtOqL.sdf")
tnt_sdf@ID <- hits_df$Library_name[2:6]
# now convert sdf file to apset file 
hits_ap <- sdf2ap(tnt_sdf)

# removing duplicate chemicals 

# 
ap_set_hit_dups <- cmp.duplicated(hits_ap, type=1)
tnt_sdf[which(!ap_set_hit_dups)]
hits_ap<-sdf2ap(tnt_sdf[which(!ap_set_hit_dups)])

# making an fp set from sdf
hits_fp <- desc2fp(hits_ap, descnames=1024, type="FPset")


# making a heatmap
# first make distance matrix with atom pairs
dummy <- cmp.cluster(db=hits_ap, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
load("distmat.rda")
hc <- hclust(as.dist(distmat), method="single")

?cmp.cluster
################# Can plot here ###########################
################# Can plot here ###########################
################# Can plot here ###########################
hc[["labels"]]<- paste0(hits_df$Library_name[2:6]) # Assign correct item labels
png(filename = "TNT_chems_rooted_heatmap_labeled.png",bg = "transparent",width = 5,height = 5, units = "in",res=600)
heatmap.2(1-distmat, # the distance matrix is just the numeric information on distances
          dendrogram="row",
          # keysize=1,
          Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), # hc determies the order of the leaves of the tree
          col=colorpanel(40, "#d73027", "#ffffbf", "#4575b4"), # change colors here, should take hexcodes
          density.info="none", 
          trace="none", 
          cexRow = 0.4,
          offsetRow = -0.20,
          labCol = "", # label for columns and rows below
          labRow = hc[["labels"]]
          )
dev.off() # 

png(filename = "TNT_chems_rooted_heatmap_labeled_vibe_check.png",bg = "transparent",width = 5,height = 5, units = "in",res=600)
heatmap.2(1-distmat, # the distance matrix is just the numeric information on distances
          dendrogram="row",
          # keysize=1,
          Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), # hc determies the order of the leaves of the tree
          col=colorpanel(40, "#d73027", "#ffffbf", "#4575b4"), # change colors here, should take hexcodes
          density.info="none", 
          trace="none", 
          cexRow = 0.4,
          offsetRow = -0.20,
          labCol = "", # label for columns and rows below
          labRow = hc[["labels"]]
)
dev.off() # 




# reading in data files
all_df<-read.csv(file = "FDA_Sensor_Sequences_supporting_File 1.xlsx - combined.csv")
tnt_chems_df <- read.csv(file = "tnt_chemicals.csv")
tnt_library<-c("2,4-Dinitrophenol","2,6-Dihydroxybenzoic Acid","2,4,6-Trihydroxybenzoic Acid","Chloroxylenol")
hits_df<-all_df %>% 
  filter(Hit. == "Yes") %>% 
  filter(Library_name %in% tnt_library) %>% 
  rbind(tnt_chems_df) %>% 
  filter(Library_name!="")


################# Using chemmineR for clustering ############### 
# see link for good overview https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Clustering

# when first making this, you need an sdf file for your chemical list, you can get one at the link below based on the smiles in the hit list
smiles_hits_df <- hits_df %>%
  select(canonical_smiles)

write(smiles_hits_df$canonical_smiles,file = "smiles_sensors.txt", sep = "\n") # writes file with all smiles in it
# I used https://cactus.nci.nih.gov/translate/ to translate smiles to SDF
tnt_sdf <- read.SDFset("tnt_sensor_chems.sdf")
tnt_sdf@ID <- hits_df$Library_name
# now convert sdf file to apset file 
hits_ap <- sdf2ap(tnt_sdf)

# removing duplicate chemicals 

# 
ap_set_hit_dups <- cmp.duplicated(hits_ap, type=1)
tnt_sdf[which(!ap_set_hit_dups)]
hits_ap<-sdf2ap(tnt_sdf[which(!ap_set_hit_dups)])

# making an fp set from sdf
hits_fp <- desc2fp(hits_ap, descnames=1024, type="FPset")


# making a heatmap
# first make distance matrix with atom pairs
dummy <- cmp.cluster(db=hits_ap, cutoff=0, save.distances="distmat_sensors.rda", quiet=TRUE)
load("distmat_sensors.rda")
hc <- hclust(as.dist(distmat), method="single")


################# Can plot here ###########################
################# Can plot here ###########################
################# Can plot here ###########################
hc[["labels"]]<- paste0(hits_df$Library_name) # Assign correct item labels
png(filename = "TNT_sensor_chems_rooted_heatmap_labeled.png",bg = "transparent",width = 5,height = 10, units = "in",res=600)
heatmap.2(1-distmat, # the distance matrix is just the numeric information on distances
          dendrogram="row",
          # keysize=1,
          Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), # hc determies the order of the leaves of the tree
          col=colorpanel(40, "#d73027", "#ffffbf", "#4575b4"), # change colors here, should take hexcodes
          density.info="none", 
          trace="none", 
          cexRow = 0.4,
          offsetRow = -0.20,
          labCol = "", # label for columns and rows below
          labRow = hc[["labels"]]
)
dev.off() # 

png(filename = "TNT_sensor_chems_rooted_heatmap_labeled_vibe_check.png",bg = "transparent",width = 5,height = 5, units = "in",res=600)
heatmap.2(1-distmat, # the distance matrix is just the numeric information on distances
          dendrogram="row",
          # keysize=1,
          Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), # hc determies the order of the leaves of the tree
          col=colorpanel(40, "#d73027", "#ffffbf", "#4575b4"), # change colors here, should take hexcodes
          density.info="none", 
          trace="none", 
          cexRow = 0.4,
          offsetRow = -0.20,
          labCol = "", # label for columns and rows below
          labRow = hc[["labels"]]
)
dev.off() # 

