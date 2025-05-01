library(tidyverse) # data wrangling
library(ChemmineR) # clustering
library(gplots)

# the goal of this script is just to cluster all of the hits from the general screen data using chemmineR

# reading in data files
hits_df<-read.csv(file = "../../data_files/PAIRs_final_data - PAIRS_final_data.csv")


################# Using chemmineR for clustering ############### 
# see link for good overview https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Clustering

# when first making this, you need an sdf file for your chemical list, you can get one at the link below based on the smiles in the hit list
smiles_hits_df <- hits_df %>%
  select(canonical_smiles) %>% 
  distinct()

write(smiles_hits_df$canonical_smiles,file = "smiles_hits.txt", sep = "\n") # writes file with all smiles in it
# I used https://cactus.nci.nih.gov/translate/ to translate smiles to SDF
all_hits_sdf <- read.SDFset("hits_sdf_file.sdf")

# now convert sdf file to atom pair set file 
hits_ap <- sdf2ap(all_hits_sdf)


# making the naming key 
names_df <- hits_df %>% 
  group_by(canonical_smiles) %>% 
  distinct(canonical_smiles,.keep_all = T) %>% 
  select(canonical_smiles, act_ingr)

names_df$AP_IDs <-hits_ap@ID 

# removing duplicate chemicals 
ap_set_hit_dups <- cmp.duplicated(hits_ap, type=1)
all_hits_sdf[which(!ap_set_hit_dups)]
hits_ap<-sdf2ap(all_hits_sdf[which(!ap_set_hit_dups)])

# identifying similar chemicals 
names_df$act_ingr[names_df$AP_IDs %in% all_hits_sdf@ID[which(ap_set_hit_dups)]]

# making an fp set from sdf
hits_fp <- desc2fp(hits_ap, descnames=1024, type="FPset")

#JP clustering ### 
# clustering used for binning
cut<-nearestNeighbors(hits_ap,numNbrs = 6)
clustering<-jarvisPatrick(cut, k=4, mode="b", linkage = "single") 
clustering_df <- data.frame(AP_IDs = names(clustering),cluster = as.integer(clustering))

# appending Jarvis Patrick cluster information to original data frame
hits_df_clust <- hits_df %>%
  left_join(names_df, by = c("canonical_smiles","act_ingr")) %>% 
  left_join(clustering_df,by=c("AP_IDs")) #%>% 

# making a heatmap
# first make distance matrix with atom pairs
dummy <- cmp.cluster(db=hits_ap, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
load("distmat.rda")
hc <- hclust(as.dist(distmat), method="single")

labels_df <- hits_df_clust %>% 
  filter(!is.na(cluster)) %>% 
  distinct(act_ingr, AP_IDs,cluster)

hc[["labels"]]<- names_df$act_ingr[names_df$AP_IDs %in% hits_ap@ID] # Assign correct item labels
pdf(file = "sensor_chems_rooted_heatmap_labeled.pdf")
heatmap.2(1-distmat, # the distance matrix is just the numeric information on distances
          dendrogram="row",
          Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), # hc determies the order of the leaves of the tree
          col=colorpanel(40, "#d73027", "#ffffbf", "#4575b4"), # change colors here, should take hexcodes
          density.info="none", 
          trace="none", 
          labCol = "", # label for columns and rows below
          labRow = hc[["labels"]],
          cexRow = 0.1,
          offsetRow = -0.2#,
          # margins = c(11,11)
)
dev.off() 
