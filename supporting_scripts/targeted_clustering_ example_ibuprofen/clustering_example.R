library(tidyverse) # data wrangling
library(ChemmineR) # clustering
library(gplots) # heatmap
# the goal of this script is to give an example of how you could use the PAIR datasets to inform new libraries, use this script to get a list of compounds you can select in the PAIR-Design app to make targeted libraries 
# alternatively to this script you can also use Chem mine tools online to cluster your target compounds smiles witht he PAIR dataset https://chemminetools.ucr.edu/

# before you get started make sure to set your working direcetory to this folder

# reading in general screen biosensor data file
hits_df<-read.csv(file = "../../data_files/PAIRs_final_data - PAIRS_final_data.csv")


# appending your target chemical to the generalscreen biosensor data
binding_df <- hits_df %>% 
  filter(is.na(act_ingr))
binding_df[1,] <- NA
binding_df$act_ingr <- "Ibuprofen" # add in the name of your target compound here
binding_df$canonical_smiles <- "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" # add in the smiles of your target compound here
hits_df <- hits_df %>% 
  rbind(binding_df)

################# Using chemmineR for clustering ############### 
# see link for good overview https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Clustering

# when first making this, you need an sdf file for your chemical list, you can get one at the link below based on the smiles in the hit list
smiles_hits_df <- hits_df %>%
  select(canonical_smiles,act_ingr) %>% 
  distinct(canonical_smiles,.keep_all = T)
smiles_list <-smiles_hits_df$canonical_smiles
names(smiles_list)<-smiles_hits_df$act_ingr
write(smiles_hits_df$canonical_smiles,file = "smiles_hits.txt", sep = "\n") # writes file with all smiles in it
# I used https://cactus.nci.nih.gov/translate/ to translate smiles to an SDF file, this step can also be done with ChemmineR OB but requries additional dependancies so I find this a bit easier
all_hits_sdf <- read.SDFset("ibuprofen.sdf") # download and name your .sdf file, and place it in this folder, then read the sdf file here 

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

# if you want you can make a quick heatmap to see if there are similar groups of chemcials, this can lead to different results than just finding the closest distance compounds
# first make distance matrix with atom pairs
dummy <- cmp.cluster(db=hits_ap, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
load("distmat.rda")
hc <- hclust(as.dist(distmat), method="single")

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
dev.off() # 

### If you want to just get the 5 most similar compounds based on clustering use the following code

distances_from_ibuprofen_list<-distmat[,178]

labels_df <- names_df$act_ingr[names_df$AP_IDs %in% hits_ap@ID]

number_of_compounds <- 5 # assigns the number of compounds to slice
smallest_distances_from_ibuprofen_df <- tibble(chem_name = labels_df,distance = distances_from_ibuprofen_list) %>% 
  slice_min(distances_from_ibuprofen_list, n = number_of_compounds+1) %>% # Get x number of rows rows with smallest values
  slice(2:(number_of_compounds+1))


