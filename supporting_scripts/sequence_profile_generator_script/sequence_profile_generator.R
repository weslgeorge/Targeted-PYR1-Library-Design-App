library(tidyverse)

chemical_hit_data<-read.csv("../../data_files/PAIRs_final_data - PAIRS_final_data.csv")
chemical_hit_data[,sapply(chemical_hit_data,class) == "logical"] <-sapply(chemical_hit_data[,sapply(chemical_hit_data,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters

# input your desired compunds to the filter below
sequence_profile_name <- "DEET_and_friends_sequence_profile"

filtered_chemicals <- chemical_hit_data %>% 
  group_by(number,chem_name) %>% 
  filter(act_ingr %in% c("DEET","Butylparaben","Anthrone")) %>%  # enter your desired chemicals here!
  pivot_longer(cols = colnames(chemical_hit_data)[9:26],names_to = "aa_position") %>% # pivots the sequences to be long format for counting
  filter(value!="") %>% # removes empty sequence spaces
  mutate(aa_position_mutation = paste0(aa_position,value)) %>% 
  group_by(aa_position_mutation) %>% 
  reframe(chem_name,aa_position_mutation,value, count = n()) %>% # assigns a count value
  filter(count >1|chem_name == "DEET") %>% # filters for amino acid changes that showed up more than once and all the amino acid changes of a target compund ("DEET" in this example)
  distinct(aa_position_mutation,.keep_all = T)
  


# makes a dummy sequence profile table  
supplimenetal_row<-tibble(wt_aa = "x", amino_acid_position = -100, mut_aas = c("A","R","N","D","C","E","Q",
                                                                                     "G","H","I","L","K","M","F","P","S",
                                                                                     "T","W","Y","V"))
  
change_to_good_aas_format <- filtered_chemicals %>%
    ungroup() %>% 
    mutate(wt_aa = str_extract(aa_position_mutation, "[:alpha:](?=[:digit:])")) %>% # pulls aa info
    mutate(amino_acid_position = parse_number(aa_position_mutation)) %>% # pulls position info
    mutate(mut_aas = str_extract(aa_position_mutation, "(?<=[:digit:])[:alpha:]")) %>% # pulls mutation info
    select(wt_aa, amino_acid_position, mut_aas) %>%
    rbind(supplimenetal_row) %>% # binds formatting information (helps ensure 22 columns)
    pivot_wider(names_from = "mut_aas",values_from = "mut_aas",names_sort  = TRUE,values_fill = "") %>%
    filter(amino_acid_position != "-100") %>% # filters out changes form formatting df
    arrange(amino_acid_position)
  
  # last check to make sure things did not go wrong
if(length(change_to_good_aas_format) == 22){
    write.csv(change_to_good_aas_format,file = paste0(as.character(sequence_profile_name),".csv"),row.names = F)
}


