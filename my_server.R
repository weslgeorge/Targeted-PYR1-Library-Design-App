#### Library in required packages #### 
library(shiny)
library(doParallel)
library(spatstat.utils)
# library(tibble)
library(bioseq)
library(DT)
library(tidyverse)
library(openxlsx)

#### Read in working directory files ####
# first load in yeast favored codons table

sensors_FDA_suppliment_df <- read.csv(file = "data_files/PAIRs_final_data - PAIRS_final_data.csv")


yeast_favored_codons <-   read.csv("data_files/yeast_codons.csv") %>% 
  group_by(residue) %>% 
  arrange(residue, -freq) %>% 
  # pick the most frequently used codon
  slice(1) %>% 
  filter(residue!="*") %>% 
  select(-freq) %>% # get rid of frequency column
  ungroup() %>% 
  tibble::column_to_rownames("residue")
# 


# now the example output 
example_good_aas <- read.csv(file = "data_files/example_input_file.csv")

# next barcoding primer set options (only validated ones for now)
barcoding_primers<-read.csv(file = "data_files/pop_experiment.xlsx - 10K_primers.csv") %>% 
                  # filter(str_detect(`validated.`, "ok|OK|Ok") ) %>% 
  mutate(set_number = str_extract(X,"(\\d)+")) %>% 
  mutate(set_number = str_pad(set_number,width = 5,pad = "0",side = "left")) %>% 
  mutate(set_name = paste0("set",set_number)) %>% 
  select(!c("X","set_number","validated.")) %>% 
  arrange(set_name)

manual_option<- tibble(F_seq = "", R_seq ="",R_seq_rc="",F_Tm="",R_Tm="",Tm_avg="",F_Dg="",R_Dg="",Dg_avg="",set_name="Manual Entry")
barcoding_primers <- rbind(manual_option,barcoding_primers)

# now lets load in a list of pyr1 mutatnts that are known to be constitutive in expression
constitutive_mutant_combinations_pyr1 <- read.csv(file = "data_files/1113_3665DSMs_with_realFP.csv")

constitutive_combinations <- constitutive_mutant_combinations_pyr1 %>% # grabs just constitutive combinations in regex format
  mutate(mutation_1 = str_extract_all(Mutation,"^.*(?=_)"),mutation_2 = str_extract_all(Mutation,"(?<=_).*")) %>%
  mutate(constitutive_combination = paste0("^",mutation_1,"_",mutation_2,"_NA_NA")) %>% # use _ instead of _.* for just double mutants
  select(constitutive_combination) 
constitutive_combinations_list <- constitutive_combinations$constitutive_combination

constitutive_combinations_regex <- paste(unlist(constitutive_combinations_list), collapse = "|")

my_server <- function(input,output,session) {
  
  observeEvent(input$showSidebar, {
    shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
  })
  observeEvent(input$hideSidebar, {
    shinyjs::addClass(selector = "body", class = "sidebar-collapse")
  })
  
  
  ######## Input generation tab #################
  
  # clusters_allowed <- input$clusterchooser
  
  
  # Observe the "submit" button and create the list of selected inputs
  observeEvent(input$submit, {
    
    # # Output the selected items as a list
    # output$selected_clusters <- renderPrint({
    #   # If nothing is selected, return "No items selected"
    #   if (length(input$cluster_choices) == 0) {
    #     "No items selected"
    #   } else {
    #     # Return the selected items
    #     selected_items <- input$cluster_choices
    #     selected_items  # This is where the list of selected items is created
    #   }
    # })
    
    
    
    
    # Output the selected items as a list
    # output$selected_chemnames <- renderPrint({
    #   # If nothing is selected, return "No items selected"
    #   if (length(input$chemname_choices) == 0) {
    #     "No items selected"
    #   } else {
    #     # Return the selected items
    #     selected_items <- input$chemname_choices
    #     selected_items  # This is where the list of selected items is created
    #   }
    # })
    
    if (length(input$chemname_choices) == 0 & length(input$cluster_choices) == 0) {
      # Show a warning message as a notification
      showNotification("No items selected! Please select an item", type = "warning")
    }else{
      
      
      # creates filtered subset 
      sensors_table<-sensors_FDA_suppliment_df
      if(is.null(input$chemname_choices) == F & is.null(input$cluster_choices) == F){
        sensors_table_filtered <- sensors_table %>% 
          filter( chem_cluster %in% input$cluster_choices)
        binding_table_filtered<-sensors_table %>%
          filter(chem_name %in% input$chemname_choices)
        sensors_table_filtered<-sensors_table_filtered %>%
          rbind(binding_table_filtered) %>%
          distinct(.keep_all = T)
      }else if (is.null(input$chemname_choices) == T &is.null(input$cluster_choices) == F){
        sensors_table_filtered <- sensors_table %>% 
          filter( chem_cluster %in% input$cluster_choices)
      }else if (is.null(input$chemname_choices) == F &is.null(input$cluster_choices) == T){
        sensors_table_filtered <- sensors_table %>% 
          filter( chem_name %in% input$chemname_choices)
      }else if (is.null(input$chemname_choices) == T &is.null(input$cluster_choices) == T){
        sensors_table_filtered <- sensors_table
      }
      
      
      # pulls amino acid changes out 
      sensors_table_filtered_long<- sensors_table_filtered %>%
        group_by(clone,chem_name) %>%
        pivot_longer(cols = colnames(sensors_table_filtered)[9:26],names_to = "aa_position") %>%
        filter(value!="") %>%
        mutate(aa_position_mutation = paste0(aa_position,value)) %>%
        group_by(aa_position_mutation) %>%
        distinct(aa_position_mutation,.keep_all = T)
      
      # dummy row generation ensures correct column length in final output
      supplimenetal_row<-data_frame(wt_aa = "x", amino_acid_position = -100, mut_aas = c("A","R","N","D","C","E","Q",
                                                                                         "G","H","I","L","K","M","F","P","S",
                                                                                         "T","W","Y","V"))
      # actual input file format generation
      change_to_good_aas_format <- sensors_table_filtered_long %>%
        ungroup() %>%
        mutate(wt_aa = str_extract(aa_position_mutation, "[:alpha:](?=[:digit:])")) %>%
        mutate(amino_acid_position = parse_number(aa_position_mutation)) %>%
        mutate(mut_aas = str_extract(aa_position_mutation, "(?<=[:digit:])[:alpha:]")) %>%
        select(wt_aa, amino_acid_position, mut_aas) %>%
        rbind(supplimenetal_row) %>%
        pivot_wider(names_from = "mut_aas",values_from = "mut_aas",names_sort  = TRUE,values_fill = c("")) %>%
        filter(amino_acid_position != "-100") %>%
        arrange(amino_acid_position)
      
      
      output$sensors_FDA_suppliment_df_filtered <- renderDT({
        return(datatable(sensors_table_filtered,options = list(pageLength = nrow(table)),editable = F,class = 'cell-border stripe', rownames = F))
      })# end render DT
      
    } # end else statement
  }) # end observe event
  
  amino_acid_filtered_df <- eventReactive(input$submit, {
    
    sensors_table<-sensors_FDA_suppliment_df
    if(is.null(input$chemname_choices) == F & is.null(input$cluster_choices) == F){
      sensors_table_filtered <- sensors_table %>% 
        filter( chem_cluster %in% input$cluster_choices)
      binding_table_filtered<-sensors_table %>%
        filter(chem_name %in% input$chemname_choices)
      sensors_table_filtered<-sensors_table_filtered %>%
        rbind(binding_table_filtered) %>%
        distinct(.keep_all = T)
    }else if (is.null(input$chemname_choices) == T &is.null(input$cluster_choices) == F){
      sensors_table_filtered <- sensors_table %>% 
        filter( chem_cluster %in% input$cluster_choices)
    }else if (is.null(input$chemname_choices) == F &is.null(input$cluster_choices) == T){
      sensors_table_filtered <- sensors_table %>% 
        filter( chem_name %in% input$chemname_choices)
    }else if (is.null(input$chemname_choices) == T &is.null(input$cluster_choices) == T){
      sensors_table_filtered <- sensors_table
    }
      # pulls amino acid changes out 
      sensors_table_filtered_long<- sensors_table_filtered %>%
        group_by(clone,chem_name) %>%
        pivot_longer(cols = colnames(sensors_table_filtered)[9:26],names_to = "aa_position") %>%
        filter(value!="") %>%
        mutate(aa_position_mutation = paste0(aa_position,value)) %>%
        group_by(aa_position_mutation) %>%
        distinct(aa_position_mutation,.keep_all = T)
      
      # dummy row generation ensures correct column length in final output
      supplimenetal_row<-data_frame(wt_aa = "x", amino_acid_position = -100, mut_aas = c("A","R","N","D","C","E","Q",
                                                                                         "G","H","I","L","K","M","F","P","S",
                                                                                         "T","W","Y","V"))
      # actual input file format generation
      change_to_good_aas_format <- sensors_table_filtered_long %>%
        ungroup() %>%
        mutate(wt_aa = str_extract(aa_position_mutation, "[:alpha:](?=[:digit:])")) %>%
        mutate(amino_acid_position = parse_number(aa_position_mutation)) %>%
        mutate(mut_aas = str_extract(aa_position_mutation, "(?<=[:digit:])[:alpha:]")) %>%
        select(wt_aa, amino_acid_position, mut_aas) %>%
        rbind(supplimenetal_row) %>%
        pivot_wider(names_from = "mut_aas",values_from = "mut_aas",names_sort  = TRUE,values_fill = c("")) %>%
        filter(amino_acid_position != "-100") %>%
        arrange(amino_acid_position)
      
      
        return(change_to_good_aas_format)
     # returns good aas format if submit has been clicked 
    
  }, ignoreNULL = TRUE)
  
  
  
  # download button for example input
  output$example_input_file_download <- downloadHandler(
    filename = function() {
      paste("example_input",".csv", sep="")
    },
    content = function(file) {
      write.csv(example_good_aas, file,row.names = FALSE)
    }
  )

  ##### Now load in input files ####
  # this takes in the upload information, and if there is none, will take in the example library file
  good_aas_upload<-reactive({
    if (is.null(input$good_aas)) {
      example_good_aas[,sapply(example_good_aas,class) == "logical"] <-sapply(example_good_aas[,sapply(example_good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
      colnames(example_good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                                       "X13", "X14", "X15", "X16","X17","X18","X19","X20")
      example_good_aas[is.na(example_good_aas) == T] <- ""
      return(example_good_aas)
    } else if (!is.null(input$good_aas)){
      # actually read the file
      good_aas<-read.csv(file = input$good_aas$datapath, stringsAsFactors=FALSE)
      good_aas[,sapply(good_aas,class) == "logical"] <-sapply(good_aas[,sapply(good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
      colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                               "X13", "X14", "X15", "X16","X17","X18","X19","X20")
      good_aas[is.na(good_aas) == T] <- ""
      return(good_aas)}
  })
  
  # Generates data table of example file for readme tab
  output$good_aas_table_example <- renderDT({
      good_aas <- good_aas_upload() 
      good_aas[,sapply(good_aas,class) == "logical"] <-sapply(good_aas[,sapply(good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
      good_aas[is.na(good_aas) == T] <- ""
      colnames(good_aas) <-  c("wt_aa","amino_acid_position","A","C","D","E","F","G","H","I","K","L","M","N",
                               "P", "Q", "R", "S","T","V","W","Y")
      table <- good_aas
      return(datatable(table,options = list(pageLength = nrow(table)),
                       class = 'cell-border stripe', rownames = F))
      })
  
  
  
  
  ### Reactible data tablehttps://stackoverflow.com/questions/70155520/how-to-make-datatable-editable-in-r-shiny ####
  # makes reactive values tabel, and sets example_good_aas as the default (this will change in observe event if new file is uploaded)
  rv <- reactiveValues(data = example_good_aas, orig=example_good_aas)
  
  observeEvent(is.null(input$good_aas), {
    good_aas <- good_aas_upload() 
    good_aas[,sapply(good_aas,class) == "logical"] <-sapply(good_aas[,sapply(good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    good_aas[is.na(good_aas) == T] <- ""
    rv$data <- good_aas
  })
  
  observeEvent(input$good_aas, {
    good_aas <- good_aas_upload()
    good_aas[,sapply(good_aas,class) == "logical"] <-sapply(good_aas[,sapply(good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    good_aas[is.na(good_aas) == T] <- ""
    rv$data <- good_aas
  })
  
  observeEvent(input$submit, {
    good_aas<-amino_acid_filtered_df()
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    good_aas[,sapply(good_aas,class) == "logical"] <-sapply(good_aas[,sapply(good_aas,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
    good_aas[is.na(good_aas) == T] <- ""
    rv$data <- good_aas
  })
  
  # takes in any edits from the data table tab and updates the rv$data object with the new information
  observeEvent(input$reactive_values_dt_cell_edit, {
    row  <- input$reactive_values_dt_cell_edit$row
    clmn <- input$reactive_values_dt_cell_edit$col+1 # needs a +1 for some reason (I think bc they start counting at 0?)
    rv$data[row, clmn] <- input$reactive_values_dt_cell_edit$value
  })
  
  output$reactive_values_dt <- renderDT({
      table <- rv$data
      colnames(table) <-  c("WT AA","AA Position","A","C","D","E","F","G","H","I","K","L","M","N",
                               "P", "Q", "R", "S","T","V","W","Y")
      return(datatable(table,options = list(pageLength = nrow(table)),editable = T,class = 'cell-border stripe', rownames = F))
  })
  
  output$input_download <- downloadHandler(
    filename = function() {
      paste("input_table_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      table <- rv$data
      colnames(table) <-  c("WT AA","AA Position","A","C","D","E","F","G","H","I","K","L","M","N",
                            "P", "Q", "R", "S","T","V","W","Y")
      write.csv(table, file,row.names = FALSE)
    })
  
  ###### PYR1 TYPE #####
  pyr1_cds_type<-reactive({
    if(input$pyr1_type == "PYR1 WT"){
      return("ATGCCTTCGGAGTTAACACCAGAAGAACGATCGGAACTAAAAAACTCAATCGCCGAGTTCCACACATACCAACTCGATCCAGGAAGCTGTTCATCACTCCACGCGCAACGAATCCACGCGCCTCCGGAACTCGTCTGGTCAATCGTACGACGATTCGACAAACCACAAACATACAAACACTTCATCAAATCCTGCTCCGTCGAACAAAACTTCGAGATGCGCGTCGGATGCACGCGCGACGTGATCGTCATCAGTGGATTACCGGCGAACACATCAACGGAAAGACTCGATATACTCGACGACGAACGGAGAGTTACCGGATTCAGTATCATCGGAGGCGAACATAGGCTGACGAATTACAAATCCGTTACGACGGTGCATCGGTTCGAGAAAGAGAATCGGATCTGGACGGTGGTTTTGGAATCTTACGTCGTTGATATGCCGGAAGGTAACTCGGAGGATGATACTCGTATGTTTGCTGATACGGTTGTGAAGCTTAATTTGCAGAAACTCGCGACGGTTGCTGAAGCTATGGCTCGTAACTCCGGTGACGGAAGTGGTTCTCAGGTGACGTGAGGTCGA")
    }else if (input$pyr1_type == "PYR1*"){
      return("ATGCCTTCGGAGTTAACACCAGAAGAACGATCGGAACTAAAAAACTCAATCGCCGAGTTCCACACATACCAACTCGATCCAGGAAGCTGTTCATCACTCCACGCGCAACGAATCCACGCGCCTCCGGAACTCGTCTGGTCAATCGTACGACGATTCGACAAACCACAAACATACAAACACTTCATCAAATCCTGCTCCGTCGAACAAAACTTCGAGATGCGCGTCGGATGCACGCGCGACGTGATCGTCATCAGTGGATTACCGGCGAACACATCAACGGAAAGACTCGATATACTCGACGACGAACGGAGAGTTACCGGATTCAGTATCATCGGAGGCGAACATAGGCTGACGAATTACAAATCCGTTACGACGGTGCATCGGTTCGAGAAAGAGAATCGGATCTGGACGGTGGTTTTGGAATCTTACGTCGTTGATATGCCGGAAGGTAACTCGGAGGATGATACTCGTATGTTTGCTGATGATGTTGTGAAGCTTAATTTGCAGAAACTCGCGACGGTTGCTGAAGCTATGGCTCGTAACTCCGGTGACGGAAGTGGTTCTCAGGTGACGTGAGGTCGA")
    }
  })
  
  
  ###### Block start stop positions ####
  output$box_block1_start <- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block1_start", label = h5("Block 1 bp position start"), value = 153)
  })
  box_block1_start <-153
  output$box_block1_end <- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block1_end", label = h5("Block 1 bp position end"), value = 306)
  })
  box_block1_end <-306
  output$box_block2_start <- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block2_start", label = h5("Block 2 bp position start"), value = 303)
  })
  box_block2_start <-303
  output$box_block2_end<- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block2_end", label = h5("Block 2 bp position end"), value = 449)
  })
  box_block2_end <-449
  output$box_block3_start <- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block3_start", label = h5("Block 3 bp position start"), value = 446)
  })
  box_block3_start <-446
  output$box_block3_end <- renderUI({
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    y<-nchar(pyr1_cds_type())
    numericInput("box_block3_end", label = h5("Block 3 bp position end"), value = 572)
  })
  box_block3_end <-572

  
  ### PYR1 Constitutive filter ###
  pyr1_constitutive_filter <-TRUE
  

  
  ######## Single site mutation library ##########
  PYR1_lib_single_df<-reactive({
  req(good_aas_upload())
  good_aas <- rv$data # reads in reactive values data 
  colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                           "X13", "X14", "X15", "X16","X17","X18","X19","X20")
  
  good_aas[is.na(good_aas) == T] <- "" # ensures no NA values in table
  
  pyr1_cds <- pyr1_cds_type() # takes in cds input
  pyr1_prot<-seq_translate(dna(pyr1_cds_type())) # translates this input into protein sequence
  
  pyr1_codons <-  pyr1_cds %>% #pyr1的密码子序列 # breaks pyr1 sequence into all of its codons and makes it a dataframe
    str_extract_all("[A-Z][A-Z][A-Z]") %>%
    as.data.frame() %>%
    dplyr::select(codons=1)
  
  pyr1_codons$position <- c(1:length(pyr1_codons$codons)) # gives all the codons a position
  
  pyr1_aas <-  pyr1_prot %>% #pyr1的氨基酸序列 # aa of pyr1
    str_extract_all("[A-Z*]") %>%
    as.data.frame() %>%
    select(residues=1)
  
  
  numbering <- good_aas$amino_acid_position
  
  wt_aas <- good_aas$wt_aa
  
  # can change to be on the input file with a block column to be a bit more generalizaable, otherwise ask for an input in the UI with how many of each block there should be
  good_aas_block_specifyer <- good_aas %>% # assigns block information based on ui input 
    mutate(amino_acid_position = as.numeric(good_aas$amino_acid_position)) %>% 
    mutate(block = case_when(inside.range(amino_acid_position,c(as.numeric(box_block1_start)/3,as.numeric(box_block1_end)/3)) == T ~ 1,
                             inside.range(amino_acid_position,c(as.numeric(box_block2_start)/3,as.numeric(box_block2_end)/3)) == T ~ 2,
                             inside.range(amino_acid_position,c(as.numeric(box_block3_start)/3,as.numeric(box_block3_end)/3)) == T ~ 3))

  block1_positions <-which(good_aas_block_specifyer$block == 1) # designates what rows contain which blocks from input file
  block2_positions <-which(good_aas_block_specifyer$block == 2)
  block3_positions <-which(good_aas_block_specifyer$block == 3)
  
  #### counting block position mutant numbers
  # if statements ask if blocks exist and if ssm1 is checked from ui
  good_aas_for_counting <- good_aas[,3:22] # loads in aas of interest from cumarins table # should make into the format of a list of X##X
  if(length(block1_positions) > 0 & input$ssm1 == T){
  single_mut_number = 0
  for (i in block1_positions) {
    for (k in good_aas_for_counting[i,]) {
      if (k=="")  {next}
      single_mut_number <- single_mut_number + 1
    }}
  block1_single_number <- single_mut_number
  block1_single_mutants <- vector("list", block1_single_number) # resets the blockx_single_mutant_df
  } else {block1_single_number <- ""}
  
  if(length(block2_positions) > 0 & input$ssm2 == T){
  single_mut_number = 0
  for (i in block2_positions) {
    for (k in good_aas_for_counting[i,]) {
      if (k=="")  {next}
      single_mut_number <- single_mut_number + 1
    }}
  block2_single_number <- single_mut_number
  block2_single_mutants <- vector("list", block2_single_number)
  } else {block2_single_number <- ""}
  
  if(length(block3_positions) > 0 & input$ssm3 == T){
  single_mut_number = 0
  for (i in block3_positions) {
    for (k in good_aas_for_counting[i,]) {
      if (k=="")  {next}
      single_mut_number <- single_mut_number + 1
    }}
  block3_single_number <- single_mut_number
  block3_single_mutants <- vector("list", block3_single_number)
  } else {block3_single_number <- ""}
  
  
  
  ###
  
  ### making df with mutation combinations
  if(length(block1_positions) > 0 & input$ssm1 == T){
  counter <-1
  for (i in block1_positions) {
    for (k in good_aas[i,c(3:22)]) {
      if (k=="")  {next}
      block1_single_mutants[[counter]]  <-  c(mutant = counter,
                                              mut_1  = paste(wt_aas[i], numbering[i], k, sep="")) # mut_1 = wtaa position good aa
      counter <- counter + 1
    }}
  block1_single_mutants %<>% bind_rows() %>% mutate(mutant=as.character(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) # adds in mutant information to each mutant type
  } else {block1_single_mutants <- ""}
  
  if(length(block2_positions) > 0 & input$ssm2 == T){
  counter <-1
  for (i in block2_positions) {
    for (k in good_aas[i,c(3:22)]) {
      if (k=="")  {next}
      block2_single_mutants[[counter]]  <-  c(mutant = counter,
                                              mut_1  = paste(wt_aas[i], numbering[i], k, sep="")) # mut_1 = wtaa position good aa
      counter <- counter + 1
    }}
  block2_single_mutants %<>% bind_rows() %>% mutate(mutant=as.character(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F)
  } else {block2_single_mutants <- ""}

  if(length(block3_positions) > 0 & input$ssm3 == T){
  counter <-1
  for (i in block3_positions) {
    for (k in good_aas[i,c(3:22)]) {
      if (k=="")  {next}
      block3_single_mutants[[counter]]  <-  c(mutant = counter,
                                              mut_1  = paste(wt_aas[i], numbering[i], k, sep="")) # mut_1 = wtaa position good aa
      counter <- counter + 1
    }}
  block3_single_mutants %<>% bind_rows() %>%
    mutate(mutant=as.character(mutant)) %>%
    filter(str_detect(mut_1, "FALSE") == F)
  } else {block3_single_mutants <- ""}
  ###
  
  ### Making list of mutations 
  pyr1_codons$position <- as.character(pyr1_codons$position)
  # sets up single site mutation fill ins
  if(length(block1_positions) > 0 & input$ssm1 == T){
  block1_single_mutants$mutant<- as.character(block1_single_mutants$mutant)
  } else {block1_single_mutants <- ""}
  if(length(block2_positions) > 0 & input$ssm2 == T){
  block2_single_mutants$mutant<- as.character(block2_single_mutants$mutant)
  } else {block2_single_mutants <- ""}
  if(length(block3_positions) > 0 & input$ssm3 == T){
  block3_single_mutants$mutant<- as.character(block3_single_mutants$mutant)
  } else {block3_single_mutants <- ""}
  
  # function for making single mutant sequences
  make_single_mut_seq <- function(df, codons_list) {
    out <- df %>% 
      pivot_longer(mut_1, names_to="mut_num", values_to="mut") %>% 
      na.omit() %>% 
      tidyr::extract(mut, "(\\w)(\\d+)(\\w)", into=c("wt","position","mut")) %>%
      mutate(codons= yeast_favored_codons[mut,]) %>% #生成相应的密码子
      select(position, codons) %>% 
      rows_update(codons_list, ., by="position") %>% #注意position都需要class相同
      pull(codons) %>% str_c(collapse="")
    return(out)
  }
  if(length(block1_positions) > 0 & input$ssm1 == T){
  output1 <- foreach(i=1:nrow(block1_single_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
    make_single_mut_seq(block1_single_mutants[i,], pyr1_codons)}
  } else {block1_single_mutants <- ""}
  if(length(block2_positions) > 0 & input$ssm2 == T){
  output2 <- foreach(i=1:nrow(block2_single_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
    make_single_mut_seq(block2_single_mutants[i,], pyr1_codons)}
  } else {block2_single_mutants <- ""}
  if(length(block3_positions) > 0 & input$ssm3 == T){
  output3 <- foreach(i=1:nrow(block3_single_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
    make_single_mut_seq(block3_single_mutants[i,], pyr1_codons)} # generates sequences for each mutation
  } else {block3_single_mutants <- ""}
  
  if(length(block1_positions) > 0 & input$ssm1 == T){
    if(length(output1)>1){ # if statement fixes bug where if the output is only a single line then the output[,1] calls for a dimension that does not exist
  single_block1 <- cbind(block1_single_mutants, output1[,1]) %>% tibble()
  colnames(single_block1)[3] <- "mut_cds" #修改第3列的列名, adds column information
    } else{
      single_block1 <- cbind(block1_single_mutants, output1) # reattaches the sequence information back to block information sheet
      colnames(single_block1)[3] <- "mut_cds" #修改第3列的列名
    }
  } else {single_block1 <- ""}
  if(length(block2_positions) > 0 & input$ssm2 == T){
    if(length(output2)>1){
  single_block2 <- cbind(block2_single_mutants, output2[,1]) %>% tibble()
  colnames(single_block2)[3] <- "mut_cds" #修改第3列的列名
    } else{
      single_block2 <- cbind(block2_single_mutants, output2) # reattaches the sequence information back to block information sheet
      colnames(single_block2)[3] <- "mut_cds" #修改第3列的列名
    }
  } else {single_block2 <- ""}
  if(length(block3_positions) > 0 & input$ssm3 == T){
    if(length(output3)>1){
  single_block3 <- cbind(block3_single_mutants, output3[,1]) %>% tibble() # reattaches the sequence information back to block information sheet
  colnames(single_block3)[3] <- "mut_cds" #修改第3列的列名
    } else {
      single_block3 <- cbind(block3_single_mutants, output3) # reattaches the sequence information back to block information sheet
      colnames(single_block3)[3] <- "mut_cds" #修改第3列的列名
      }
  } else {single_block3 <- ""}
  
  if(length(block1_positions) > 0 & input$ssm1 == T){
  single_oligo_block1 <- single_block1 %>%
    mutate(block1= str_sub(mut_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end)  )) %>% # Need to talk to hao about this, I think it has something to do with the size of the blocks and I can use that information pretty well to make it generalizable but I will keep static for now
    pivot_longer(block1,names_to="block", values_to="oligo")
  } else {single_oligo_block1 <- "0"}
  if(length(block2_positions) > 0 & input$ssm2 == T){
  single_oligo_block2 <- single_block2 %>%
    mutate(block2= str_sub(mut_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end)  )) %>% #生成新列block2
    pivot_longer(block2,names_to="block", values_to="oligo")
  } else {single_oligo_block2 <- "0"}
  if(length(block3_positions) > 0 & input$ssm3 == T){
  single_oligo_block3 <- single_block3 %>%
    mutate(block3= str_sub(mut_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end)  )) %>% #生成新列block3
    pivot_longer(block3,names_to="block", values_to="oligo") # creates oligos for each sequence
  } else {single_oligo_block3 <- "0"}
  
  # probably could have done the if statments here much more cleanly but basically it just ensures the output has all blocks that are present from ui input
  if (length(single_oligo_block1) >= 5 & length(single_oligo_block2) >= 5 & length(single_oligo_block3) >= 5){
    PYR1_lib_single <- rbind(single_oligo_block1, single_oligo_block2,
                             single_oligo_block3) %>%  # binds all blocks into one file
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  } else if(length(single_oligo_block1) < 5 & length(single_oligo_block2) >= 5 & length(single_oligo_block3) >= 5){
    PYR1_lib_single <- rbind(single_oligo_block2,
                             single_oligo_block3) %>%  # binds all blocks into one file
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  }  else if(length(single_oligo_block1) >= 5 & length(single_oligo_block2) < 5 & length(single_oligo_block3) >= 5){
    PYR1_lib_single <- rbind(single_oligo_block1,
                             single_oligo_block3) %>%  # binds all blocks into one file
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  } else if(length(single_oligo_block1) >= 5 & length(single_oligo_block2) >= 5 & length(single_oligo_block3) < 5){
    PYR1_lib_single <- rbind(single_oligo_block1,
                             single_oligo_block2) %>%  # binds all blocks into one file
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  }  else if(length(single_oligo_block1) < 5 & length(single_oligo_block2) < 5 & length(single_oligo_block3) >= 5){
    PYR1_lib_single <- single_oligo_block3 %>%
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  } else if(length(single_oligo_block1) < 5 & length(single_oligo_block2) >= 5 & length(single_oligo_block3) < 5){
    PYR1_lib_single <- single_oligo_block2 %>%
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  } else if(length(single_oligo_block1) >= 5 & length(single_oligo_block2) < 5 & length(single_oligo_block3) < 5){
    PYR1_lib_single <- single_oligo_block1 %>%
      add_column(mut_2 = NA, mut_3 = NA,mut_4 = NA) %>% 
      select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
  } else if(length(single_oligo_block1) < 5 & length(single_oligo_block2) < 5 & length(single_oligo_block3) < 5){
    PYR1_lib_single <- data.frame(warning = "No selected AA changes in any block")
  }
  
  return(PYR1_lib_single)
  })
  
  ######## Doulble site mutation ###########
  # see Single site mutatoin for comment examples, its a very similar chunk of code
  PYR1_lib_double_df <- reactive({
    req(good_aas_upload())
    good_aas <- rv$data
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    
    good_aas[is.na(good_aas) == T] <- F
    
    pyr1_cds <- pyr1_cds_type()
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    
    pyr1_codons <-  pyr1_cds %>% #pyr1的密码子序列 # breaks pyr1 sequence into all of its codons and makes it a dataframe
      str_extract_all("[A-Z][A-Z][A-Z]") %>%
      as.data.frame() %>%
      dplyr::select(codons=1)
    
    pyr1_codons$position <- c(1:length(pyr1_codons$codons)) # gives all the codons a position
    
    pyr1_aas <-  pyr1_prot %>% #pyr1的氨基酸序列 # aa of pyr1
      str_extract_all("[A-Z*]") %>%
      as.data.frame() %>%
      select(residues=1)
    
    # can change numbering to unique(good_aas$position)
    numbering <- good_aas$amino_acid_position
    # can change wt aas to good_aas %>% group_by(wt_aa, position) %>% distinct(position, .keep_all = T) %>% select(wt_aas) %>% as.list()
    wt_aas <- good_aas$wt_aa
    
    good_aas_block_specifyer <- good_aas %>%
      mutate(amino_acid_position = as.numeric(good_aas$amino_acid_position),
             block = case_when(inside.range(amino_acid_position,c(as.numeric(box_block1_start)/3,as.numeric(box_block1_end)/3)) == T ~ 1,
                               inside.range(amino_acid_position,c(as.numeric(box_block2_start)/3,as.numeric(box_block2_end)/3)) == T ~ 2,
                               inside.range(amino_acid_position,c(as.numeric(box_block3_start)/3,as.numeric(box_block3_end)/3)) == T ~ 3))
    
    block1_positions <-which(good_aas_block_specifyer$block == 1)
    block2_positions <-which(good_aas_block_specifyer$block == 2)
    block3_positions <-which(good_aas_block_specifyer$block == 3)
    
    if(length(block1_positions) > 1 & input$dsm1 == T){
    double_mut_number = 0
    for (i in block1_positions) { 
      for (j in block1_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              double_mut_number <- double_mut_number + 1
            }}}}}
    block1_double_number <- double_mut_number
    block1_double_mutants <- vector("list", block1_double_number)
    } else {block1_double_number <- ""}
    
    if(length(block2_positions) > 1 & input$dsm2 == T){
    double_mut_number = 0
    for (i in block2_positions) { 
      for (j in block2_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              double_mut_number <- double_mut_number + 1
            }}}}}
    block2_double_number <- double_mut_number
    block2_double_mutants <- vector("list", block2_double_number)
    } else {block2_double_number <- ""}
    
    if(length(block3_positions) > 1 & input$dsm3 == T){
    double_mut_number = 0
    for (i in block3_positions) { 
      for (j in block3_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              double_mut_number <- double_mut_number + 1
            }}}}}
    block3_double_number <- double_mut_number
    block3_double_mutants <- vector("list", block3_double_number)
    } else {block3_double_number <- ""}
    
    if(length(block1_positions) > 1 & input$dsm1 == T){
    counter <-1
    for (i in block1_positions) {
      for (j in block1_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              block1_double_mutants[[counter]]  <-  c(
                mutant = counter,
                mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                mut_2  = paste(wt_aas[j], numbering[j], l, sep=""))
              counter <- counter + 1
            }}}}}
    block1_double_mutants %<>% bind_rows() %>% mutate(mutant=as.character(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F)
    } else {block1_double_mutants <- ""}
    
    if(length(block2_positions) > 1 & input$dsm2 == T){
    counter <-1
    for (i in block2_positions) {
      for (j in block2_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              block2_double_mutants[[counter]]  <-  c(
                mutant = counter,
                mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                mut_2  = paste(wt_aas[j], numbering[j], l, sep=""))
              counter <- counter + 1
            }}}}}
    block2_double_mutants %<>% bind_rows() %>% mutate(mutant=as.character(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F)
    } else {block2_double_mutants <- ""}
    
    if(length(block3_positions) > 1 & input$dsm3 == T){
    counter <-1
    for (i in block3_positions) {
      for (j in block3_positions) {
        if (j>i) {
          for (k in good_aas[i,c(3:22)]) {
            if (k=="")  {next}
            for (l in good_aas[j,c(3:22)]) {
              if (l=="") {next}
              block3_double_mutants[[counter]]  <-  c(
                mutant = counter,
                mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                mut_2  = paste(wt_aas[j], numbering[j], l, sep=""))
              counter <- counter + 1
            }}}}}
    block3_double_mutants %<>% bind_rows() %>% 
      mutate(mutant=as.character(mutant))
    } else {block3_double_mutants <- ""}
    
    make_double_mut_seq <- function(df, codons_list) {
      out <- df %>% 
        pivot_longer(mut_1:mut_2, names_to="mut_num", values_to="mut") %>% 
        na.omit() %>% 
        tidyr::extract(mut, "(\\w)(\\d+)(\\w)", into=c("wt","position","mut")) %>%
        #拆分为三个部分w = word；d+means 数字
        #每次只读取一行，所以一次就是2个突变
        mutate(codons= yeast_favored_codons[mut,]) %>% #生成相应的密码子
        select(position, codons) %>% 
        rows_update(codons_list, ., by="position") %>% #注意position都需要class相同
        pull(codons) %>% str_c(collapse="")
      return(out)
    }
    
    pyr1_codons$position <- as.character(pyr1_codons$position)
    if(length(block1_positions) > 1 & input$dsm1 == T){
    block1_double_mutants$mutant<- as.character(block1_double_mutants$mutant)
    output1 <- foreach(i=1:nrow(block1_double_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_double_mut_seq(block1_double_mutants[i,], pyr1_codons)}
    if(length(output1)>1){
    double_block1 <- cbind(block1_double_mutants, output1[,1]) %>% tibble()
    colnames(double_block1)[4] <- "mut_cds" #修改第4列的列名
    } else{
      double_block1 <- cbind(block1_double_mutants, output1) %>% tibble()
      colnames(double_block1)[4] <- "mut_cds" #修改第4列的列名
    }
    double_oligo_block1 <- double_block1 %>% 
      mutate(block1= str_sub(mut_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end)  )) %>% #生成新列block1
      pivot_longer(block1,names_to="block", values_to="oligo")
    } else {double_oligo_block1 <- "0"}
    
    if(length(block2_positions) > 1 & input$dsm2 == T){
    block2_double_mutants$mutant<- as.character(block2_double_mutants$mutant)
    output2 <- foreach(i=1:nrow(block2_double_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_double_mut_seq(block2_double_mutants[i,], pyr1_codons)}
    if(length(output2)>1){
    double_block2 <- cbind(block2_double_mutants, output2[,1]) %>% tibble()
    colnames(double_block2)[4] <- "mut_cds" #修改第4列的列名
    } else{
      double_block2 <- cbind(block2_double_mutants, output2) %>% tibble()
      colnames(double_block2)[4] <- "mut_cds" #修改第4列的列名
    }
    double_oligo_block2 <- double_block2 %>% 
      mutate(block2= str_sub(mut_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end)  )) %>% #生成新列block2
      pivot_longer(block2,names_to="block", values_to="oligo")
    } else {double_oligo_block2 <- "0"}
    
    if(length(block3_positions) > 1 & input$dsm3 == T){
    block3_double_mutants$mutant<- as.character(block3_double_mutants$mutant)
    output3 <- foreach(i=1:nrow(block3_double_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_double_mut_seq(block3_double_mutants[i,], pyr1_codons)}
    if(length(output3)>1){
    double_block3 <- cbind(block3_double_mutants, output3[,1]) %>% tibble()
    colnames(double_block3)[4] <- "mut_cds" #修改第4列的列名
    } else{
      double_block3 <- cbind(block3_double_mutants, output3) %>% tibble()
      colnames(double_block3)[4] <- "mut_cds" #修改第4列的列名
    }
    double_oligo_block3 <- double_block3 %>% 
      mutate(block3= str_sub(mut_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end)  )) %>% #生成新列block3
      pivot_longer(block3,names_to="block", values_to="oligo")
    } else {double_oligo_block3 <- "0"}

    if (length(double_oligo_block1) >= 5 & length(double_oligo_block2) >= 5 & length(double_oligo_block3) >= 5){
      PYR1_lib_double <- rbind(double_oligo_block1, double_oligo_block2,
                               double_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_3 = NA,mut_4 = NA) %>% 
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(double_oligo_block1) < 5 & length(double_oligo_block2) >= 5 & length(double_oligo_block3) >= 5){
      PYR1_lib_double <- rbind(double_oligo_block2,
                               double_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_3 = NA,mut_4 = NA) %>%          
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(double_oligo_block1) >= 5 & length(double_oligo_block2) < 5 & length(double_oligo_block3) >= 5){
      PYR1_lib_double <- rbind(double_oligo_block1,
                               double_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_3 = NA,mut_4 = NA) %>%          
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(double_oligo_block1) >= 5 & length(double_oligo_block2) >= 5 & length(double_oligo_block3) < 5){
      PYR1_lib_double <- rbind(double_oligo_block1,
                               double_oligo_block2) %>%  # binds all blocks into one file
        add_column(mut_3 = NA,mut_4 = NA) %>%          
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(double_oligo_block1) < 5 & length(double_oligo_block2) < 5 & length(double_oligo_block3) >= 5){
      PYR1_lib_double <- double_oligo_block3 %>%
        add_column(mut_3 = NA,mut_4 = NA) %>%       
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(double_oligo_block1) < 5 & length(double_oligo_block2) >= 5 & length(double_oligo_block3) < 5){
      PYR1_lib_double <- double_oligo_block2 %>%
        add_column(mut_3 = NA,mut_4 = NA) %>%          
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(double_oligo_block1) >= 5 & length(double_oligo_block2) < 5 & length(double_oligo_block3) < 5){
      PYR1_lib_double <- double_oligo_block1 %>%
        add_column(mut_3 = NA,mut_4 = NA) %>%          
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(double_oligo_block1) < 5 & length(double_oligo_block2) < 5 & length(double_oligo_block3) < 5){
      PYR1_lib_double <- data.frame(warning = "No selected AA changes in any block")
    }
    return(PYR1_lib_double)
  })
  
  ######## Triple site mutation library #######
  # see Single site mutation library for comments as it is very similar in format
  PYR1_lib_triple_df <- reactive({
    req(good_aas_upload())
    good_aas <- rv$data
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    
    good_aas[is.na(good_aas) == T] <- ""
    
    pyr1_cds <- pyr1_cds_type()
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    
    pyr1_codons <-  pyr1_cds %>% #pyr1的密码子序列 # breaks pyr1 sequence into all of its codons and makes it a dataframe
      str_extract_all("[A-Z][A-Z][A-Z]") %>%
      as.data.frame() %>%
      dplyr::select(codons=1)
    
    pyr1_codons$position <- c(1:length(pyr1_codons$codons)) # gives all the codons a position
    
    pyr1_aas <-  pyr1_prot %>% #pyr1的氨基酸序列 # aa of pyr1
      str_extract_all("[A-Z*]") %>%
      as.data.frame() %>%
      select(residues=1)
    
    # can change numbering to unique(good_aas$position)
    numbering <- good_aas$amino_acid_position
    # can change wt aas to good_aas %>% group_by(wt_aa, position) %>% distinct(position, .keep_all = T) %>% select(wt_aas) %>% as.list()
    wt_aas <- good_aas$wt_aa
    
    good_aas_block_specifyer <- good_aas %>%
      mutate(amino_acid_position = as.numeric(good_aas$amino_acid_position),
             block = case_when(inside.range(amino_acid_position,c(as.numeric(box_block1_start)/3,as.numeric(box_block1_end)/3)) == T ~ 1,
                               inside.range(amino_acid_position,c(as.numeric(box_block2_start)/3,as.numeric(box_block2_end)/3)) == T ~ 2,
                               inside.range(amino_acid_position,c(as.numeric(box_block3_start)/3,as.numeric(box_block3_end)/3)) == T ~ 3))
    
    block1_positions <-which(good_aas_block_specifyer$block == 1)
    block2_positions <-which(good_aas_block_specifyer$block == 2)
    block3_positions <-which(good_aas_block_specifyer$block == 3)
    
    #_______________work well===generate triple mutants
    
    #生成3突CDS——————————————————————————————
    make_mutant_sequence <- function(df, codons_list) {
      out <- df %>% 
        pivot_longer(mut_1:mut_3, names_to="mut_num", values_to="mut") %>% 
        na.omit() %>% 
        tidyr::extract(mut, "(\\w)(\\d+)(\\w)", into=c("wt","position","mut")) %>%
        #拆分为三个部分w = word；d+means 数字
        #每次只读取一行，所以一次就是三个突变
        mutate(codons= yeast_favored_codons[mut,]) %>% #生成相应的密码子
        select(position, codons) %>% 
        rows_update(codons_list, ., by="position") %>% #注意position都需要class相同
        pull(codons) %>% str_c(collapse="")
      return(out)
    }
    
    pyr1_codons$position <- as.character(pyr1_codons$position)
    
    if(length(block1_positions) > 2 & input$tsm1 == T){
    triple_mut_number = 0
    for (i in block1_positions) { 
      for (j in block1_positions) {
        if (j>i) {
          for (m in block1_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    block1_triple_number <- triple_mut_number + 1
                  }}}}}}}}
    counter <-1
    block1_triple_mutants <- vector("list", block1_triple_number)# build a vector with many elements.
    for (i in block1_positions) { 
      for (j in block1_positions) {
        if (j>i) {
          for (m in block1_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    
                    block1_triple_mutants[[counter]]  <-  c(
                      mutant = counter,
                      mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                      mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                      mut_3  = paste(wt_aas[m], numbering[m], n, sep="") 
                    )
                    counter <- counter + 1
                    
                  }}}}}}}}
    block1_triple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant))%>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
    
    block1_triple_mutants$mutant<- as.character(block1_triple_mutants$mutant)
    
    output1 <- foreach(i=1:nrow(block1_triple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_mutant_sequence(block1_triple_mutants[i,], pyr1_codons)
    }
    
    if(length(output1)>1){
      triple_block1 <- cbind(block1_triple_mutants, output1[,1]) %>% tibble()
      colnames(triple_block1)[5] <- "mut_cds" #修改第5列的列名
    } else{
      triple_block1 <- cbind(block1_triple_mutants, output1) %>% tibble()
      colnames(triple_block1)[5] <- "mut_cds" #修改第5列的列名
    }
    
    triple_oligo_block1 <- triple_block1 %>% 
      mutate(block1= str_sub(mut_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) )) %>% #生成新列block1
      pivot_longer(block1,names_to="block", values_to="oligo")
    
    } else {triple_oligo_block1 <- "0"}
    
    
    if(length(block2_positions) > 2 & input$tsm2 == T){
    triple_mut_number = 0
    for (i in block2_positions) { 
      for (j in block2_positions) {
        if (j>i) {
          for (m in block2_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    block2_triple_number <- triple_mut_number + 1
                  }}}}}}}}
    counter <-1
    block2_triple_mutants <- vector("list", block2_triple_number)# build a vector with many elements.
    for (i in block2_positions) { 
      for (j in block2_positions) {
        if (j>i) {
          for (m in block2_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    
                    block2_triple_mutants[[counter]]  <-  c(
                      mutant = counter,
                      mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                      mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                      mut_3  = paste(wt_aas[m], numbering[m], n, sep="") 
                    )
                    counter <- counter + 1
                    
                  }}}}}}}}
    block2_triple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
    
    block2_triple_mutants$mutant<- as.character(block2_triple_mutants$mutant)
    
    output2 <- foreach(i=1:nrow(block2_triple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_mutant_sequence(block2_triple_mutants[i,], pyr1_codons)
    }
    
    if(length(output2)>1){
      triple_block2 <- cbind(block2_triple_mutants, output2[,1]) %>% tibble()
      colnames(triple_block2)[5] <- "mut_cds" #修改第5列的列名
    } else{
      triple_block2 <- cbind(block2_triple_mutants, output2) %>% tibble()
      colnames(triple_block2)[5] <- "mut_cds" #修改第5列的列名
    }
    
    triple_oligo_block2 <- triple_block2 %>% 
      mutate(block2= str_sub(mut_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end)  )) %>%
      pivot_longer(block2,names_to="block", values_to="oligo")
    
    } else {triple_oligo_block2 <- "0"}
    
    if(length(block3_positions) > 2 & input$tsm3 == T){
    triple_mut_number = 0
    for (i in block3_positions) { 
      for (j in block3_positions) {
        if (j>i) {
          for (m in block3_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    block3_triple_number <- triple_mut_number + 1
                  }}}}}}}}

    counter <-1
    block3_triple_mutants <- vector("list", block3_triple_number)# build a vector with many elements.
    for (i in block3_positions) { 
      for (j in block3_positions) {
        if (j>i) {
          for (m in block3_positions) {
            if (m>j) {
              for (k in good_aas[i,c(3:22)]) {
                if (k=="")  {next}
                for (l in good_aas[j,c(3:22)]) {
                  if (l=="") {next}
                  for (n in good_aas[m,c(3:22)]) {
                    if (n=="") {next}
                    
                    block3_triple_mutants[[counter]]  <-  c(
                      mutant = counter,
                      mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                      mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                      mut_3  = paste(wt_aas[m], numbering[m], n, sep="") 
                    )
                    counter <- counter + 1
                    
                  }}}}}}}}
    block3_triple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
    #here to modify the position column to as.numeric
    # create coding sequences for all of the triple mutants
    #把pyr1_codons以及block_triple_mutants中的position都改成as.character
    block3_triple_mutants$mutant<- as.character(block3_triple_mutants$mutant)

    output3 <- foreach(i=1:nrow(block3_triple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
      make_mutant_sequence(block3_triple_mutants[i,], pyr1_codons)
    }

    
    if(length(output3)>1){
      triple_block3 <- cbind(block3_triple_mutants, output3[,1]) %>% tibble()
      colnames(triple_block3)[5] <- "mut_cds" #修改第5列的列名
    } else{
      triple_block3 <- cbind(block3_triple_mutants, output3) %>% tibble()
      colnames(triple_block3)[5] <- "mut_cds" #修改第5列的列名
    }
    #先生成oligo再合并是最好的
    #生成3突oligos——————————————————————————————
    triple_oligo_block3 <- triple_block3 %>% 
      mutate(block3= str_sub(mut_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end)  )) %>%
      pivot_longer(block3,names_to="block", values_to="oligo")
    
    } else {triple_oligo_block3 <- "0"}
    
    if (length(triple_oligo_block1) >= 5 & length(triple_oligo_block2) >= 5 & length(triple_oligo_block3) >= 5){
      PYR1_lib_triple <- rbind(triple_oligo_block1, triple_oligo_block2,
                               triple_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(triple_oligo_block1) < 5 & length(triple_oligo_block2) >= 5 & length(triple_oligo_block3) >= 5){
      PYR1_lib_triple <- rbind(triple_oligo_block2,
                               triple_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(triple_oligo_block1) >= 5 & length(triple_oligo_block2) < 5 & length(triple_oligo_block3) >= 5){
      PYR1_lib_triple <- rbind(triple_oligo_block1,
                               triple_oligo_block3) %>%  # binds all blocks into one file
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(triple_oligo_block1) >= 5 & length(triple_oligo_block2) >= 5 & length(triple_oligo_block3) < 5){
      PYR1_lib_triple <- rbind(triple_oligo_block1,
                               triple_oligo_block2) %>%  # binds all blocks into one file
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(triple_oligo_block1) < 5 & length(triple_oligo_block2) < 5 & length(triple_oligo_block3) >= 5){
      PYR1_lib_triple <- triple_oligo_block3 %>%
        add_column(mut_4 = NA) %>% 
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(triple_oligo_block1) < 5 & length(triple_oligo_block2) >= 5 & length(triple_oligo_block3) < 5){
      PYR1_lib_triple <- triple_oligo_block2 %>%
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(triple_oligo_block1) >= 5 & length(triple_oligo_block2) < 5 & length(triple_oligo_block3) < 5){
      PYR1_lib_triple <- triple_oligo_block1 %>%
        add_column(mut_4 = NA) %>%
        select(block, "mutation_1"=mut_1, "mutation_2"= mut_2,"mutation_3" = mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(triple_oligo_block1) < 5 & length(triple_oligo_block2) < 5 & length(triple_oligo_block3) < 5){
      PYR1_lib_triple <- data.frame(warning = "No selected AA changes in any block")
    }
    
    
    return(PYR1_lib_triple)
  })
  
  ######## Quadruple site mutation library #######
  # see Single site mutation library for comments as it is very similar in format
  PYR1_lib_Quadruple_df <- reactive({
    req(good_aas_upload())
    good_aas <- rv$data
    colnames(good_aas) <-  c("wt_aa","amino_acid_position","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12",
                             "X13", "X14", "X15", "X16","X17","X18","X19","X20")
    
    good_aas[is.na(good_aas) == T] <- ""
    
    pyr1_cds <- pyr1_cds_type()
    pyr1_prot<-seq_translate(dna(pyr1_cds_type()))
    
    pyr1_codons <-  pyr1_cds %>% #pyr1的密码子序列 # breaks pyr1 sequence into all of its codons and makes it a dataframe
      str_extract_all("[A-Z][A-Z][A-Z]") %>%
      as.data.frame() %>%
      dplyr::select(codons=1)
    
    pyr1_codons$position <- c(1:length(pyr1_codons$codons)) # gives all the codons a position
    
    pyr1_aas <-  pyr1_prot %>% #pyr1的氨基酸序列 # aa of pyr1
      str_extract_all("[A-Z*]") %>%
      as.data.frame() %>%
      select(residues=1)
    
    # can change numbering to unique(good_aas$position)
    numbering <- good_aas$amino_acid_position
    # can change wt aas to good_aas %>% group_by(wt_aa, position) %>% distinct(position, .keep_all = T) %>% select(wt_aas) %>% as.list()
    wt_aas <- good_aas$wt_aa
    
    good_aas_block_specifyer <- good_aas %>%
      mutate(amino_acid_position = as.numeric(good_aas$amino_acid_position),
             block = case_when(inside.range(amino_acid_position,c(as.numeric(box_block1_start)/3,as.numeric(box_block1_end)/3)) == T ~ 1,
                               inside.range(amino_acid_position,c(as.numeric(box_block2_start)/3,as.numeric(box_block2_end)/3)) == T ~ 2,
                               inside.range(amino_acid_position,c(as.numeric(box_block3_start)/3,as.numeric(box_block3_end)/3)) == T ~ 3))
    
    block1_positions <-which(good_aas_block_specifyer$block == 1)
    block2_positions <-which(good_aas_block_specifyer$block == 2)
    block3_positions <-which(good_aas_block_specifyer$block == 3)
    
    #_______________work well===generate Quadruple mutants
    
    #生成3突CDS——————————————————————————————
    make_quadruple_sequence <- function(df, codons_list) {
      out <- df %>% 
        pivot_longer(mut_1:mut_4, names_to="mut_num", values_to="mut") %>% 
        na.omit() %>% 
        tidyr::extract(mut, "(\\w)(\\d+)(\\w)", into=c("wt","position","mut")) %>%
        #拆分为三个部分w = word；d+means 数字
        #每次只读取一行，所以一次就是三个突变
        mutate(codons= yeast_favored_codons[mut,]) %>% #生成相应的密码子
        select(position, codons) %>% 
        rows_update(codons_list, ., by="position") %>% #注意position都需要class相同
        pull(codons) %>% str_c(collapse="")
      return(out)
    }
    
    pyr1_codons$position <- as.character(pyr1_codons$position)
    
    if(length(block1_positions) > 3 & input$qsm1 == T){
      quadruple_mut_number = 0
      for (i in block1_positions) { 
        for (j in block1_positions) {
          if (j>i) {
            for (m in block1_positions) {
              if (m>j) {
                for(o in block1_positions){
                  if(o>m){
                    for (k in good_aas[i,c(3:22)]) {
                      if (k=="")  {next}
                        for (l in good_aas[j,c(3:22)]) {
                          if (l=="") {next}
                            for (n in good_aas[m,c(3:22)]) {
                              if (n=="") {next}
                                for (p in good_aas[o,c(3:22)]) {
                                  if (p==""){next}
                                  block1_quadruple_number <- quadruple_mut_number + 1
                    }}}}}}}}}}}
      counter <-1
      block1_quadruple_mutants <- vector("list", block1_quadruple_number)# build a vector with many elements.
      for (i in block1_positions) { 
        for (j in block1_positions) {
          if (j>i) {
            for (m in block1_positions) {
              if (m>j) {
                for(o in block1_positions){
                  if(o>m){
                for (k in good_aas[i,c(3:22)]) {
                  if (k=="")  {next}
                  for (l in good_aas[j,c(3:22)]) {
                    if (l=="") {next}
                    for (n in good_aas[m,c(3:22)]) {
                      if (n=="") {next}
                      for (p in good_aas[o,c(3:22)]) {
                        if (p==""){next}
                      
                      block1_quadruple_mutants[[counter]]  <-  c(
                        mutant = counter,
                        mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                        mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                        mut_3  = paste(wt_aas[m], numbering[m], n, sep=""),
                        mut_4  = paste(wt_aas[o], numbering[o], p, sep="") 
                      )
                      counter <- counter + 1
                      
                      }}}}}}}}}}}
      block1_quadruple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant))%>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
      
      block1_quadruple_mutants$mutant<- as.character(block1_quadruple_mutants$mutant)
      
      output1 <- foreach(i=1:nrow(block1_quadruple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
        make_quadruple_sequence(block1_quadruple_mutants[i,], pyr1_codons)
      }
      
      if(length(output1)>1){
        quadruple_block1 <- cbind(block1_quadruple_mutants, output1[,1]) %>% tibble()
        colnames(quadruple_block1)[6] <- "mut_cds" #修改第5列的列名
      } else{
        quadruple_block1 <- cbind(block1_quadruple_mutants, output1) %>% tibble()
        colnames(quadruple_block1)[6] <- "mut_cds" #修改第5列的列名
      }
      
      
      quadruple_oligo_block1 <- quadruple_block1 %>% 
        mutate(block1= str_sub(mut_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) )) %>% #生成新列block1
        pivot_longer(block1,names_to="block", values_to="oligo")
      
    } else {quadruple_oligo_block1 <- "0"}
    
    
    if(length(block2_positions) > 3 & input$qsm2 == T){
      quadruple_mut_number = 0
      for (i in block2_positions) { 
        for (j in block2_positions) {
          if (j>i) {
            for (m in block2_positions) {
              if (m>j) {
                for(o in block2_positions){
                  if(o>m){
                for (k in good_aas[i,c(3:22)]) {
                  if (k=="")  {next}
                  for (l in good_aas[j,c(3:22)]) {
                    if (l=="") {next}
                    for (n in good_aas[m,c(3:22)]) {
                      if (n=="") {next}
                      for (p in good_aas[o,c(3:22)]) {
                        if (p==""){next}
                          block2_quadruple_number <- quadruple_mut_number + 1
                    }}}}}}}}}}}
      counter <-1
      block2_quadruple_mutants <- vector("list", block2_quadruple_number)# build a vector with many elements.
      for (i in block2_positions) { 
        for (j in block2_positions) {
          if (j>i) {
            for (m in block2_positions) {
              if (m>j) {
                for(o in block2_positions){
                  if(o>m){
                for (k in good_aas[i,c(3:22)]) {
                  if (k=="")  {next}
                  for (l in good_aas[j,c(3:22)]) {
                    if (l=="") {next}
                    for (n in good_aas[m,c(3:22)]) {
                      if (n=="") {next}
                      for (p in good_aas[o,c(3:22)]) {
                        if (p==""){next}
                      
                      block2_quadruple_mutants[[counter]]  <-  c(
                        mutant = counter,
                        mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                        mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                        mut_3  = paste(wt_aas[m], numbering[m], n, sep=""),
                        mut_4  = paste(wt_aas[o], numbering[o], p, sep="") 
                      )
                      counter <- counter + 1
                      
                    }}}}}}}}}}}
      block2_quadruple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
      
      block2_quadruple_mutants$mutant<- as.character(block2_quadruple_mutants$mutant)
      
      output2 <- foreach(i=1:nrow(block2_quadruple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
        make_quadruple_sequence(block2_quadruple_mutants[i,], pyr1_codons)
      }
      if(length(output2)>1){
        quadruple_block2 <- cbind(block2_quadruple_mutants, output2[,1]) %>% tibble()
        colnames(quadruple_block2)[6] <- "mut_cds" #修改第5列的列名
      } else{
        quadruple_block2 <- cbind(block2_quadruple_mutants, output2) %>% tibble()
        colnames(quadruple_block2)[6] <- "mut_cds" #修改第5列的列名
      }
      
      quadruple_oligo_block2 <- quadruple_block2 %>% 
        mutate(block2= str_sub(mut_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end)  )) %>%
        pivot_longer(block2,names_to="block", values_to="oligo")
      
    } else {quadruple_oligo_block2 <- "0"}
    
    if(length(block3_positions) > 3 & input$qsm3 == T){
      quadruple_mut_number = 0
      for (i in block3_positions) { 
        for (j in block3_positions) {
          if (j>i) {
            for (m in block3_positions) {
              if (m>j) {
                for(o in block3_positions){
                  if(o>m){
                for (k in good_aas[i,c(3:22)]) {
                  if (k=="")  {next}
                  for (l in good_aas[j,c(3:22)]) {
                    if (l=="") {next}
                    for (n in good_aas[m,c(3:22)]) {
                      if (n=="") {next}
                      for (p in good_aas[o,c(3:22)]) {
                        if (p==""){next}
                      block3_quadruple_number <- quadruple_mut_number + 1
                    }}}}}}}}}}}
      
      counter <-1
      block3_quadruple_mutants <- vector("list", block3_quadruple_number)# build a vector with many elements.
      for (i in block3_positions) { 
        for (j in block3_positions) {
          if (j>i) {
            for (m in block3_positions) {
              if (m>j) {
                for(o in block3_positions){
                  if(o>m){
                for (k in good_aas[i,c(3:22)]) {
                  if (k=="")  {next}
                  for (l in good_aas[j,c(3:22)]) {
                    if (l=="") {next}
                    for (n in good_aas[m,c(3:22)]) {
                      if (n=="") {next}
                      for (p in good_aas[o,c(3:22)]) {
                        if (p==""){next}
                      
                      block3_quadruple_mutants[[counter]]  <-  c(
                        mutant = counter,
                        mut_1  = paste(wt_aas[i], numbering[i], k, sep=""),
                        mut_2  = paste(wt_aas[j], numbering[j], l, sep=""),
                        mut_3  = paste(wt_aas[m], numbering[m], n, sep=""),
                        mut_4  = paste(wt_aas[o], numbering[o], p, sep="")
                      )
                      counter <- counter + 1
                      
                    }}}}}}}}}}}
      block3_quadruple_mutants %<>% bind_rows() %>% mutate(mutant=as.numeric(mutant)) %>% filter(str_detect(mut_1, "FALSE") == F) %>% filter(str_detect(mut_2, "FALSE") == F) %>% filter(str_detect(mut_3, "FALSE") == F)
      #here to modify the position column to as.numeric
      # create coding sequences for all of the quadruple mutants
      #把pyr1_codons以及block_quadruple_mutants中的position都改成as.character
      block3_quadruple_mutants$mutant<- as.character(block3_quadruple_mutants$mutant)
      
      output3 <- foreach(i=1:nrow(block3_quadruple_mutants), .packages=c("tidyverse"), .combine="rbind") %dopar% {
        make_quadruple_sequence(block3_quadruple_mutants[i,], pyr1_codons)
      }
      
      if(length(output3)>1){
        quadruple_block3 <- cbind(block3_quadruple_mutants, output3[,1]) %>% tibble()
        colnames(quadruple_block3)[6] <- "mut_cds" #修改第5列的列名
      } else{
        quadruple_block3 <- cbind(block3_quadruple_mutants, output3) %>% tibble()
        colnames(quadruple_block3)[6] <- "mut_cds" #修改第5列的列名
      }

      #先生成oligo再合并是最好的
      #生成3突oligos——————————————————————————————
      quadruple_oligo_block3 <- quadruple_block3 %>% 
        mutate(block3= str_sub(mut_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end)  )) %>%
        pivot_longer(block3,names_to="block", values_to="oligo")
      
    } else {quadruple_oligo_block3 <- "0"}
    
    if (length(quadruple_oligo_block1) >= 5 & length(quadruple_oligo_block2) >= 5 & length(quadruple_oligo_block3) >= 5){
      PYR1_lib_quadruple <- rbind(quadruple_oligo_block1, quadruple_oligo_block2,
                               quadruple_oligo_block3) %>%  # binds all blocks into one file
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(quadruple_oligo_block1) < 5 & length(quadruple_oligo_block2) >= 5 & length(quadruple_oligo_block3) >= 5){
      PYR1_lib_quadruple <- rbind(quadruple_oligo_block2,
                               quadruple_oligo_block3) %>%  # binds all blocks into one file
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(quadruple_oligo_block1) >= 5 & length(quadruple_oligo_block2) < 5 & length(quadruple_oligo_block3) >= 5){
      PYR1_lib_quadruple <- rbind(quadruple_oligo_block1,
                               quadruple_oligo_block3) %>%  # binds all blocks into one file
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(quadruple_oligo_block1) >= 5 & length(quadruple_oligo_block2) >= 5 & length(quadruple_oligo_block3) < 5){
      PYR1_lib_quadruple <- rbind(quadruple_oligo_block1,
                               quadruple_oligo_block2) %>%  # binds all blocks into one file
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    }  else if(length(quadruple_oligo_block1) < 5 & length(quadruple_oligo_block2) < 5 & length(quadruple_oligo_block3) >= 5){
      PYR1_lib_quadruple <- quadruple_oligo_block3 %>%
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(quadruple_oligo_block1) < 5 & length(quadruple_oligo_block2) >= 5 & length(quadruple_oligo_block3) < 5){
      PYR1_lib_quadruple <- quadruple_oligo_block2 %>%
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(quadruple_oligo_block1) >= 5 & length(quadruple_oligo_block2) < 5 & length(quadruple_oligo_block3) < 5){
      PYR1_lib_quadruple <- quadruple_oligo_block1 %>%
        select(block, "mutation_1"=mut_1, "mutation_2" = mut_2,"mutation_3"=mut_3,"mutation_4"=mut_4,oligo,"mutant_full_length_cds"=mut_cds)
    } else if(length(quadruple_oligo_block1) < 5 & length(quadruple_oligo_block2) < 5 & length(quadruple_oligo_block3) < 5){
      PYR1_lib_quadruple <- data.frame(warning = "No selected AA changes in any block")
    }
    
    
    return(PYR1_lib_quadruple)
  })
  
  #### Pooled Data Table ####
  combined_data_table_df<-reactive({
    # read in librarys
    lib_tab_1<-PYR1_lib_single_df()
    lib_tab_2<-PYR1_lib_double_df()
    lib_tab_3<-PYR1_lib_triple_df()
    lib_tab_4<-PYR1_lib_Quadruple_df()
    
    # series of if statements binding librarys if they exist 
    if(length(lib_tab_1) > 1 & length(lib_tab_2) > 1 & length(lib_tab_3) > 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_1,lib_tab_2,lib_tab_3,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>% 
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) > 1 & length(lib_tab_3) > 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_2,lib_tab_3,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) > 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_1,lib_tab_3,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) > 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_1,lib_tab_2,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) > 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_3,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) > 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_2,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_1,lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) > 1 & length(lib_tab_3) > 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_1,lib_tab_2,lib_tab_3) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) > 1 & length(lib_tab_3) > 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_2,lib_tab_3) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) > 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_1,lib_tab_3) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) > 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_1,lib_tab_2) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) > 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_3) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) > 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_2) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) > 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)<=1){
      combined_table<-rbind(lib_tab_1) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>% 
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)>1){
      combined_table<-rbind(lib_tab_4) %>% 
        mutate(`site_mutations_in_block` = case_when(is.na(mutation_2) == T ~ "SSM",
                                                     is.na(mutation_3) == T~ "DSM",
                                                     is.na(mutation_4) == T~ "TSM",
                                                     is.na(mutation_4) == F ~ "QSM")) %>% 
        mutate(oligo = case_when(block == "block1"~ str_sub(mutant_full_length_cds, as.numeric(box_block1_start) ,as.numeric(box_block1_end) ),
                                 block == "block2"~ str_sub(mutant_full_length_cds, as.numeric(box_block2_start) ,as.numeric(box_block2_end) ),
                                 block == "block3"~ str_sub(mutant_full_length_cds, as.numeric(box_block3_start) ,as.numeric(box_block3_end) ))) %>%
        select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
      
    }
    if(length(lib_tab_1) <= 1 & length(lib_tab_2) <= 1 & length(lib_tab_3) <= 1 & length(lib_tab_4)<=1){
      warning_table<-data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)")
      return(warning_table)
    }
    ###### Constitutive filter ########
    if (pyr1_constitutive_filter == T){
    combined_table_filtered <- combined_table %>%
      mutate(list_of_mutations = paste0(mutation_1,"_",mutation_2,"_",mutation_3,"_",mutation_4)) %>%
      filter(!str_detect(list_of_mutations,constitutive_combinations_regex)) %>% # known/predicted constitutive mutants for all double blocks, see prior to server for regex creation
      filter(!str_detect(list_of_mutations,"V83F|F159V")) %>%  # constitutive mutations for all block types so take them out
      filter(!str_detect(list_of_mutations,"I62")) # was never found in the original DSM library so take it out
    
    # example of how to filter,, will make the regex before the server and call for it here
    # combined_table_filtered <- combined_table %>% 
    #   mutate(list_of_mutations = paste0(mutation_1,"_",mutation_2,"_",mutation_3,"_",mutation_4)) %>% 
    #   filter(str_detect(list_of_mutations, "(?=.*K59A)(?=.*V81I)|(?=.*K59A)(?=.*V81R)")==F )
    return(combined_table_filtered)} else{
      return(combined_table)}
  })

  # renders table for combined data table output 
  output$SSM_DSM_TSM_table <- renderDT({
    combined_data_table <- combined_data_table_df()
    combined_data_table[,sapply(combined_data_table,class) == "logical"] <-sapply(combined_data_table[,sapply(combined_data_table,class) == "logical"],function(i) substr(as.character(i),1,1)) # fixes bug where F and T values are treated as False and True values rather than charecters
    table <- combined_data_table
    return(datatable(table,options = list(pageLength = 100),editable = F,class = 'cell-border stripe', rownames = F))
  })

  # download file for combined data table output
  output$SSM_DSM_TSM_library <- downloadHandler(
    filename = function() {
      paste("oligo_library",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(combined_data_table_df(), file,row.names = FALSE)
    })
  
  fasta_file_df<-reactive({
    fasta_file_df <- combined_data_table_df()
    fasta_file_df <- fasta_file_df %>% 
      mutate(fasta_name = paste0(">",site_mutations_in_block,":",block,":",mutation_1,":",mutation_2,":",mutation_3,":",mutation_4)) %>% 
      select(fasta_name, oligo) %>% 
      rowwise() %>% 
      reframe(fasta_name = c(fasta_name,oligo))
    return(fasta_file_df)
  })
  
  # download file for combined data table output
  output$SSM_DSM_TSM_library_fasta <- downloadHandler(
    filename = function() {
      paste0("oligo_library_fasta",Sys.Date(),".fasta")
    },
    content = function(file) {
      write_delim(fasta_file_df(), file,delim = "\n",col_names = F)
    })
  
  #### Getting summary information #####
  # counts how many unique oligos are in each block and puts this information into a table
  count_data_table_df<-reactive({
    combined_table<-combined_data_table_df()
    if(length(combined_table) == length(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))){
      return(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))
    }else {
    combined_table_summary <- combined_table %>% 
      group_by(`site_mutations_in_block`,block) %>% 
      reframe(site_mutations_in_block,block, mutation_count_per_block = n() ) %>% 
      distinct(site_mutations_in_block,block,.keep_all = T) %>% 
      mutate(ID = paste0(site_mutations_in_block,block)) %>% 
      select(ID, mutation_count_per_block)
    return(combined_table_summary)}
  })
  
  # renders table for basic oligo counts
  output$SSM_DSM_TSM_count_table <- renderDT({
    table <- count_data_table_df()
    sum_of_oligos<-sum(table$mutation_count_per_block,na.rm=T)
    sum_of_oligos_df<-tibble(ID = "Total oligos", number_of_oligos_per_block = sum_of_oligos)
    table <- table %>% 
      select(ID, "number_of_oligos_per_block"=mutation_count_per_block)
    new_table <- rbind(table,sum_of_oligos_df) %>% 
      select(`Block Type`=ID,`Count of Oligos`=number_of_oligos_per_block)
    return(datatable(new_table,options = list(pageLength = nrow(new_table)),editable = F,class = 'cell-border stripe', rownames = F))
  })
  
  # total oligo count text
  output$total_oligo_text <- renderText({
    table <- count_data_table_df()
    sum_of_oligos<-sum(table$mutation_count_per_block,na.rm=T)
    sum_of_oligos_df<-tibble(ID = "Total oligos", number_of_oligos_per_block = sum_of_oligos)
    table <- table %>% 
      select(ID, "number_of_oligos_per_block"=mutation_count_per_block)
    new_table <- rbind(table,sum_of_oligos_df) %>% 
      select(`Block Type`=ID,`Count of Oligos`=number_of_oligos_per_block)
    
    
    total_num<-max(new_table$`Count of Oligos`,na.rm=T)
    return(paste0("Total Oligos:"," ",total_num))
  })
  
  # adds download button for oligo counts
  output$SSM_DSM_TSM_count_download <- downloadHandler(
    filename = function() {
      paste("SSM_DSM_TSM_count",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(count_data_table_df(), file,row.names = FALSE)
    })
  
  # adds counts for all the mutation combinations possible
  mutation_site_table_df<-reactive({
    combined_table<-combined_data_table_df()
    if(length(combined_table) == length(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))){
      return(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))
    }else{
    
    count_summary_tabel<-count_data_table_df()
    summary_of_counts<-data.frame(ID = c("SSMblock1","SSMblock2", "SSMblock3", "DSMblock1","DSMblock2", "DSMblock3", "TSMblock1", "TSMblock2", "TSMblock3","QSMblock1", "QSMblock2", "QSMblock3"),mutation_count_per_block = 0 )
    count_summary_tabel <- rbind(count_summary_tabel,summary_of_counts) %>% 
      group_by(ID) %>% 
      filter(mutation_count_per_block == max(mutation_count_per_block)) # fixes bug where there was no value for types of blocks rather than there being 0 value for types of blocks
    
    # 1
    `100` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `010` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `001` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    
    # 2
    `110` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `101` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `011` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `200` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `020` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `002` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    
    # 3
    `111` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `120` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `102` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `210` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `012` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `201` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `021` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `300` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
    `030` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `003` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    
    # 4 
    `220` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `202` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `022` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `112` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `121` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `211` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `103` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `130` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `310` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `013` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `301` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `031` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `400` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
    `040` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
    `004` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
    
    # 5
    `113` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `131` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `311` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
    `221` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `212` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `122` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `203` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `230` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `320` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `023` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `302` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `032` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `104` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `140` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `410` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `014` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `401` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `041` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    
    # 6
    `222` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `330` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `303` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `033` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `123` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `213` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `132` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `231` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `312` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `321` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `204` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `240` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `420` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `024` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `402` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `042` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `114` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `141` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `411` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    
    # 7
    `133` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `313` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `331` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `322` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `232` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `223` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `304` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
    `340` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
    `430` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `034` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `403` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `043` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `124` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `214` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `142` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `241` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `412` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `421` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    
    # 8
    `233` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `323` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `332` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `422` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `242` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `224` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `440` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
    `404` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
    `044` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
    `431` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `413` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `341` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `143` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `314` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `134` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    
    # 9 
    `333` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `441` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
    `414` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
    `144` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
    `432` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `423` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `342` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `243` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    `324` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `234` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    
    # 10
    `334` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
    `343` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
    `433` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
    `442` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
    `424` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
    `244` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
    
    # 11
    `443` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
    `434` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
    `344` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
   
    # 12
    `444` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
    
    
    number_list <- c(`100`,`010`, `001`, # 1
                     `110`,`101`,`011`,`200`,`020`,`002`, # 2 
                     `111`,`120`,`102`,`210`,`012`,`201`,`021`,`300`,`030`,`003`, # 3 
                     `220`,`202`,`022`,`112`,`121`,`211`,`130`,`310`,`013`,`301`,`103`,`031`,`400`,`040`,`004`, # 4 
                     `221`,`212`,`122`,`113`,`131`,`311`,`203`,`023`,`230`,`032`,`302`,`320`,`140`,`104`,`410`,`014`,`401`,`041`, # 5 
                     `222`,`330`,`303`,`033`,`123`,`213`,`231`,`132`,`312`,`321`,`204`,`240`,`420`,`024`,`402`,`042`,`114`,`141`,`411`, # 6
                     `133`,`313`,`331`,`223`,`232`,`322`,`304`,`340`,`430`,`034`,`403`,`043`,`124`,`214`,`241`,`142`,`412`,`421`, # 7
                     `233`,`323`,`332`,`422`,`242`,`224`,`440`,`404`,`044`,`431`,`413`,`341`,`143`,`314`,`134`, # 8
                     `333`,`441`,`414`,`144`,`432`,`423`,`342`,`243`,`324`,`234`, # 9 
                     `433`,`343`,`334`,`442`,`424`,`244`, # 10
                     `443`,`434`,`344`, # 11
                     `444`) # 12
    
    mutation_site_1 <- sum(`100`,`010`, `001`) # 1
    mutation_site_2 <- sum(`110`,`101`,`011`,`200`,`020`,`002`) # 2
    mutation_site_3 <- sum(`111`,`120`,`102`,`210`,`012`,`201`,`021`,`300`,`030`,`003`) # 3
    mutation_site_4 <- sum(`220`,`202`,`022`,`112`,`121`,`211`,`130`,`310`,`013`,`301`,`103`,`031`,`400`,`040`,`004`) # 4
    mutation_site_5 <- sum(`221`,`212`,`122`,`113`,`131`,`311`,`203`,`023`,`230`,`032`,`302`,`320`,`140`,`104`,`410`,`014`,`401`,`041`) # 5
    mutation_site_6 <- sum(`222`,`330`,`303`,`033`,`123`,`213`,`231`,`132`,`312`,`321`,`204`,`240`,`420`,`024`,`402`,`042`,`114`,`141`,`411`) # 6
    mutation_site_7 <- sum(`133`,`313`,`331`,`223`,`232`,`322`,`304`,`340`,`430`,`034`,`403`,`043`,`124`,`214`,`241`,`142`,`412`,`421`) # 7
    mutation_site_8 <- sum(`233`,`323`,`332`,`422`,`242`,`224`,`440`,`404`,`044`,`431`,`413`,`341`,`143`,`314`,`134`) # 8
    mutation_site_9 <- sum(`333`,`441`,`414`,`144`,`432`,`423`,`342`,`243`,`324`,`234`) # 9
    mutation_site_10 <- sum(`433`,`343`,`334`,`442`,`424`,`244`) # 10
    mutation_site_11 <- sum(`443`,`434`,`344`) # 11
    mutation_site_12 <- sum(`444`) # 12
    
    total_mutation_combinations<-sum(mutation_site_1,mutation_site_2,mutation_site_3,mutation_site_4,mutation_site_5,mutation_site_6,mutation_site_7,mutation_site_8,mutation_site_9,mutation_site_10,mutation_site_11,mutation_site_12)
    
    summary_of_counts<-tibble(`Number of substituitons` = c("1","2", "3", "4","5", "6", "7", "8", "9","10","11","12","Total Potential Unique PYR1s"), # binds dataframe from cacluated data
                              `Count of unique PYR1s` =  prettyNum(c(mutation_site_1,mutation_site_2,mutation_site_3,mutation_site_4,mutation_site_5,mutation_site_6,mutation_site_7,mutation_site_8,mutation_site_9,mutation_site_10,mutation_site_11,mutation_site_12,total_mutation_combinations),big.mark = ",", scientific = FALSE) , 
                              count = c(mutation_site_1,mutation_site_2,mutation_site_3,mutation_site_4,mutation_site_5,mutation_site_6,mutation_site_7,mutation_site_8,mutation_site_9,mutation_site_10,mutation_site_11,mutation_site_12,total_mutation_combinations))
    
    summary_of_counts_filtered <- summary_of_counts %>% # removes 0s
      filter(count != 0) %>% 
      select(-count)
      
    
    return(summary_of_counts_filtered)}
  })
  
  # renders mutation combination table
  output$mutation_site_table <- renderDT({
    table <- mutation_site_table_df()
    return(datatable(table,options = list(pageLength = nrow(table)),editable = F,class = 'cell-border stripe', rownames = F))
  })
  
  # total PYR1 count
  output$total_pyr1_text <- renderText({
    table <- mutation_site_table_df()
    total_pyr1s <- table %>% 
      filter(`Number of substituitons`=="Total Potential Unique PYR1s") %>% 
      select(`Count of unique PYR1s`)
    
    total_num<-total_pyr1s[1,1]
    return(paste0("Total PYR1s:"," ",total_num))
  })
  
  
  # adds download button 
  output$mutation_site_download <- downloadHandler(
    filename = function() {
      paste("mutation_site_table_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(mutation_site_table_df(), file,row.names = FALSE)
    })
  
  # calculates mutation combination counts with oligo tables 
  golden_gate_table_df<-reactive({
    combined_table<-combined_data_table_df()
    if(length(combined_table) == length(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))){
      return(data.frame(warning = "No mutation sites checked/selected (SSM/DSM/TSM)"))
    }else if (input$long_format ==F){
      
      count_summary_tabel<-count_data_table_df()
      
      # 1
      `100` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `010` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `001` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      # 2
      `110` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `101` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `011` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `200` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `020` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `002` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      
      # 3
      `111` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `120` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `102` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `210` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `012` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `201` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `021` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `300` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
      `030` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `003` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      # 4 
      `220` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `202` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `022` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `112` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `121` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `211` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `103` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `130` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `310` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `013` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `301` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `031` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `400` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
      `040` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
      `004` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
      # 5
      `221` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `212` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `122` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `113` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `131` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `311` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
      `203` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `230` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `320` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `023` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `302` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `032` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `104` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `140` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `410` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `014` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `401` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `041` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      # 6
      `222` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `330` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `303` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `033` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `123` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `213` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `132` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `231` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `312` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `321` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `204` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `240` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `420` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `024` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `402` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `042` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `114` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `141` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `411` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      # 7
      `133` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `313` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `331` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `322` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `232` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `223` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `304` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
      `340` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
      `430` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `034` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `403` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `043` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `124` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `214` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `142` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `241` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `412` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `421` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      # 8
      `233` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `323` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `332` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `422` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `242` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `224` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `440` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
      `404` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
      `044` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
      `431` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `413` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `341` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `143` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `314` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `134` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      # 9 
      `333` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `441` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
      `414` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
      `144` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
      `432` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `423` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `342` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `243` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      `324` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `234` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      # 10
      `334` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
      `343` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
      `433` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
      `442` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
      `424` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
      `244` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
      # 11
      `443` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
      `434` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
      `344` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
      # 12
      `444` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
      
      number_list <- c(`100`,`010`, `001`, # 1
                       `110`,`101`,`011`,`200`,`020`,`002`, # 2
                       `111`,`120`,`102`,`210`,`012`,`201`,`021`,`300`,`030`,`003`, # 3
                       `220`,`202`,`022`,`112`,`121`,`211`,`130`,`310`,`013`,`301`,`103`,`031`,`400`,`040`,`004`, # 4
                       `221`,`212`,`122`,`113`,`131`,`311`,`203`,`023`,`230`,`032`,`302`,`320`,`140`,`104`,`410`,`014`,`401`,`041`, # 5
                       `222`,`330`,`303`,`033`,`123`,`213`,`231`,`132`,`312`,`321`,`204`,`240`,`420`,`024`,`402`,`042`,`114`,`141`,`411`, # 6
                       `133`,`313`,`331`,`223`,`232`,`322`,`304`,`340`,`430`,`034`,`403`,`043`,`124`,`214`,`241`,`142`,`412`,`421`, # 7
                       `233`,`323`,`332`,`422`,`242`,`224`,`440`,`404`,`044`,`431`,`413`,`341`,`143`,`314`,`134`, # 8
                       `333`,`441`,`414`,`144`,`432`,`423`,`342`,`243`,`324`,`234`, # 9
                       `334`,`343`,`433`,`442`,`424`,`244`, # 10
                       `443`,`434`,`344`, # 11
                       `444`) # 12

      mmm <- sum(`111`, # 3
               `112`,`121`,`211`, # 4
               `221`,`212`,`122`,`113`,`131`,`311`, # 5
               `222`,`123`,`213`,`231`,`132`,`312`,`321`,`114`,`141`,`411`, # 6
               `133`,`313`,`331`,`223`,`232`,`322`,`124`,`214`,`241`,`142`,`412`,`421`, # 7
               `233`,`323`,`332`,`422`,`242`,`224`,`431`,`413`,`341`,`143`,`314`,`134`, # 8
               `333`,`441`,`414`,`144`,`432`,`423`,`342`,`243`,`324`,`234`, # 9
               `442`,`424`,`244`,`433`,`343`,`334`, # 10
               `443`,`434`,`344`, # 11
               `444`) # 12
      mmw <- sum(`110`, # 2
               `120`,`210`, # 3
               `220`,`130`,`310`, # 4
               `230`,`320`,`140`,`410`, # 5
               `330`,`240`,`420`, # 6
               `340`,`430`, # 7
               `440`) # 8
      mwm <- sum(`101`, # 2
               `102`,`201`, # 3
               `202`,`103`,`301`, # 4
               `203`,`302`,`104`,`401`, # 5
               `303`,`204`,`402`, # 6
               `304`,`403`, # 7
               `404`)
      wmm <- sum(`011`, # 2
               `012`,`021`, # 3
               `022`,`013`,`031`, # 4
               `023`,`032`,`014`,`041`, # 5
               `033`,`024`,`042`, # 6
               `034`,`043`, # 7
               `044`) # 8

      mww <- sum(`100`, # 1
               `200`, # 2
               `300`, # 3
               `400`) # 4
      wmw <- sum(`010`, # 1
               `020`, # 2
               `030`, # 3
               `040`) # 4
      wwm <- sum(`001`, # 1
               `002`, # 2
               `003`, # 3
               `004`) # 4

      
      # 
      # mmm <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block3")], na.rm=T))
      # mmw <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block2")], na.rm=T))
      # mwm <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block3")], na.rm=T))
      # wmm <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block3")], na.rm=T))
      # mww <- sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block1")], na.rm=T)
      # wmw <- sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block2")], na.rm=T)
      # wwm <- sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"block3")], na.rm=T)

      
      summary_of_gg_counts<-tibble(block_combinations = c("MMM","MMW", "MWM", "WMM","MWW", "WMW", "WWM"),
                                   count_of_mutation_combinations = prettyNum(c(mmm,mmw,mwm,wmm,mww,wmw,wwm),big.mark = ",", scientific = FALSE))
      
      
      return(summary_of_gg_counts)} else if(input$long_format == T){
        count_summary_tabel<-count_data_table_df()
        # 1
        `100` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `010` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `001` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        # 2
        `110` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `101` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `011` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `200` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `020` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `002` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        
        # 3
        `111` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `120` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `102` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `210` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `012` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `201` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `021` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `300` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
        `030` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `003` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        # 4 
        `220` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `202` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `022` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `112` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `121` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `211` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `103` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `130` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `310` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `013` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `301` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `031` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `400` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
        `040` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
        `004` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
        # 5
        `221` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `212` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `122` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `113` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `131` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `311` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
        `203` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `230` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `320` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `023` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `302` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `032` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `104` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `140` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `410` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `014` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `401` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `041` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        # 6
        `222` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `330` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `303` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `033` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `123` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `213` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `132` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `231` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `312` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `321` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `204` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `240` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `420` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `024` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `402` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `042` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `114` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `141` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `411` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        # 7
        `133` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `313` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `331` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `322` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `232` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `223` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `304` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
        `340` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
        `430` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `034` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `403` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `043` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `124` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `214` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `142` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `241` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `412` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `421` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        # 8
        `233` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `323` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `332` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `422` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `242` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `224` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `440` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
        `404` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
        `044` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
        `431` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `413` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `341` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `143` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `314` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `134` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        # 9 
        `333` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `441` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock3")], na.rm=T))
        `414` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock2")], na.rm=T))
        `144` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"SSMblock1")], na.rm=T))
        `432` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `423` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `342` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `243` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        `324` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `234` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        # 10
        `334` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
        `343` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T))
        `433` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T))
        `442` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock3")], na.rm=T))
        `424` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock2")], na.rm=T))
        `244` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"DSMblock1")], na.rm=T))
        # 11
        `443` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock3")], na.rm=T))
        `434` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock2")], na.rm=T))
        `344` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"TSMblock1")], na.rm=T))
        # 12
        `444` <- prod(sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock1")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock2")], na.rm=T),sum(count_summary_tabel$mutation_count_per_block[str_detect(count_summary_tabel$ID,"QSMblock3")], na.rm=T))
        
        number_list <- c(`100`,`010`, `001`, # 1
                         `110`,`101`,`011`,`200`,`020`,`002`, # 2 
                         `111`,`120`,`102`,`210`,`012`,`201`,`021`,`300`,`030`,`003`, # 3 
                         `220`,`202`,`022`,`112`,`121`,`211`,`130`,`310`,`013`,`301`,`103`,`031`,`400`,`040`,`004`, # 4 
                         `221`,`212`,`122`,`113`,`131`,`311`,`203`,`023`,`230`,`032`,`302`,`320`,`140`,`104`,`410`,`014`,`401`,`041`, # 5 
                         `222`,`330`,`303`,`033`,`123`,`213`,`231`,`132`,`312`,`321`,`204`,`240`,`420`,`024`,`402`,`042`,`114`,`141`,`411`, # 6
                         `133`,`313`,`331`,`223`,`232`,`322`,`304`,`340`,`430`,`034`,`403`,`043`,`124`,`214`,`241`,`142`,`412`,`421`, # 7
                         `233`,`323`,`332`,`422`,`242`,`224`,`440`,`404`,`044`,`431`,`413`,`341`,`143`,`314`,`134`, # 8
                         `333`,`441`,`414`,`144`,`432`,`423`,`342`,`243`,`324`,`234`, # 9 
                         `433`,`343`,`334`,`442`,`424`,`244`, # 10
                         `443`,`434`,`344`, # 11
                         `444`) # 12
        total_mutations<-sum(number_list,na.rm = T)
        
        summary_of_gg_counts<-tibble(block_combinations = c("100","010", "001", # 1
                                                            "110","101","011","200","020","002", # 2 
                                                            "111","120","102","210","012","201","021","300","030","003", # 3
                                                            "220","202","022","112","121","211","130","310","103","013","301","031","400","040","004", # 4
                                                            "221","212","122","113","131","311","203","023","230","032","302","320","140","104","410","014","401","041", # 5
                                                            "222","330","303","033","123","213","231","132","312","321","204","240","420","024","402","042","114","141","411", # 6
                                                            "133","313","331","223","232","322","304","340","430","034","403","043","124","214","241","142","412","421", # 7
                                                            "233","323","332","422","242","224","440","404","044","431","413","341","143","314","134", # 8
                                                            "333","441","414","144","432","423","342","243","324","234", # 9
                                                            "442","424","244","433","343","334", # 10
                                                            "443","434","344", # 11
                                                            "444","total_mutations"),
                                     count_of_mutation_combinations = prettyNum(c(number_list,total_mutations),big.mark = ",", scientific = FALSE))
        
        
        return(summary_of_gg_counts)
        
      }
  })
  
  # adds table information depending on which check box used
  output$gg_table_message <-renderText({if(input$long_format == F)
    {paste0("Under `block_combinations` M = block synthesized in the oligo library, W = WT block based on native coding sequence entered in `Coding sequence used (DNA)` input. The order at which they appear is the order of block they represent, (ex: MMW = mutant block 1, mutant block 2, wild type block 3).")
  }else if (input$long_format == T){
    paste0("Under `block_combinations` 1 = SSM block, 2 = DSM block, 3 = TSM block, 4 = QSM block, 0 = wild type block based on native coding sequence entered in `Coding sequence used (DNA)` input. The order at which they appear is the order of block they represent, (ex: 120 = SSM block 1, DSM block 2, wild type block 3).")
    }
    })
  
  output$golden_gate_table <- renderDT({
    table <- golden_gate_table_df()
    return(datatable(table,options = list(pageLength = nrow(table)),editable = F,class = 'cell-border stripe', rownames = F))
  })
  
  # adds download button
  output$golden_gate_download <- downloadHandler(
    filename = function() {
      paste("golden_gate_table_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(golden_gate_table_df(), file,row.names = FALSE)
    })
  
  ##### Attaching barcodes to pooled data ########
  
  # allows the selection of sets
  # block 1
  output$barcode_block1_name <- renderUI({
    barcode_used<-input$block1_set 
    y<-barcoding_primers$set_name[barcode_used]
    return(textInput(inputId = "block1_barcode_name",
                       label = "Block 1 set name",
                       value = barcode_used))
  })
  output$barcode_block1_forward_primer <- renderUI({
    barcode_used<-input$block1_set 
    y<-barcoding_primers$F_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block1_barcode_F",
                     label = "Block 1 forward barcode primer 5'-3'",
                     value = y))
  })
  output$barcode_block1_reverse_primer <- renderUI({
    barcode_used<-input$block1_set 
    y<-barcoding_primers$R_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block1_barcode_R",
                     label = "Block 1 reverse barcode primer 5'-3'",
                     value = y))
  })
  # block 2
  output$barcode_block2_name <- renderUI({
    barcode_used<-input$block2_set 
    y<-barcoding_primers$set_name[barcode_used]
    return(textInput(inputId = "block2_barcode_name",
                     label = "Block 2 set name",
                     value = barcode_used))
  })
  output$barcode_block2_forward_primer <- renderUI({
    barcode_used<-input$block2_set 
    y<-barcoding_primers$F_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block2_barcode_F",
                     label = "Block 2 forward barcode primer 5'-3'",
                     value = y))
  })
  output$barcode_block2_reverse_primer <- renderUI({
    barcode_used<-input$block2_set 
    y<-barcoding_primers$R_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block2_barcode_R",
                     label = "Block 2 reverse barcode primer 5'-3'",
                     value = y))
  })
  # block 3
  output$barcode_block3_name <- renderUI({
    barcode_used<-input$block3_set 
    y<-barcoding_primers$set_name[barcode_used]
    return(textInput(inputId = "block3_barcode_name",
                     label = "Block 3 set name",
                     value = barcode_used))
  })
  output$barcode_block3_forward_primer <- renderUI({
    barcode_used<-input$block3_set 
    y<-barcoding_primers$F_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block3_barcode_F",
                     label = "Block 3 forward barcode primer 5'-3'",
                     value = y))
  })
  output$barcode_block3_reverse_primer <- renderUI({
    barcode_used<-input$block3_set 
    y<-barcoding_primers$R_seq[barcoding_primers$set_name ==barcode_used]
    return(textInput(inputId = "block3_barcode_R",
                     label = "Block 3 reverse barcode primer 5'-3'",
                     value = y))
  })
  
  ## generating barcoded oligos
  barcode_data_file<-reactive({
    # reading in inputs from `Attaching barcoding tags`
    b1f<-input$block1_barcode_F
    b1r<-as.character(seq_reverse(seq_complement(dna(input$block1_barcode_R)))) # reverse compliment
    b1name<-input$block1_barcode_name
    
    b2f<-input$block2_barcode_F
    b2r<-as.character(seq_reverse(seq_complement(dna(input$block2_barcode_R))))
    b2name<-input$block2_barcode_name
    
    b3f<-input$block3_barcode_F
    b3r<-as.character(seq_reverse(seq_complement(dna(input$block3_barcode_R))))
    b3name<-input$block3_barcode_name
    
    # reading pooled oligo table
    combined_table<-combined_data_table_df()
    combined_table<-combined_table %>% 
      select(`site_mutations_in_block`, block,mutation_1,mutation_2,mutation_3,mutation_4,oligo,mutant_full_length_cds)
    
    # adding BsaI sites to cleave off barcodes during golden gate reaction https://www.neb.com/en-us/products/r3733-bsai-hf-v2

    # # overhangs ATTC b1 CGAA # correct ones
    # # overhangs CGAA b2 AAGG
    # # overhangs AAGG b3 GGTC (after stop codon)
    
    
    restriction_site_forward <- "GGTCTCG" # GGTCTCN_your_sequence is needed
    bsmbI_forward <- "CGTCTCG" # CGTCTCN_your sequence
    bbsI_forward <- "GAAGACAC" # GAAGACNN_your sequence 
    
    restriction_site_forward<- reactive({
      if(input$restriction_site == "BsaI"){
        return("GGTCTCG")
      }else if (input$restriction_site == "BsmBI" ){
        return("CGTCTCG")
      }else if (input$restriction_site == "BbsI"){
        return("GAAGACAC")
      }else if (input$restriction_site == "None"){
        return("")
      }
    })
    
    # set up for attaching code blocks

    # adds flanking sites
    restriction_site_forward<-restriction_site_forward()
    bsaI_reverse_compliment <- as.character(seq_reverse(seq_complement(dna(restriction_site_forward)))) # reverse compliment
    combined_table_with_names <- combined_table %>% 
      mutate(oligo_barcode_primer_F = case_when(block == "block1" ~ b1f,
                                         block == "block2" ~ b2f,
                                         block == "block3" ~ b3f),
             oligo_barcode_primer_R = case_when(block == "block1" ~ input$block1_barcode_R,
                                         block == "block2" ~ input$block2_barcode_R,
                                         block == "block3" ~ input$block3_barcode_R),
             oligo_barcode_set = case_when(block == "block1" ~ b1name,
                                           block == "block2" ~ b2name,
                                           block == "block3" ~ b3name)) %>% # attaches barcode name and primer information
      mutate(oligo_translation = case_when(block == "block1"~ str_sub(seq_translate(dna(mutant_full_length_cds)), as.numeric(box_block1_start)/3,as.numeric(box_block1_end)/3),
                                           block == "block2"~ str_sub(seq_translate(dna(mutant_full_length_cds)), as.numeric(box_block2_start)/3,as.numeric(box_block2_end)/3),
                                           block == "block3"~ str_sub(seq_translate(dna(mutant_full_length_cds)), as.numeric(box_block3_start)/3,as.numeric(box_block3_end)/3))) %>%  # attaches a translation to double check the sequence information is correct
      mutate(oligo = case_when(block == "block1" ~ paste0(b1f,restriction_site_forward,oligo,bsaI_reverse_compliment,b1r ),
                               block == "block2" ~ paste0(b2f,restriction_site_forward,oligo,bsaI_reverse_compliment,b2r ),
                               block == "block3" ~ paste0(b3f,restriction_site_forward,oligo,bsaI_reverse_compliment,b3r )) # attaches barcode sequence to either side of the oligos
             ) %>%
      select(`site_mutations_in_block`, block,oligo_barcode_set,oligo_barcode_primer_F,oligo_barcode_primer_R,mutation_1,mutation_2,mutation_3,mutation_4,oligo_translation,oligo,mutant_full_length_cds)
    
    
    return(combined_table_with_names)
  })
  
  output$barcode_table <- renderDT({
    table <- barcode_data_file()
  })
  
  # adds download button
  # output$barcode_table_download <- downloadHandler(
  #   filename = function() {
  #     paste("barcoded_oligo_table_",Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(barcode_data_file(), file,row.names = FALSE)
  #   })
  
  # barcoded fasta file
  extract_numbers <- function(string) {
    str_extract(string, "\\d+") %>% str_replace("^0+", "")
  }
  barcoded_fasta_file_df<-reactive({
    library_name <-input$library_name
    pyr1_type<-input$pyr1_type
    fasta_file_df <- barcode_data_file()
    if(pyr1_type == "PYR1 WT"){
      fasta_file_df <- fasta_file_df %>%
        mutate(fasta_name = paste0(">","LIBRARY_",library_name,"_",
                                   "PYR1","_WT","_FRAG_",block,":",mutation_1,":",mutation_2,":",mutation_3,":",
                                   mutation_4,"_PRI_",str_extract(oligo_barcode_set,"set"),
                                   extract_numbers(oligo_barcode_set) )) %>%
        select(fasta_name, oligo) %>%
        rowwise() %>%
        reframe(fasta_name = c(fasta_name,oligo))
    }else if(pyr1_type == "PYR1*"){
      fasta_file_df <- fasta_file_df %>%
        mutate(fasta_name = paste0(">","LIBRARY_",library_name,"_",
                                   "PYR1","_STAR","_FRAG_",block,":",mutation_1,":",mutation_2,":",mutation_3,
                                   ":",mutation_4,"_PRI_",str_extract(oligo_barcode_set,"set"),extract_numbers(oligo_barcode_set) )) %>%
        select(fasta_name, oligo) %>%
        rowwise() %>%
        reframe(fasta_name = c(fasta_name,oligo))
    }else if(pyr1_type == "double_daggar"){
      fasta_file_df <- fasta_file_df %>%
        mutate(fasta_name = paste0(">","LIBRARY_",library_name,"_",
                                   "PYR1","_DOUBLE_DAGGAR","_FRAG_",block,":",mutation_1,":",mutation_2,":",
                                   mutation_3,":",mutation_4,"_PRI_",str_extract(oligo_barcode_set,"set",extract_numbers(oligo_barcode_set) ))) %>%
        select(fasta_name, oligo) %>%
        rowwise() %>%
        reframe(fasta_name = c(fasta_name,oligo))
    }else{
      fasta_file_df <- fasta_file_df %>%
        mutate(fasta_name = paste0(">","LIBRARY_",library_name,"_",
                                   "_FRAG_",block,":",mutation_1,":",mutation_2,":",mutation_3,":",mutation_4,
                                   "_PRI_",str_extract(oligo_barcode_set),extract_numbers(oligo_barcode_set) )) %>%
        select(fasta_name, oligo) %>%
        rowwise() %>%
        reframe(fasta_name = c(fasta_name,oligo))
    }
    barcoded_fasta_file_df<-fasta_file_df
    return(fasta_file_df)
  })

  
  # download file for combined data table output
  output$barcode_library_fasta <- downloadHandler(
    filename = function() {
      paste0("oligos_fasta_",input$library_name,".fasta")
    },
    content = function(file) {
      write_delim(barcoded_fasta_file_df(), file,delim = "\n",col_names = F)
    })
  
  # adding primer summary information
  primer_summary_file <- reactive({
    # Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    bc_df<-barcode_data_file()
    bc_primers <- bc_df %>% 
      select(block,oligo_barcode_set,oligo_barcode_primer_F,oligo_barcode_primer_R) %>% 
      distinct(block,.keep_all = T) %>%
      add_column(restriction_site = input$restriction_site) %>% 
      add_column(notes = "") 
      
      # group_by(oligo_barcode_set) %>% 
      # mutate(number_of_nuc_F = as.double(str_count(oligo_barcode_primer_F,"A|T|C|G"))) %>% 
      # mutate(number_of_gc_F = as.double(str_count(oligo_barcode_primer_F,"C|G"))) %>% 
      # mutate(number_of_nuc_R = as.double(str_count(oligo_barcode_primer_R,"A|T|C|G"))) %>% 
      # mutate(number_of_gc_R = as.double(str_count(oligo_barcode_primer_R,"C|G"))) %>% 
      # mutate(forward_primer_tm = (64.9+41*((number_of_gc_F)-16.4)/(number_of_nuc_F))) %>% 
      # mutate(reverse_primer_tm = (64.9+prod(41,((number_of_gc_R)-16.4))/(number_of_nuc_R))) %>% 
      # mutate(average_primer_tm = mean(c(forward_primer_tm,reverse_primer_tm),na.rm=T)) %>% 
      # mutate(test_count = str_count(oligo_barcode_primer_R,"G"))
  })
  
  ##### Oligo Summary Download #####
  
  output$oligo_summary_data_download <- downloadHandler(
    filename = function() {
      paste0("oligo_summary_data_",input$library_name,".xlsx")
    
  },
  
    content = function(file) {
      wb <- openxlsx::createWorkbook()
      
      amino_acid_table <- rv$data
      colnames(amino_acid_table) <-  c("WT AA","AA Position","A","C","D","E","F","G","H","I","K","L","M","N",
                                       "P", "Q", "R", "S","T","V","W","Y")
      openxlsx::addWorksheet(wb, "Sequence profile")
      openxlsx::writeData(wb, "Sequence profile", amino_acid_table)
      
      block_substituions_data<- tribble(
        ~Block, ~`1 Substitution`, ~`2 Substitutions`, ~`3 Substitutions`, ~`4 Substitutions`,
        "Block 1", as.character(input$ssm1), as.character(input$dsm1), as.character(input$tsm1), as.character(input$qsm1),
        "Block 2", as.character(input$ssm2), as.character(input$dsm2), as.character(input$tsm2), as.character(input$qsm2),
        "Block 3", as.character(input$ssm3), as.character(input$dsm3), as.character(input$tsm3), as.character(input$qsm3),
        "PYR1 Type Selected", input$pyr1_type, "", "", ""
      )
      
      openxlsx::addWorksheet(wb, "Block substituions selected")
      openxlsx::writeData(wb, "Block substituions selected", block_substituions_data)
      
      mutation_site_table <- mutation_site_table_df()
      openxlsx::addWorksheet(wb, "Unique PYR1s")
      openxlsx::writeData(wb, "Unique PYR1s", mutation_site_table)
      
      oligo_count_data <- count_data_table_df()
      openxlsx::addWorksheet(wb, "Oligo counts")
      openxlsx::writeData(wb, "Oligo counts", oligo_count_data)
      
      golden_gate_table_data <- golden_gate_table_df()
      openxlsx::addWorksheet(wb, "Cloning guide")
      openxlsx::writeData(wb, "Cloning guide", golden_gate_table_data)
      
      oligo_data <- barcode_data_file()
      openxlsx::addWorksheet(wb, "Oligos")
      openxlsx::writeData(wb, "Oligos", oligo_data)
      
      primer_data <- primer_summary_file()
      openxlsx::addWorksheet(wb, "Barcoding primers")
      openxlsx::writeData(wb, "Barcoding primers", primer_data)
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
   
  } # end conent
  
  
  ) # end download handeler
  
  
  
}

