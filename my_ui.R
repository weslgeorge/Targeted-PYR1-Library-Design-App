library(shiny)
library(shinyjs)
source("data_files/chooser.R")

my_ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),  # <--- Specify theme here
  titlePanel("Targeted-PYR1-Library-Design-App"),
  h5(""),
  sidebarLayout(
    sidebarPanel( # use to control figures
      useShinyjs(),
      textInput(inputId = "library_name",
                label = "Library Name",
                value = ""
      ),
      # adds in good amino acids file
      fileInput(inputId = "good_aas",
                label = "Desired sequence profile CSV",
                accept = ".csv",
                placeholder = "custom_sequence_profile.csv"
                
      ),
      selectInput(inputId = "pyr1_type",label = "PYR1 type", choices = c("PYR1 WT","PYR1*","PYR1*mandi","PYR1doubledaggar","PYR1doubledaggar-mandi")), # adds in pyr1_cds information
      checkboxInput("pyr1_constitutive_filter", "PYR1 Constitutive mutations filter", TRUE), # uncomment to add option to not apply constitutive filter
      # h4("Block 1"), # uncomment to add option to change block positions
      # uiOutput("box_block1_start"), # numbox reactive in server # uncomment to add option to change block positions
      # uiOutput("box_block1_end"), # numbox reactive in server # uncomment to add option to change block positions
      h5("Number of substitutions allowed for Block 1"),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("ssm1", "1", TRUE)),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("dsm1", "2", TRUE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("tsm1", "3", FALSE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("qsm1", "4", FALSE)),
      # h4("Block 2"), # uncomment to add option to change block positions
      # uiOutput("box_block2_start"), # numbox reactive in server # uncomment to add option to change block positions
      # uiOutput("box_block2_end"), # numbox reactive in server # uncomment to add option to change block positions
      h5("Number of substitutions allowed for Block 2"),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("ssm2", "1", TRUE)),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("dsm2", "2", TRUE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("tsm2", "3", FALSE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("qsm2", "4", FALSE)),
      # h4("Block 3"), # uncomment to add option to change block positions
      # uiOutput("box_block3_start"), # numbox reactive in server # uncomment to add option to change block positions
      # uiOutput("box_block3_end"), # numbox reactive in server # uncomment to add option to change block positions
      h5("Number of substitutions allowed for Block 3"),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("ssm3", "1", TRUE)),
      div(style = "display: inline-block;vertical-align:top;",
          checkboxInput("dsm3", "2", TRUE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("tsm3", "3", FALSE)),
      div(style = "display: inline-block;vertical-align:bottom;",
          checkboxInput("qsm3", "4", FALSE)),
      
      h5("Quick Summary Information"),
      textOutput(outputId = "total_pyr1_text"),
      textOutput(outputId = "total_oligo_text"),
      
      h5("Data Download"),
      downloadButton("oligo_summary_data_download", "Oligo summary data"),
      downloadButton("barcode_library_fasta", "Oligos.fasta")
      
      
    ), # end sidebar panel # illustration
    
    ##### Main Panel #########
    mainPanel( # use to display figures 
      tabsetPanel( type = "tabs" ,
                   ##### READ ME ########
                   tabPanel("Readme",
                            h2("Getting started with the Targeted-PYR1-Library-Design-App oligo generator"),
                            h3("Description:"),
                            p("    The goal of this app is to generate oligos that can be combined into a mutant PYR1 library using designed mutations."),
                            h3("Targeted-PYR1-Library-Design-App Workflow:"),
                            p("    The general workflow for designing to generate PYR1 mutant libraries is as follows. A sequence profile can first be generated in the app or manually uploaded to the app. The sequence profile is just a set of substitution mutations allowed for oligo design, and wild type amino acids will always be included regardless of if they are specified in the profile. Following the upload or creation of the sequence profile, the app will break the PYR1 coding sequence into three distinct oligo blocks; block 1 being from nucleic acid position 153 to 306; block 2 being from nucleic acid position 303 to 449; and block 3 being from nucleic acid position 446 to 572. Due to this the sequence profile will only affect amino acids from amino acid position K59 to Q189 which is sufficient to mutate the PYR1 binding pocket. These oligo blocks are then appended with flanking barcoding sequences to allow for amplification, and internal BsaI restriction enzyme cut sites that allow cloning with the pBD-PYR1-156bp-F1 plasmid backbone to be used in yeast-2-hybrid chemical screens in combination with MAV99 harboring both pACT-HAB1 and the plasmid library constructed. "),
                            div(
                              style = "display: flex; align-items: center; justify-content: space-between;",
                              
                              # Left side: header + paragraph stacked
                              div(
                                style = "flex: 1;",
                                h4("Sequence profile generation:"),
                                p("    First target small molecule interactors are identified, and similarly structured compounds are found in the screening data. This can be done with external chemical clustering methods, or by using programs like ",tags$a(href="https://chemminetools.ucr.edu/", "https://chemminetools.ucr.edu/")," to cluster a target compound with the screening data SMILES. Identified similar small molecules used to inform the library will be selected in the `Sequence profile generation` tab to generate a sequence profile that informs the app what mutations at what positions in PYR1 to allow in the library. Selecting desired small molecules will filter for biosensor mutations that were identified in biosensors for the selected small molecules that then populate the sequence profile. Sequence profiles can be made manually and uploaded in the sidebar using the `Sequence profile CSV` upload button, as well as made in the `Sequence profile generation` tab. See the `Sequence profile` sheet of the `Oligo summary data` download button for a template of the file, new rows can be added as long as the proper WT amino acid and position are added. ")
                              ),
                              
                              # Right side: image
                              img(src = "subset-mutationmatrix.png",style= "margin-left: 20px;", height="30%", width="30%")
                            ),
                            div(
                              style = "display: flex; align-items: center; justify-content: space-between;",
                              
                              # Left side: header + paragraph stacked
                              div(
                                style = "flex: 1;",
                                h4("Library parameter selection:"),
                                p("    Next use the sidebar to inform the parameters for the library. Input a library name to help organize oligos in a larger Twist order, then select what PYR1 scaffold you wish to use and the number of mutations per block you wish to allow for. (note that increasing the number of mutations per block may make the app slow for a time while it generates oligos). The constitutive filter will filter out single and double mutation site blocks with known constitutive activity including two key residues that are commonly constitutive, and can help reduce library size. Navigate to the `Oligo summary information` tab to see details on library size, and expected number of PYR1s generated by the library, and oligos needed for its construction. As a reference we will normally aim to make libraries between 10,000 and 500,000 unique PYR1s to have a realistically screenable library size. ")
                              ),
                              
                              # Right side: image
                              img(src="subs-per-block.png", height="20%", width="20%"),
                              img(src = "block-params.png", height="30%", width="30%")
                            ),
                            div(
                              style = "display: flex; align-items: center; justify-content: space-between;",
                              
                              # Left side: header + paragraph stacked
                              div(
                                style = "flex: 1;",
                                h4("Barcoding oligo blocks:"),
                                p("    The third step is to use the `Oligo barcoding` tab to select the barcodes and cut sites you want to use to flank each oligo block. This allows you to make multiple libraries and order them with a single oligo pool order from Twist biosciences, cutting down on library costs. Typically each block in a library gets a unique oligo barcode which we have defaulted to “set####” primers which we use in our own lab. These primer sequences will be in the downloadable report, and additionally can be browsed in the `pop_experiment.xlsx - 10K_primers.csv` file (Not all of these primers have been tested). The primers are designed to be amplifeid using Q5 polymerase at an annealing temp of 60C. Barcodes can also be entered manually if your lab has existing primers commonly used. Our library approach clones into the pBD-PYR1-153bp-F1 backbone using golden gate assembly; other cutsites can be selected but only BsaI cut sites will work with this backbone.  ")
                              ),
                              
                              # Right side: image
                              img(src = "flanking-images.png", height="30%", width="30%")
                            ),
                            h4("Downloading files:"),
                            p("    Finally make sure to download both the `Oligos.fasta` file and `Oligo summary data` file in the from the sidebar once all tabs have been opened and the library is to your liking. The oligos fasta file will generate a list of oligos that can be submitted to Twist bioscience orders and have the naming convention of >LIBRARY_”input library name”_PYR1_”PYR1 type”_FRAG_”Block number”_mutation1:mutation2:mutation3:mutation4_PRI_”barcode primers”. The `Oligo summary data` will download an excel file that contains sheets including the sequence profile used to generate the library, oligo summary information, and a table of the oligos and the barcodes used. "),
                            h3("Citations "),
                            p("    We used R version 4.4.3 (R Core Team 2025) and the following R packages: bioseq v. 0.1.4 (Keck 2020), doParallel v. 1.0.17 (Corporation and Weston 2022), DT v. 0.33 (Xie, Cheng, and Tan 2024), openxlsx v. 4.2.8 (Schauberger and Walker 2025), shiny v. 1.10.0 (Chang et al. 2024), shinyjs v. 2.1.0 (Attali 2021), shinythemes v. 1.2.0 (Chang 2021), spatstat.utils v. 3.1.3 (Baddeley, Turner, and Rubak 2025), tidyverse v. 2.0.0 (Wickham et al. 2019)."),
                            p("Attali, Dean. 2021. shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. https://CRAN.R-project.org/package=shinyjs."),
                            p("Baddeley, Adrian, Rolf Turner, and Ege Rubak. 2025. spatstat.utils: Utility Functions for “spatstat”. https://CRAN.R-project.org/package=spatstat.utils."),
                            p("Chang, Winston. 2021. shinythemes: Themes for Shiny. https://CRAN.R-project.org/package=shinythemes."),
                            p("Chang, Winston, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert, and Barbara Borges. 2024. shiny: Web Application Framework for r. https://CRAN.R-project.org/package=shiny."),
                            p("Corporation, Microsoft, and Steve Weston. 2022. doParallel: Foreach Parallel Adaptor for the “parallel” Package. https://CRAN.R-project.org/package=doParallel."),
                            p("Keck, Francois. 2020. “Handling Biological Sequences in r with the Bioseq Package.” Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13490."),
                            p("R Core Team. 2025. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/."),
                            p("Schauberger, Philipp, and Alexander Walker. 2025. openxlsx: Read, Write and Edit Xlsx Files. https://CRAN.R-project.org/package=openxlsx."),
                            p("Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019. “Welcome to the tidyverse.” Journal of Open Source Software 4 (43): 1686. https://doi.org/10.21105/joss.01686."),
                            p("Xie, Yihui, Joe Cheng, and Xianying Tan. 2024. DT: A Wrapper of the JavaScript Library “DataTables”. https://CRAN.R-project.org/package=DT."),
                            
                            p("Grateful was used to generate citations for this app"),
                            p("Rodriguez-Sanchez F, Jackson C (2024). _grateful: Facilitate citation of R packages_.<https://pakillo.github.io/grateful/>.")
                   ),
                   ###### SEQUENCE PROFILE GENERATOR TAB #######
                   tabPanel("Sequence profile generation",
                            h3("Input generation from screen data"),
                            verbatimTextOutput("selected_clusters"),
                            p("    First make sure to have a list of desired chemicals, then use the dropdown or paste in a list of the desired chemicals to the input box to filter for biosensors binding to the small molecules listed and generating a Sequence profile for a new library."),
                            p("
                              "),
                            h4("Chemical Active Ingredient Filter"),
                            # Dropdown for multiple selections (with Selectize) for chemical choices
                            selectizeInput("chemname_choices", 
                                           "Select items:",
                                           choices = sort(unique(sensors_FDA_suppliment_df$act_ingr)),
                                           multiple = TRUE,
                                           options = list(placeholder = 'Select or remove items')),
                            
                            # Action button to collect the selections
                            actionButton("submit", "Submit"),
                            p("
                              "),
                            h4("Biosensors filtered from screening data"),
                            DTOutput(outputId = "sensors_FDA_suppliment_df_filtered"),
                            p("    The Biosensors filtered from screening data table shows the biosensor sequences and corresponding information filtered by the chemical(s) selected. The substitution mutations in the filtered biosensors are what inform the sequence profile below."),
                            p("
                              "),
                            h4("Sequence profile generated"),
                            DTOutput(outputId = "reactive_values_dt"),
                            p("    The Sequence profile generated table shown above details the amino acids substitutions allowed in the oligos generated. The table above is determined by the small molecules selected, and can be changed by interacting with the table above for quick changes to the library. **Changing the filter selection will remove manual changes to this table, it is recommended to make major manual changes on a csv file separately then upload it using the sidebar**. ")
                   ),
                   ####### SUMMARY INFORMATION TAB ########
                   tabPanel("Oligo summary information",
                            h3("Unique PYR1s in Assembled Library"),
                            DTOutput(outputId = "mutation_site_table"),
                            p("    The Unique PYR1s in Assembled Library table details the unique PYR1s that can be generated by the sequence profile and number of allowed substitutions per block parameters submitted. Use this to get information on how many PYR1 mutants can be expected in your library, and to play around with the sidebar parameters to see if you can get a large enough library to screen that is not so large it is untenable to screen. Typically we go for libraries between 10,000 unique PYR1s to 500,000. "),
                            p("
                              "),
                            h3("Block oligo count"),
                            DTOutput(outputId = "SSM_DSM_TSM_count_table"),
                            p("    Use the Block oligo count table to see how many oligos will be required for a library's construction. SSMblockX denotes an oligo with one designed mutation, DSMblockX denotes an oligo with two designed mutations, TSMblockX three mutations, QSMblockX four mutations.The more oligos ordered the more expensive the Twist order may become. Typical library sizes are between 100-400 oligos, but this varies greatly depending on how spread out the sequence profile is among the different oligo blocks and typically allowing more mutations per block will greatly increase the oligo count. "),
                            p("
                              "),
                            
                            h3("Cloning guide"),
                            DTOutput(outputId = "golden_gate_table"),
                            checkboxInput("long_format", "Long Format", FALSE),
                            textOutput(outputId = "gg_table_message"),
                            p("
                              "),
                            p("    Use the Cloning guide table to help give expected PYR1s for different block combinations. For targeted libraries we have typically pooled our SSM, DSM, TSM,...  blocks into the same oligo pools. This pool of mutant blocks is notated as M in this table while the position notest the block number, ie MMM is block 1 2 and 3 as pooled mutant blocks while MMW is blocks 1 and 2 as pooled mutant blocks and a WT PYR1 block 3. When we clone our libraries together we assemble with golden gate assembly, and for targeted libraries with pooled mutant blocks we will do seven golden gate reactions similar to the table depicted above. The following table will help to give an idea of the individual coverage goal for each golden gate assembly transformation, these transformations can be pooled, if they are simply add up the pooled assemblies. "),
                            p("    If you are ordering separate mutant oligo blocks on different barcodes by substitution and block, use the long format checkbox to make the guide show the long format cloning guide. This will split the assemblies so that they are divided by each block having a unique barcode to amplify, and can allow you to make what might be incredibly large libraries, smaller by restricting higher order mutation blocks (ie TSM, QSM blocks) to be on separate assemblies, and thus not create as many unique PYR1s.")
                            
                   ),
                   ###### ATTACHING BARCODING TAGS/Oligo data table ########
                   tabPanel("Oligo barcoding",
                            h3("Flanking barcoding primer selection"),
                            p("    Use this tab to select and create flanking sequence information to append to your oligo blocks. The flanking sequences can be used to uniquely amplify your oligos from the Twist oligo pool, and allow for ordering of multiple libraries at once, while still selecting what blocks to amplify out from the pool. Additionally this section will append BsaI restriction sites to the sequence so that the PYR1 blocks can be assembled into the pBD-PYR1-151bp-F1 backbone forming the unique PYR1s that make up your designed library. "),
                            p("    Simply use the dropdown to select what barcode to add to your oligo sequence, or manually enter a sequence and sequence name. If using the \"set\" primers they are designed to be amplifeid using Q5 polymerase at an annealing temp of 60C. "),
                            p("
                              "),
                            h4("Block 1 Barcode"),
                            selectInput(inputId = "block1_set",label = "Pre-made set primer", choices = barcoding_primers$name, selected = "ind_0001"),
                            tags$div(uiOutput("barcode_block1_name"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block1_forward_primer"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block1_reverse_primer"),style="display:inline-block"),
                            p("
                              "),
                            h4("Block 2 Barcode"),
                            selectInput(inputId = "block2_set",label = "Pre-made set primer", choices = barcoding_primers$name, selected = "ind_0002"),
                            tags$div(uiOutput("barcode_block2_name"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block2_forward_primer"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block2_reverse_primer"),style="display:inline-block"),
                            p("
                              "),
                            h4("Block 3 Barcode"),
                            selectInput(inputId = "block3_set",label = "Pre-made set primer", choices = barcoding_primers$name, selected = "ind_0003"),
                            tags$div(uiOutput("barcode_block3_name"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block3_forward_primer"),style="display:inline-block"),
                            tags$div(uiOutput("barcode_block3_reverse_primer"),style="display:inline-block"),
                            p("
                              "),
                            selectInput(inputId = "restriction_site",label = "Flanking Restriction Enzyme Site", choices = c("BsaI","BsmBI","BbsI","None")),
                            p("
                              "),
                            h3("Designed Oligos"),
                            DTOutput(outputId = "barcode_table") # add fasta file output oligo download button
                   )
                   
                   
      ) # end tabPanel 
    ) # end main panel
    
  ) # end sidebar layout
  
  
  
) # end ui
