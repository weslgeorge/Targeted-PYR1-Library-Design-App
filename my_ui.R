library(shiny)
library(shinyjs)
source("chooser.R")

my_ui <- fluidPage(
  theme = shinythemes::shinytheme("lumen"),  # <--- Specify theme here
  titlePanel("PAIR-D"),
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
                label = "Desired mutation matrix file CSV",
                accept = ".csv",
                placeholder = "mutation_table.csv"
                
      ),
      selectInput(inputId = "pyr1_type",label = "PYR1 type", choices = c("PYR1 WT","PYR1*")), # adds in pyr1_cds information
      # checkboxInput("pyr1_constitutive_filter", "PYR1 Constitutive mutations filter", TRUE), # uncomment to add option to not apply constitutive filter
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
                         h2("Getting started with the PAIR-Design oligo generator"),
                         h3("Description:"),
                         p("The goal of this app is to generate oligos that can be combined into a mutant PYR1 library using designed mutations."),
                         h3("Getting setup locally: WILL REMOVE FOR PUBLISHED VERSION"),
                         p("If you want to use the PAIR-D app locally on your computer simply..."),
                         tags$pre(
                           style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px; font-family: 'Courier New', monospace; font-size: 90%; white-space: pre-wrap; margin: 0;",
                           tags$code('some code to download dependancies')
                         ),
                         h3("PAIR-D Workflow:"),
                         p("The general workflow for designing to generate PYR1 mutant libraries is as follows. "),
                         div(
                           style = "display: flex; align-items: center; justify-content: space-between;",
                           
                           # Left side: header + paragraph stacked
                           div(
                             style = "flex: 1;",
                             h4("Mutation matrix generation:"),
                             p("First target small molecule interactors are identified, and similarly structured compounds are found in the PAIR data. This can be done with external chemical clustering methods, or by using programs like ",tags$a(href="https://chemminetools.ucr.edu/", "https://chemminetools.ucr.edu/")," to cluster a target compound with the PAIR data SMILES. Identified similar small molecules used to inform the library will be selected in the `Amino Acid Input Generation` tab to generate a mutation matrix that informs the app what mutations at what positions in PYR1 to allow in the library. Mutational matrices can be made manually and uploaded in the sidebar using the `Desired mutation matrix file CSV` upload button, as well as made in the `Amino Acid Input Generation` tab. See the `Mutation matrix` sheet of the `Oligo summary data` download button for a template of the file, new rows can be added as long as the proper WT amino acid and position are added. ")
                           ),
                           
                           # Right side: image
                           img(src = "subset-mutationmatrix.png", height="30%", width="30%")
                         ),
                         div(
                           style = "display: flex; align-items: center; justify-content: space-between;",
                           
                           # Left side: header + paragraph stacked
                           div(
                             style = "flex: 1;",
                             h4("Library parameter selection:"),
                             p("Next use the sidebar to inform the parameters for the library. Input a library name to help organize oligos in a larger Twist order, then select what PYR1 scaffold you wish to use and the number of mutations per block you wish to allow for. (note that increasing the number of mutations per block may make the app slow for a time while it generates oligos). Navigate to the `Oligo summary information` tab to see details on library size, and expected number of PYR1s generated by the library, and oligos needed for its construction. As a reference we will normally aim to make libraries between 10,000 and 500,000 unique PYR1s to have a realistically screenable library size. ")
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
                             p("The third step is to use the `Oligo Sequence and Barcoding Information` tab to select the barcodes and cut sites you want to use to flank each oligo block. This allows you to make multiple libraries and order them with a single oligo pool order from Twist biosciences, cutting down on library costs. Typically each block in a library gets a unique oligo barcode which we have defaulted to “set####” primers which we use in our own lab. These primer sequences will be in the downloadable report, and additionally can be browsed in the `pop_experiment.xlsx - 10K_primers.csv` file (Not all of these primers have been tested). Barcodes can also be entered manually if your lab has existing primers commonly used. Our library approach clones into the pBD-PYR1-153bp-F1 backbone using golden gate assembly; other cutsites can be selected but only BsaI cut sites will work with this backbone. ")
                           ),
                           
                           # Right side: image
                           img(src = "flanking-images.png", height="30%", width="30%")
                         ),
                         h4("Downloading files:"),
                         p("Finally make sure to download both the `Oligos.fasta` file and `Oligo summary data` file in the from the sidebar once all tabs have been opened and the library is to your liking. The oligos fasta file will generate a list of oligos that can be submitted to Twist bioscience orders and have the naming convention of >LIBRARY_”input library name”_PYR1_”PYR1 type”_FRAG_”Block number”_mutation1:mutation2:mutation3:mutation4_PRI_”barcode primers”. The `Oligo summary data` will download an excel file that contains sheets including the mutation matrix used to generate the library, oligo summary information, and a table of the oligos and the barcodes used."),
                         h3("Citations "),
                         p("We used R version 4.4.3 (R Core Team 2025) and the following R packages: bioseq v. 0.1.4 (Keck 2020), doParallel v. 1.0.17 (Corporation and Weston 2022), DT v. 0.33 (Xie, Cheng, and Tan 2024), openxlsx v. 4.2.8 (Schauberger and Walker 2025), shiny v. 1.10.0 (Chang et al. 2024), shinyjs v. 2.1.0 (Attali 2021), shinythemes v. 1.2.0 (Chang 2021), spatstat.utils v. 3.1.3 (Baddeley, Turner, and Rubak 2025), tidyverse v. 2.0.0 (Wickham et al. 2019)."),
                         p("Attali, Dean. 2021. shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. https://CRAN.R-project.org/package=shinyjs.
Baddeley, Adrian, Rolf Turner, and Ege Rubak. 2025. spatstat.utils: Utility Functions for “spatstat”. https://CRAN.R-project.org/package=spatstat.utils.
Chang, Winston. 2021. shinythemes: Themes for Shiny. https://CRAN.R-project.org/package=shinythemes.
Chang, Winston, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert, and Barbara Borges. 2024. shiny: Web Application Framework for r. https://CRAN.R-project.org/package=shiny.
Corporation, Microsoft, and Steve Weston. 2022. doParallel: Foreach Parallel Adaptor for the “parallel” Package. https://CRAN.R-project.org/package=doParallel.
Keck, Francois. 2020. “Handling Biological Sequences in r with the Bioseq Package.” Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13490.
R Core Team. 2025. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/.
Schauberger, Philipp, and Alexander Walker. 2025. openxlsx: Read, Write and Edit Xlsx Files. https://CRAN.R-project.org/package=openxlsx.
Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019. “Welcome to the tidyverse.” Journal of Open Source Software 4 (43): 1686. https://doi.org/10.21105/joss.01686.
Xie, Yihui, Joe Cheng, and Xianying Tan. 2024. DT: A Wrapper of the JavaScript Library “DataTables”. https://CRAN.R-project.org/package=DT.


Grateful was used to generate citations for this app
Rodriguez-Sanchez F, Jackson C (2024). _grateful: Facilitate citation of R packages_.
  <https://pakillo.github.io/grateful/>.
")
                         ),
                ###### AMINO ACID INPUT GENERATOR TAB #######
                tabPanel("Amino Acid Input Generation",
                         h3("Input Generation from PAIR data"),
                         verbatimTextOutput("selected_clusters"),
                         h4("Chemical Active Ingredient Filter"),
                         # Dropdown for multiple selections (with Selectize) for chemical choices
                         selectizeInput("chemname_choices", 
                                        "Select items:",
                                        choices = sort(unique(sensors_FDA_suppliment_df$act_ingr)),
                                        multiple = TRUE,
                                        options = list(placeholder = 'Select or remove items')),
                         
                         # Action button to collect the selections
                         actionButton("submit", "Submit"),
                         h3("Filtered Sensors Selected"),
                         DTOutput(outputId = "sensors_FDA_suppliment_df_filtered"),
                         h3("Mutation Matrix Generated"),
                         DTOutput(outputId = "reactive_values_dt")
                ),
                ####### SUMMARY INFORMATION TAB ########
                tabPanel("Oligo summary information",
                         h3("Unique PYR1s in Assembled Library"),
                         DTOutput(outputId = "mutation_site_table"),
                         
                         h3("Block oligo count table"),
                         DTOutput(outputId = "SSM_DSM_TSM_count_table"),
                        
                         
                         h3("Cloning Guide Table"),
                         DTOutput(outputId = "golden_gate_table"),
                         checkboxInput("long_format", "Long Format", FALSE),
                         textOutput(outputId = "gg_table_message"),
                         
                ),
                ###### ATTACHING BARCODING TAGS/Oligo data table ########
                tabPanel("Oligo Sequence and Barcoding Information",
                         h3("Flanking barcoding primer selection"),
                         h4("Block 1 Barcode"),
                         selectInput(inputId = "block1_set",label = "Pre-made set primer", choices = barcoding_primers$set_name),
                         tags$div(uiOutput("barcode_block1_name"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block1_forward_primer"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block1_reverse_primer"),style="display:inline-block"),
                         h4("Block 2 Barcode"),
                         selectInput(inputId = "block2_set",label = "Pre-made set primer", choices = barcoding_primers$set_name),
                         tags$div(uiOutput("barcode_block2_name"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block2_forward_primer"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block2_reverse_primer"),style="display:inline-block"),
                         h4("Block 3 Barcode"),
                         selectInput(inputId = "block3_set",label = "Pre-made set primer", choices = barcoding_primers$set_name),
                         tags$div(uiOutput("barcode_block3_name"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block3_forward_primer"),style="display:inline-block"),
                         tags$div(uiOutput("barcode_block3_reverse_primer"),style="display:inline-block"),
                         selectInput(inputId = "restriction_site",label = "Flanking Restriction Enzyme Site", choices = c("BsaI","BsmBI","BbsI","None")),
                         h3("Designed Oligos"),
                         DTOutput(outputId = "barcode_table") # add fasta file output oligo download button
                )
                
        
      ) # end tabPanel 
    ) # end main panel
    
  ) # end sidebar layout
  
  
  
) # end ui
