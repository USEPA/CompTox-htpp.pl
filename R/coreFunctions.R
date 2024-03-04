#' Summarize, clean and format the raw data, and store it as a Mongo collection.  Run basic QC before other steps
#'
#' @param InputPath character string: Path to the Harmony file
#' @param PlateID character string: the plate_id for the well plate
#' @param mongoUrl character string: The URL of the mongo database holding the collection, with user credentials to access it
#' @param Cell_Type character string: the cell line used; null by default.
#' @param CellArea.Limit list: The size range for cells examined.  If the Cell_Type is "U-2 OS", "MCF7", "A549", "ARPE-19", "HepG2" or "HTB-9", this will autofill with their ranges.  Otherwise list(c(0, 99999999)) by default.
#' @param NucleiArea.Limit list: The nucleus area range for cells examined.  If the Cell_Type is "U-2 OS", "MCF7", "A549", "ARPE-19", "HepG2" or "HTB-9", this will autofill with their ranges.  Otherwise list(c(0,99999999)) by default.
#' @param n_max numeric: The maximum dimensions of the table
#' @param SType character string: The type of plate used, for instance, whether it is control
#'
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#' @import readr
#'
#' @export Raw2Level4

Raw2Level4 <- function(InputPath,  PlateID, mongoUrl, Cell_Type=NULL, CellArea.Limit=NULL, NucleiArea.Limit=NULL, n_max=2000*384,
                       SType = "vehicle control") {

  message(paste("InputPath: ", InputPath))

  ######### that information could also be stored in a collection or given as input to the function
  if(!is.null(Cell_Type)){
    if(Cell_Type=="U-2 OS"){
      CellArea.Limit <- c("U-2 OS" = list(c(100,6700)))
      NucleiArea.Limit <- c("U-2 OS" = list(c(20,900)))
    }else if(Cell_Type=="MCF7"){
      CellArea.Limit <- c("MCF7" = list(c(75,4500)))
      NucleiArea.Limit <- c("MCF7" = list(c(30,800)))
    }else if(Cell_Type=="A549"){
      CellArea.Limit <- c("A549" = list(c(75,10000)))
      NucleiArea.Limit <- c("A549" = list(c(30,800)))
    }else if(Cell_Type=="ARPE-19"){
      CellArea.Limit <- c("ARPE-19" = list(c(75,10000)))
      NucleiArea.Limit <- c("ARPE-19" = list(c(30,800)))
    }else if(Cell_Type=="ARPE-19"){
      CellArea.Limit <- c("ARPE-19" = list(c(150,13000)))
      NucleiArea.Limit <- c("ARPE-19" = list(c(35,900)))
    }else if(Cell_Type=="HepG2"){
      CellArea.Limit <- c("HepG2" = list(c(50,3700)))
      NucleiArea.Limit <- c("HepG2" = list(c(30,800)))
    }else if(Cell_Type=="HTB-9"){
      CellArea.Limit <- c("HTB-9" = list(c(100,4400)))
      NucleiArea.Limit <- c("HTB-9" = list(c(35,800)))
      #if there is only a CellArea.Limit and an input but not pre-coded cell type, set NucleiArea.Limit to default
    }else if(!is.null(CellArea.Limit) & is.null(NucleiArea.Limit)){
      NucleiArea.Limit <- list(c(0,99999999))
      #and vice versa, if there is only a NucleiArea.Limit and an input but not pre-coded cell type, set CellArea.Limit to default
    }else if(is.null(CellArea.Limit) & !is.null(NucleiArea.Limit)){
      CellArea.Limit <- list(c(0,99999999))
    }else{
      CellArea.Limit <- CellArea.Limit
      NucleiArea.Limit <- NucleiArea.Limit
    }
  }else{
    CellArea.Limit <- list(c(0,99999999))
    NucleiArea.Limit <- list(c(0,99999999))
  }



  ###################### 1. Read in the data and modify columns -> Table1 ##########################################
  message("***Step1: reading in the raw file and generating new columns")

  ### read in the Harmony file
  tic()
  Table1 <- suppressMessages( read_delim(paste0(InputPath, "/Objects_Population\ -\ Cells\ Non-Border.txt"), delim="\t", col_names=T, skip=9, n_max=n_max) )
  message(paste0("Table1 has ", dim(Table1)[1], " rows"))

  if(dim(Table1)[1]<n_max){

    cell_area_col <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Shape_Cells Area ")]
    nuclei_area_col <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Shape_Nuclei Area ")]

    cell_cen_x <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Position_Cells Centroid X in Image ")]
    cell_cen_y <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Position_Cells Centroid Y in Image ")]
    nuclei_cen_x <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Position_Nuclei Centroid X in Image ")]
    nuclei_cen_y <- colnames(Table1)[str_detect(colnames(Table1), "Cells Non-Border - Position_Nuclei Centroid Y in Image ")]

    ##generate 5 new columns
    Table1 <- as_tibble(Table1) %>% select(-starts_with("X1")) %>%
      dplyr::mutate(plate_id=PlateID,
                    well_id=Well3Digit(paste0(LETTERS[Row], Column)),
                    sample_id = paste0(plate_id, "_", well_id),
                    sample_id_field = paste0(sample_id, "_f", Field),
                    #flag cells with F (false) if the are too big or too small
                    cell_keep = ifelse(dplyr::between(get(cell_area_col),  CellArea.Limit[[1]][1],   CellArea.Limit[[1]][2]) &
                                         dplyr::between(get(nuclei_area_col), NucleiArea.Limit[[1]][1], NucleiArea.Limit[[1]][2]),
                                       T, F),
                    #use pythagoras to calculate a new column
                    Position_both_Centroid_Distance=sqrt((get(cell_cen_x) - get(nuclei_cen_x))^2+
                                                           (get(cell_cen_y) - get(nuclei_cen_y))^2))


    ###################### 2. Modify the column names to fit into mongo schema ##########################################
    message("***Step2: Modify the column names to fit into mongo schema")
    # obtain the list of feature names from the DB
    htpp_feature <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))
    FeatureList <- htpp_feature$find('{}')

    #find out which columns are needed
    GoodNames <- FeatureList %>% filter(column_type=="M" | !is.na(feature_name_mongo))
    table(GoodNames$column_type, useNA = "ifany")#10 metadata columns + 1300

    #retain only the columns that are needed
    Table1 <- Table1 %>% select(one_of(GoodNames$feature_name_harmony))

    #change all the column names from feature_name_harmony to feature_name_mongo
    for(iCol in 1:dim(Table1)[2]){
      oldName <- colnames(Table1)[iCol]
      newName <- GoodNames$feature_name_mongo[which(GoodNames$feature_name_harmony==oldName)]
      colnames(Table1)[iCol] <- newName
    }

    Table1 <- Table1 %>% select(sample_id_field, sample_id, plate_id, well_id, field, object_no, position_x, position_y, bounding_box, cell_keep, everything())  %>%
      dplyr::mutate_at(.vars=c("field", "object_no", "position_x", "position_y"), .funs="as.integer")

    ###################### 3. Calculate the median for each well --> htpp_well_raw ##########################################
    message("***Step3: Calculate the median for each well")
    FeatureNames = FeatureList$feature_name_mongo[which(FeatureList$column_type=="F" & !is.na(FeatureList$feature_id))]

    Summary1 <- Table1 %>%
      group_by(sample_id, plate_id, well_id) %>%
      dplyr::summarise(n_fields = as.integer(length(unique(sample_id_field))),
                       n_cells_total = as.integer(length(cell_keep)),
                       n_cells_keep = as.integer(sum(cell_keep)))

    Summary2 <- Table1 %>% dplyr::filter(cell_keep) %>%
      group_by(sample_id, plate_id, well_id) %>%
      dplyr::summarise_at(.vars=FeatureNames, .funs="median", na.rm=T)

    Summary <- Summary1 %>% left_join(Summary2)
    rm(Summary1, Summary2)

    cat("Inserting", dim(Summary)[1], "wells into htpp_well_raw for", PlateID, ".\n")

    htpp_well_trt <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))
    htpp_well_raw <- mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))

    for(SampleID in sort(unique(Summary$sample_id))){
      #if there is already a document with this sample_id, delete it
      if(length(htpp_well_raw$find(query=mongoQuery(sample_id = SampleID), fields='{"_id" : 1}')) > 0){
        cat("Entry ", SampleID, "already existed!!!\n")
        htpp_well_raw$remove(query=mongoQuery(sample_id = SampleID) )
      }
      #now the new document can be inserted
      htpp_well_raw$insert(Summary[which(Summary$sample_id==SampleID),], auto_unbox=TRUE, na = "null")
    }
    trtCount <- htpp_well_trt$count(query=mongoQuery(plate_id = PlateID, assay = "HTPP"))
    rawCount <- htpp_well_raw$count(query=mongoQuery(plate_id = PlateID))
    if(rawCount!=trtCount){
      warning(paste("Expected", trtCount, "documents in htpp_well_raw, based on htpp_well_trt, for", PlateID, "instead there are", rawCount, "documents."))
    }

    rm(FeatureList, GoodNames, Summary)
    message("***Step4: Perform MAD normalization to Level4")

    ############## 4.1 find out which solvent control wells should be used for the normalization
    htpp_well_trt <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))
    htpp_well_raw <- mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))

    ## 1.a) find which wells  are solvent controls
    #find all DMSO sample_id's from this plate, but only pick wells of good quality (qc_flag == "OK")
    DMSO_wells <- htpp_well_trt$find(query=mongoQuery(plate_id = PlateID, stype=SType,  qc_flag="OK"),
                                     fields='{"sample_id": 1, "_id":0}')$sample_id

    # should no longer be needed as binary qcflags are corrected to characters in SampleKey Validation.  Keeping as an atavism for now
    # if(length(DMSO_wells) == 0){ #logic for when qc_flags are binary
    #   DMSO_wells <- htpp_well_trt$find(query=mongoQuery(plate_id = PlateID, stype=SType, qc_flag = 0),
    #                                    fields='{"sample_id": 1, "_id":0}')$sample_id
    # }
    ## 1.b) find how many cells were in analysed in each well

    DMSO_cell_number <- htpp_well_raw$find(query=mongoQuery(sample_id = DMSO_wells),
                                           fields='{"sample_id": 1, "n_fields":1,  "n_cells_total":1, "n_cells_keep":1 ,"_id":0}')

    ## 1.c) calculate median across wells
    Median <- median(DMSO_cell_number$n_cells_keep, na.rm=T)
    DMSO_cell_number <- DMSO_cell_number %>%
      dplyr::mutate(well_keep = ifelse(between(n_cells_keep, 0.5*Median, 2*Median), T, F))

    ## 1.d) remove the bad wells from the list of DMSO wells
    DMSO_wells <- DMSO_cell_number$sample_id[which(DMSO_cell_number$well_keep)]

    ############## 4.2 calculate median + mad for solvent control cells: Attention! only use cells with good qc
    ## 2.a) get all cell data for DMSO wells
    DMSO_cells <- Table1 %>% filter(sample_id %in% DMSO_wells)

    ## 2.b) keep only cells that are considerd good
    DMSO_cells <- DMSO_cells %>% filter(cell_keep)

    ## 2.c) calculate median and nMad
    FeatureList <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))$find('{}') %>%
      filter(!is.na(feature_id))

    DMSO_raw_Median <- DMSO_cells %>%
      dplyr::summarise_at(.vars=FeatureList$feature_name_mongo, .funs="median",  na.rm=T)

    DMSO_raw_nMad <- DMSO_cells %>%
      dplyr::summarise_at(.vars=FeatureList$feature_name_mongo, .funs="mad", constant=1.4826, na.rm=T)

    #for features with nMad=0 or NA, set the nMad to 1 to not affect the scaling
    iCol <- which(DMSO_raw_nMad[1,]==0  | is.na(DMSO_raw_nMad[1,]))
    DMSO_raw_nMad[1,iCol]=1

    ############## 4.3 summarise for each well + normalize

    ## 3.a) find which wells need to be normalized
    WellList <- htpp_well_trt$find(query=mongoQuery(plate_id = PlateID, qc_flag = "OK"))

    # should no longer be needed as binary qcflags are corrected to characters in SampleKey Validation.  Keeping as an atavism for now
    # if(length(WellList) == 0){
    #   WellList <- htpp_well_trt$find(query=mongoQuery(plate_id = PlateID, qc_flag = 0),
    #                                  fields='{"qc_flag":0, "_id":0}')
    # }

    #add logic to handle cases where a whole plate has bad wells due to QC issue
    if(length(WellList) == 0){
      message(paste("No wells for plate", PlateID, "have qc_flag == 'OK', nothing will be written to htpp_well for this plate.\n", sep = " "))
    }else{

      WellList <- dplyr::arrange(WellList, sample_id)

      ## 3.b) load the raw medians from the collection
      Chem_raw_Median <- htpp_well_raw$find(query=mongoQuery(plate_id = PlateID, sample_id=WellList$sample_id))

      ## 3.c) calculate the relative cell number
      Chem_raw_Median <- Chem_raw_Median %>% dplyr::mutate(rel_cell_count = round(100*n_cells_keep/Median, 1))

      ## 3.d) attach the Median in the first row and the nMad in the second row
      Table_raw <- bind_rows(DMSO_raw_Median, DMSO_raw_nMad) %>% bind_rows(Chem_raw_Median) %>%
        select(ends_with("_id"), starts_with("n_"), rel_cell_count, starts_with("f_"))

      ## 3.e) normalize each feature for each well
      Table_norm <- Table_raw

      iRow <- 3:dim(Table_norm)[1]#the first two lines are the median and nMad and should not be normalized!

      for(FeatureName in FeatureList$feature_name_mongo){
        iCol <- which(colnames(Table_norm)==FeatureName)

        Table_norm[iRow,iCol] <- round( (Table_norm[iRow,iCol] - as.numeric(Table_norm[1,iCol])) / as.numeric(Table_norm[2,iCol]), 3)
      }

      #remove the two first rows that contained the median and nMad
      Chem_norm_Median <- Table_norm[iRow,]

      ############## 4.4 attach the SampleKey information from htpp_well_trt and write the documents into htpp_well

      ## 4.a) Attach the SampleKey information
      Table4 <- WellList %>% dplyr::left_join(Chem_norm_Median) %>% arrange(trt_name)

      ## 4.b) Write well-level normalized data into collection (takes ~ 1min)
      htpp_well <- mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))

      #if there is already a document with this plate_id, delete it
      if(length(htpp_well$find(query=mongoQuery(plate_id = PlateID), fields='{"_id" : 1}')) > 0){
        cat("Plate ", PlateID, "already existed!!!\n")
        htpp_well$remove(query=mongoQuery(plate_id = PlateID) )
      }

      for(SampleID in Table4$sample_id){
        newDocument <- Table4 %>% filter(sample_id==SampleID)
        htpp_well$insert(newDocument, auto_unbox=TRUE, na = "null")
      }#for each well

      wellCount <- htpp_well$count(query=mongoQuery(plate_id = PlateID))
      trtCount <- htpp_well_trt$count(query=mongoQuery(plate_id = PlateID, assay = "HTPP", qc_flag = "OK"))
      if(wellCount != trtCount){
        warning(paste("Expected", trtCount, "documents in well, based on htpp_well_trt, for", PlateID, "instead there are", wellCount, "documents."))
      }

      message(paste(htpp_well$count(query=mongoQuery(plate_id = PlateID)), "wells written to htpp_well"))
    }
    ###################### 5. read out some phenix metadata information and push it to a collection ##########################################
    message("***Step5: Phenix Metadata")
    htpp_image_metadata<-mongo(collection="htpp_image_metadata", url=mongoUrl, verbose=getOption("verbose"))

    Header <- suppressMessages( read_delim(paste0(InputPath, "Objects_Population - Cells Non-Border.txt"), delim="\t", col_names=F, skip=0, n_max=8) )
    Header <- Header %>% spread(key="X1", value="X2")

    imageMetadata <- suppressMessages( read_delim(paste0(InputPath, "../indexfile.txt"), delim="\t", col_names=T, skip=0) )
    imageMetadata <- as_tibble(imageMetadata) %>% select(-starts_with("X1")) %>%
      dplyr::mutate(well_id=Well3Digit(paste0(LETTERS[Row], Column))) %>%
      select(well_id, everything()) %>% select(-Row, -Column) %>%
      dplyr::mutate_at(.vars=c( "Plane", "Timepoint", "Field", "Channel ID"), .funs="as.integer")

    for(WellID in sort(unique(imageMetadata$well_id))){
      SampID <- paste0(PlateID, "_", WellID)
      imageTable <- imageMetadata %>% filter(well_id==WellID) %>% select(-well_id)
      newDocument <- list(sample_id  = SampID,
                          plate_id   = PlateID,
                          well_id    = WellID,
                          plate_name =      Header$"Plate Name" ,
                          measurement_nr =  Header$"Measurement",
                          evaluation_nr =   Header$"Evaluation",
                          evaluation_id =   Header$"Evaluation Signature",
                          input_path =      InputPath,
                          image_metadata = imageTable)

      #if there is already a document with this sample_id, delete it
      if(length(htpp_image_metadata$find(query=mongoQuery(sample_id = SampID), fields='{"_id" : 1}')) > 0){
        cat("Entry ", SampID, "already existed!!!\n")
        htpp_image_metadata$remove(query=mongoQuery(sample_id = SampID) )
      }
      #now the new document can be inserted
      htpp_image_metadata$insert(newDocument, auto_unbox=TRUE, na = "null")
    }#for each well
    imgMetaDataCount <- htpp_image_metadata$count(query=mongoQuery(plate_id = PlateID))
    trtCount <- htpp_well_trt$count(query=mongoQuery(plate_id = PlateID, assay = "HTPP"))
    if(imgMetaDataCount!= trtCount){
      warning(paste("Expected", trtCount, "documents in image_metadata, based on htpp_well_trt, for", PlateID, "instead there are", imgMetaDataCount, "documents."))
    }else{
      message("DONE")
    }



  }else{
    warning("maximum number of lines reached --> file not written to the database!!!")
  }

}#end of function

#' Level 5 analysis on plate data
#'
#' @param PlateGroup character string: The plate group id
#' @param SType character string: What type of plate is used, for QC check if the stype field in the data agrees with it
#' @param mongoUrl character string: The MongoDB url of the database with user credentials
#'
#' @return Median, nMAD and normalized well data
#'
#' @export Level5
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'

Level5 <- function(PlateGroup, SType = "vehicle control",  mongoUrl){

  ############## 1. grab all well-level data from the corresponding plate group

  htpp_well <- mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))
  Table4 = data.table(htpp_well$find(query=mongoQuery(pg_id = PlateGroup) ,
                                     fields='{"_id":0}'))

  message(paste(dim(Table4)[1], "Level4 observations retrieved"))
  #this should typically result in ca 384x4 obs = 1536 obs


  ############## 2. calculate median + mad for solvent control wells

  #Check if cell type is specified. If so, run normalization by cell type
  if("cell_type" %in% colnames(Table4)){

    for(cell in unique(Table4$cell_type)){

      message(paste("Beginning level 5 normalization for", cell, "cells", sep = " "))

      Table4_cell <- Table4[cell_type == cell,]

      ## 2.a) find solvent wells
      DMSO_wells <- Table4_cell %>% filter(stype==SType)

      ## 2.b) find names of features
      FeatureList <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))$find('{}') %>%
        filter(!is.na(feature_id))

      ## 2.c) calculate median and nMad
      DMSO_Mean <- DMSO_wells %>%
        dplyr::summarise_at(.vars=FeatureList$feature_name_mongo, .funs="mean", na.rm=TRUE)

      DMSO_SD <- DMSO_wells %>%
        dplyr::summarise_at(.vars=FeatureList$feature_name_mongo, .funs="sd", na.rm=TRUE)

      #for features with nMad=0 or NA, set the nMad to 1 to not affect the scaling
      iCol <- which(DMSO_SD[1,] == 0  | is.na(DMSO_SD[1,]))
      DMSO_SD[1,iCol] <- 1

      message(paste0("Features that could not be scaled in ", cell, ":"))
      message(colnames(DMSO_SD)[iCol])

      ############## 3. Normalize (Scale) each well
      ## 3.a) attach the Mean in the first row and the SD in the second row
      Table4_cell <- bind_rows(DMSO_Mean, DMSO_SD) %>% bind_rows(Table4_cell)

      # reorder the columns so that metadata is first
      Table4_cell <- Table4_cell %>% select(setdiff(colnames(Table4_cell), FeatureList$feature_name_mongo ), everything() )

      ## 3.b) normalize each feature for each well
      Table5 <- Table4_cell

      iRow <- 3:dim(Table5)[1]#the first two lines are the median and nMad and should not be normalized!

      for(FeatureName in FeatureList$feature_name_mongo){
        iCol <- which(colnames(Table5) == FeatureName)

        Table5[iRow,iCol] <- round((Table5[iRow,iCol] - Table5[1,iCol]) / Table5[2,iCol], 3)
      }

      #remove the two first rows that contained the median and nMad
      Table5 <- Table5[iRow,]

      ############## 4. Write into new collection
      htpp_well_norm <- mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))

      #if there are already a documents for this plate group, delete it
      if(length(htpp_well_norm$find(query = mongoQuery(pg_id = PlateGroup, cell_type = cell), fields='{"_id" : 1}')) > 0){
        cat("Entry for plate group", PlateGroup, "and", cell, "already existed!!!\n")
        htpp_well_norm$remove(query = mongoQuery(pg_id = PlateGroup, cell_type = cell))
      }

      for(SampleID in Table5$sample_id){
        newDocument <- Table5 %>% filter(sample_id == SampleID)
        htpp_well_norm$insert(newDocument, auto_unbox=TRUE, na = "null")
      }#for each well

      message(paste(dim(htpp_well_norm$find(query = mongoQuery(pg_id = PlateGroup), fields='{"_id" : 1}'))[1], "Level5 observations written for", cell, "cells"))
    }#for each cell type
  } else{

    #Standard normalization for studies with one cell type
    ## 2.a) find solvent wells
    DMSO_wells <- Table4 %>% filter(stype == SType)

    ## 2.b) find names of features
    FeatureList <- mongo(collection = "htpp_feature", url = mongoUrl, verbose = getOption("verbose"))$find('{}') %>%
      filter(!is.na(feature_id))

    ## 2.c) calculate median and nMad
    DMSO_Mean <- DMSO_wells %>%
      dplyr::summarise_at(.vars = FeatureList$feature_name_mongo, .funs = "mean", na.rm = TRUE)

    DMSO_SD <- DMSO_wells %>%
      dplyr::summarise_at(.vars = FeatureList$feature_name_mongo, .funs = "sd", na.rm = TRUE)

    #for features with nMad=0 or NA, set the nMad to 1 to not affect the scaling
    iCol <- which(DMSO_SD[1,] == 0  | is.na(DMSO_SD[1,]))
    DMSO_SD[1,iCol] <- 1

    message("Features that could not be scaled:")
    message(colnames(DMSO_SD)[iCol])

    ############## 3. Normalize (Scale) each well
    ## 3.a) attach the Mean in the first row and the SD in the second row
    Table4 <- bind_rows(DMSO_Mean, DMSO_SD) %>% bind_rows(Table4)

    # reorder the columns so that metadata is first
    Table4 <- Table4 %>% select(setdiff(colnames(Table4), FeatureList$feature_name_mongo ), everything() )

    ## 3.b) normalize each feature for each well
    Table5 <- Table4

    iRow <- 3:dim(Table5)[1]#the first two lines are the median and nMad and should not be normalized!

    for(FeatureName in FeatureList$feature_name_mongo){
      iCol <- which(colnames(Table5) == FeatureName)

      Table5[iRow,iCol] <- round((Table5[iRow,iCol] - Table5[1,iCol]) / Table5[2,iCol], 3)
    }

    #remove the two first rows that contained the median and nMad
    Table5 <- Table5[iRow,]

    ############## 4. Write into new collection
    htpp_well_norm <- mongo(collection = "htpp_well_norm", url = mongoUrl, verbose = getOption("verbose"))

    #if there are already a documents for this plate group, delete it
    if(length(htpp_well_norm$find(query=mongoQuery(pg_id = PlateGroup), fields='{"_id" : 1}')) > 0){
      cat("Entry for plate group", PlateGroup, "already existed!!!\n")
      htpp_well_norm$remove(query=mongoQuery(pg_id = PlateGroup) )
    }

    for(SampleID in Table5$sample_id){
      newDocument <- Table5 %>% filter(sample_id == SampleID)
      htpp_well_norm$insert(newDocument, auto_unbox = TRUE, na = "null")
    }#for each well

    message(paste(dim(htpp_well_norm$find(query = mongoQuery(pg_id = PlateGroup), fields='{"_id" : 1}'))[1], "Level5 observations written"))
  }
}

#' Reformats the data into a Mongo collection, normalizes it based on the solvent control and finds the percent responder cells.  Inputs the cv_well and cv_image_metadata collections.
#'
#' @param InputPath character string: the input path to the Harmony file
#' @param PlateID character string: the plate_id value
#' @param SType character string: Defines which sample type will be used for data normalization; "vehicle control" by default
#' @param mongoUrl character string: The MongoDB host, user, password and database
#' @param minNucleiArea numeric: The minimum area for something flagged as nucleus for QC
#' @param maxNucleiArea numeric: The maximum area for something flagged as nucleus for QC
#' @param minRoundness numeric: The minimum cell roundness for something to be recognized as a cell
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#'
#' @export CVanalysis
#'
#'
CVanalysis <- function(InputPath, PlateID, SType = "vehicle control",  mongoUrl,
                       minNucleiArea=30, maxNucleiArea=1000, minRoundness=0.5){

  message(paste("*************************", InputPath, "*************************"))

  ###################### 1. Read in the data and modify columns -> Table1 ##########################################
  message("***Step1: reading in the raw file and generate new columns")

  ### read in the Harmony file
  Table1 <- suppressMessages( read_delim(paste0(InputPath, "Objects_Population - Selected Nuclei.txt"), delim="\t", col_names=T, skip=9) )


  ##generate 5 new columns
  Table1<-as_tibble(Table1) %>% select(-starts_with("X5")) %>%
    dplyr::mutate(plate_id=PlateID,
                  well_id=Well3Digit(paste0(LETTERS[Row], Column)),
                  sample_id = paste0(plate_id, "_", well_id)) %>%
    #delete columns that are not needed
    dplyr::select(-Row, -Column, -Timepoint, -Plane, -10, -11, -`Selected Nuclei - ROI No`,
                  -Compound, -Concentration, -`Cell Type`, -`Cell Count`) #10 and 11 refer to columns 10 and 11


  ###################### 2. Modify the column names to fit into mongo schema ##########################################
  message("***Step2: Modify the column names to fit into mongo schema")

  ColNames <- colnames(Table1)

  ColNames <- str_replace_all(ColNames, "Selected Nuclei - ", "")
  ColNames <- str_replace_all(ColNames, "Selected Nuclei - ", "")
  ColNames <- str_replace_all(ColNames, "\\u00b5m\\u00B2", "um2")
  ColNames <- str_replace_all(ColNames, "\\u00b5m", "um")
  ColNames <- str_replace_all(ColNames, "\\u005B", "")
  ColNames <- str_replace_all(ColNames, "\\u005D", "")
  ColNames <- str_replace_all(ColNames, "%", "")

  ColNames <- str_trim(ColNames, side="both")

  ColNames <- str_replace_all(ColNames, " ", "_")
  ColNames <- str_to_lower(ColNames)

  colnames(Table1) <- ColNames

  Table2 <- Table1 %>% rename(position_x = x, position_y  = y) %>%
    select(sample_id, plate_id, well_id, field, object_no, position_x, position_y, bounding_box, everything())  %>%
    mutate_at(.vars=c("field", "object_no", "position_x", "position_y"), .funs="as.integer")

  ###################### 3. Normalize to solvent control ##########################################
  message("***Step3. Normalize to solvent control")


  ## 3.a) find which wells need to be normalized
  htpp_well_trt <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))

  if(htpp_well_trt$count(query=mongoQuery(plate_id = PlateID, assay = "CV"))<1){
    stop("There are no cell viability records in htpp_well_trt.")
  }

  WellList <- htpp_well_trt$find(query=mongoQuery(plate_id = PlateID, qc_flag = "OK"))%>%
    arrange(sample_id)

  ## 3.b) Select valid objects
  Cell.Data <- Table2 %>% dplyr::filter(sample_id %in% WellList$sample_id) %>%
    dplyr::mutate(cell_keep = (nucleus_area_um2<maxNucleiArea & nucleus_area_um2>minNucleiArea & nucleus_roundness>minRoundness))

  ## 3.c) Identify which solvent control wells can be used for the normalization

  # i) find which wells  are solvent controls
  DMSO_wells <- WellList$sample_id[which(WellList$stype==SType)]

  # ii) find how many cells were in analysed in each well
  DMSO_cell_number <- Cell.Data %>% filter(sample_id %in% DMSO_wells) %>%
    group_by(sample_id) %>% summarise(n_cells_keep = sum(cell_keep, na.rm=T)) %>% ungroup()

  # iii) calculate median across wells
  Median <- median(DMSO_cell_number$n_cells_keep, na.rm=T)
  DMSO_cell_number <- DMSO_cell_number %>%
    dplyr::mutate(well_keep = ifelse(between(n_cells_keep, 0.5*Median, 2*Median), T, F))

  # iv) remove the bad wells from the list of DMSO wells
  DMSO_wells <- DMSO_cell_number$sample_id[which(DMSO_cell_number$well_keep)]

  rm(DMSO_cell_number)

  ## 3.d) Calculate the thresholds based on the solvent control wells
  #Threshold for PI
  Threshold <- Cell.Data %>% filter(sample_id %in% DMSO_wells) %>% filter(cell_keep) %>%
    summarise(PI = quantile(propidium_iodide_mean, probs=0.95, na.rm=T))

  cat("PI threshold: ", Threshold$PI, "\n")

  #Median cell count
  DMSO_cells <- Cell.Data %>% filter(sample_id %in% DMSO_wells) %>% filter(cell_keep) %>%
    group_by(sample_id) %>%
    summarise(n_cells_total = sum(cell_keep),
              n_fields = length(unique(field)) ) %>% ungroup() %>%
    dplyr::mutate(n_cells_keep_per_field = n_cells_total/n_fields) %>%
    summarise(Median_cell_count = median(n_cells_keep_per_field))

  cat("Median cell count per field (solvent control): ", DMSO_cells$Median_cell_count, "\n")

  ### 3.e) find the % responder cells
  Well.Data <- Cell.Data %>%
    group_by(sample_id) %>%
    dplyr::mutate(n_fields = length(unique(field)),
                  n_cells_total = length(cell_keep),
                  n_cells_keep = sum(cell_keep)) %>%
    filter(cell_keep) %>%
    dplyr::mutate(isResponder_PI = propidium_iodide_mean>Threshold$PI) %>%
    group_by(sample_id, n_fields, n_cells_total, n_cells_keep) %>%
    summarise(responder_count_pi = sum(isResponder_PI),
              percent_responder_pi = round(100*mean(isResponder_PI),   2))  %>%
    dplyr::mutate(n_cells_keep_per_field = n_cells_total/n_fields,
                  rel_cell_count = round(100*n_cells_keep_per_field/DMSO_cells$Median_cell_count, 2) ) %>%
    select(-n_cells_keep_per_field)

  ### 3.f) attach the SampleKey information from htpp_well_trt and write the documents into cv_well

  # i) Attach the SampleKey information
  Table <- WellList %>% left_join(Well.Data) %>% arrange(trt_name)

  # Subset = Table %>% filter(grepl("positive", stype)) %>% arrange(conc) %>% select(chem_id, conc, rel_cell_count, percent_responder_pi)

  ## ii) Write well-level normalized data into collection (takes ~ 1min)
  cv_well <- mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))

  #if there is already documents with this plate_id delete it
  if(length(cv_well$find(query=mongoQuery(plate_id = PlateID), fields='{"_id" : 1}')) > 0){
    cat("PlateID ", PlateID, "already existed!!!\n")
    cv_well$remove(query=mongoQuery(plate_id = PlateID) )
  }

  for(SampleID in Table$sample_id){
    newDocument = Table %>% filter(sample_id==SampleID)
    cv_well$insert(newDocument, auto_unbox=TRUE, na = "null")
  }#for each well

  rm(WellList, Cell.Data, DMSO_wells, Threshold, DMSO_cells, Well.Data, Table, newDocument)

  ###################### 4. read out some phenix metadata information and push it to a collection ##########################################
  message("***Step4: read out some phenix metadata information and push it to a collection")

  Header <-suppressMessages( read_delim(paste0(InputPath, "Objects_Population - Selected Nuclei.txt"), delim="\t", col_names=F, skip=0, n_max=8) )
  Header <- Header %>% spread(key="X1", value="X2")

  imageMetadata <-  suppressMessages( read_delim(paste0(InputPath, "../indexfile.txt"), delim="\t", col_names=T, skip=0) )
  imageMetadata <- as_tibble(imageMetadata) %>% select(-starts_with("X1")) %>%
    dplyr::mutate(well_id=Well3Digit(paste0(LETTERS[Row], Column))) %>%
    select(well_id, everything()) %>% select(-Row, -Column) %>%
    dplyr::mutate_at(.vars=c( "Plane", "Timepoint", "Field", "Channel ID"), .funs="as.integer")

  cv_image_metadata <- mongo(collection="cv_image_metadata", url=url, verbose=getOption("verbose"))

  if(length(cv_image_metadata$find(query=mongoQuery(plate_id = PlateID), fields='{"_id" : 1}')) > 0){
    cat("Entry ", PlateID, "already existed!!!")
    cv_image_metadata$remove(query=mongoQuery(plate_id = PlateID) )
  }

  for(WellID in sort(unique(imageMetadata$well_id))){
    SampleID = paste0(PlateID, "_", WellID)
    imageTable = imageMetadata %>% filter(well_id==WellID) %>% select(-well_id)
    newDocument = list(sample_id  = SampleID,
                       plate_id   = PlateID,
                       well_id    = WellID,
                       plate_name =      Header$"Plate Name" ,
                       measurement_nr =  Header$"Measurement",
                       evaluation_nr =   Header$"Evaluation",
                       evaluation_id =   Header$"Evaluation Signature",
                       input_path =      InputPath,
                       image_metadata = imageTable)


    #now the new document can be inserted
    cv_image_metadata$insert(newDocument, auto_unbox=TRUE, na = "null")
  }#for each well

  trtCount <- htpp_well_trt$count(query=mongoQuery(plate_id = PlateID, assay = "CV"))
  imgMetaCount <- cv_well$count(query=mongoQuery(plate_id = PlateID))
  if(imgMetaCount != trtCount){
    warning(paste("Expected", trtCount, "documents in cv_well, based on htpp_well_trt, for", PlateID, "instead there are", imgMetaCount, "documents."))
  }
  wellCount <- cv_image_metadata$count(query=mongoQuery(plate_id = PlateID))
  if(wellCount != trtCount){
    warning(paste("Expected", trtCount, "documents in cv_well, based on htpp_well_trt, for", PlateID, "instead there are", wellCount, "documents."))
  }

  message("plate done")
}#end of function

#' Calculate global Mahalanobis distances for a table
#'
#' @param Table1 table: A table of well data
#' @param coverVariance numeric: The known variance of the well data
#' @param minObjects numeric: Minimum number of cells for plate to pass QC filter
#' @param SType character string: Defines which sample type will be used for data normalization; "vehicle control" by default
#' @param url character string: The MongoDB host, user, password and database
#'
#' @return A list of cumulative proportion, rotation matrix and inverse covariance
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export globalMahalanobisDistances
#'
#'

globalMahalanobisDistances <- function(Table1, coverVariance, minObjects, SType = "vehicle control",  url){

  ############## 1. grab all well-level data

  #Table1 = mongo(collection="htpp_well_norm", url=url, verbose=getOption("verbose"))$find(fields='{"_id":0}')

  Table2 = Table1 %>% filter(n_cells_keep>=minObjects & rel_cell_count>50 & stype!="null")

  ############## 2. Calculate the Eigen features from the well-level data (~ 1 min)

  iDataCols = which(grepl("^f_", colnames(Table2)))
  n_features = length(iDataCols)
  print(n_features)

  PCA = prcomp(Table2[,iDataCols], center=F, scale.=F)
  RotationMatrix = PCA$rotation

  CumProportion = cumsum(PCA$sdev^2)/sum(PCA$sdev^2)

  PC = length(which(CumProportion<coverVariance))+1
  print(PC)

  #add warning and change PC to 2 if PC = 1
  if(PC == 1){
    cat("Warning, the number of PCs at coverVariance", coverVariance, "equals 1. Add one more PC so PC = 2\n")
    PC <- 2
  }

  ##############  3. Find the inverse of the covariance matrix
  ## 3.a) Model the data to get the covariance matrix
  print("3. Find the inverse of the covariance matrix")

  Model = lm(PCA$x[,1:PC] ~ 0+Table2$trt_name)

  ## 3.b)
  Cov = estVar(Model)

  ## 3.c) inverse
  invCov = solve(Cov)

  ##############  4. Calculate Mahalanobis for each well on a per plate basis

  print("4. Calculate Global Mahalanobis for each well on a per plate basis")

  ## 4.a) transform data from table1 to the eigenfeatures
  iDataCols     = which(grepl("^f_", colnames(Table1)))
  iMetadataCols = which(!grepl("^f_", colnames(Table1)))

  data = as.matrix(Table1[,iDataCols]) %*% PCA$rotation[,1:PC]
  transfTable1 = cbind(Table1[,iMetadataCols], data)
  rm(data)

  iEigenfeatureCol=length(iMetadataCols)+(1:PC)

  PlateID = transfTable1$plate_id[1]

  Result <- NULL
  for (PlateID in unique(transfTable1$plate_id)){
    print(PlateID)

    Subset = transfTable1 %>% filter(plate_id==PlateID)

    ctrMean = Subset %>% filter(stype==SType) %>%
      summarise_at(.vars=colnames(Subset)[iEigenfeatureCol], .funs="mean")

    Delta = Subset
    Delta[iEigenfeatureCol] = sweep(as.matrix(Subset[iEigenfeatureCol]), 2, as.matrix(ctrMean), "-")

    D = apply(Delta[iEigenfeatureCol], 1, function(x) sqrt((x %*% invCov %*% x)))

    Distance = Subset %>% select(1:length(iMetadataCols)) %>%
      dplyr::mutate(n_features = n_features,
                    n_pc = PC,
                    d = round(D,3))

    Result = bind_rows(Result, Distance)

    rm(Subset, ctrMean, Delta, Distance)

  }#for each plate

  ##############  5. Write results in a new collection
  print("5. Write results into collection")
  htpp_global_mah <- mongo(collection="htpp_global_mah", url=url, verbose=getOption("verbose"))


  ## 5.) write each row as a document in the collection
  print(dim(Result)[1])
  for(i in 1:dim(Result)[1]){
    newDocument = Result[i,]
    htpp_global_mah$insert(newDocument, auto_unbox=TRUE, na = "null")
  }#for each row

  Output = list(CumProportion = CumProportion,
                RotationMatrix = RotationMatrix,
                invCov = invCov)
  return(Output)
}#end function

#' Calculate category mahalanobis distances for pipeline data
#' @param Level5 table: A table of well data at "level 5" in the pipeline
#' @param FeatureList matrix: The features for each category
#' @param CategoryName character string: The category whose variances are being compared
#' @param coverVariance numeric: The known variance of the well data
#' @param minObjects numeric: The minimum number of expected objects
#' @param SType character string: Defines which sample type will be used for data normalization; "vehicle control" by default
#' @param mongoUrl character string: the database where the collections are be stored and the required credentials, generated by the mongoURL function
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#'
#' @return VarianceExplainedList a list of variances with category Mahalanobis distances calculated
#'
#' @export categoryMahalanobisDistances
#'

categoryMahalanobisDistances <- function(Level5, FeatureList, CategoryName, coverVariance, minObjects,  SType = "vehicle control",  mongoUrl){
  message(paste("********************", CategoryName, "********************"))

  ############## 1. Subset the data to only contain features from the respective category

  #remove data from wells with little cells
  FeatureList <- FeatureList %>%
    filter(category_name_r==CategoryName)

  Subset1 <- Level5 %>% select(which(!grepl("^f", colnames(Level5))), one_of(FeatureList$feature_name_mongo))
  rm(Level5)
  #remove data from wells with little cells
  Subset2 <- Subset1 %>% filter(n_cells_keep>=minObjects & stype!="null")

  ############## 2. Calculate the Eigen features from the well-level data (~ 1 min)
  iDataCols <- which(grepl("^f_", colnames(Subset2)))
  n_features <- length(iDataCols)
  message(n_features)

  PCA <- prcomp(Subset2[,iDataCols], center=F, scale.=F)

  CumProportion = cumsum(PCA$sdev^2)/sum(PCA$sdev^2)

  PC <- length(which(CumProportion<coverVariance))+1
  message(PC)

  VarianceExplainedList <- as_tibble(data.frame("CategoryName" = CategoryName, "n_features" = n_features, "n_pc" = PC, "pc_varianceExplained" = str_c(CumProportion, collapse="|")))

  if(PC == 1){
    warning(paste("The number of PCs at coverVariance", coverVariance, "equals 1. Add one more PC so PC = 2\n"))
    PC <- 2
  }#I write the true number of PC that cover x variance in the Output, that I write the actual number of PC that was used for the distance calculation in the database

  ##############  3. Find the inverse of the covariance matrix

  message("3. Find the inverse of the covariance matrix")

  ## 3.a) Model the data to get the covariance matrix
  Model <- lm(PCA$x[,1:PC] ~ 0+Subset2$trt_name)

  ## 3.b)
  Cov <- estVar(Model)

  ## 3.c) inverse
  invCov <- solve(Cov)
  rm(Subset2)

  ##############  4. Calculate Mahalanobis for each well on a per plate basis

  message("4. Calculate Category Mahalanobis for each well on a per plate basis")

  ## 4.a) transform data from table1 to the eigenfeatures
  iDataCols     <- which(grepl("^f_", colnames(Subset1)))
  iMetadataCols <- which(!grepl("^f_", colnames(Subset1)))

  data <- as.matrix(Subset1[,iDataCols]) %*% PCA$rotation[,1:PC]
  transfSubset1 <- cbind(Subset1[,iMetadataCols], data)
  rm(data)

  iEigenfeatureCol<-length(iMetadataCols)+(1:PC)

  PlateID <- transfSubset1$plate_id[1]

  Result <- NULL
  for (PlateID in unique(transfSubset1$plate_id)){

    SubSubset <- transfSubset1 %>% filter(plate_id==PlateID)

    ctrMean <- SubSubset %>% filter(stype==SType) %>%
      summarise_at(.vars=colnames(SubSubset)[iEigenfeatureCol], .funs="mean")

    Delta <- SubSubset
    Delta[iEigenfeatureCol] = sweep(as.matrix(SubSubset[iEigenfeatureCol]), 2, as.matrix(ctrMean), "-")

    D <- apply(Delta[iEigenfeatureCol], 1, function(x) sqrt((x %*% invCov %*% x)))

    Distance <- SubSubset %>% select(1:length(iMetadataCols)) %>%
      dplyr::mutate(category_name_r = CategoryName,
                    n_features = n_features,
                    n_pc = PC,
                    d = round(D,3))

    if(is.null(Result)){
      Result = Distance
    }else{
      Result = bind_rows(Result, Distance)
    }

    rm(SubSubset, ctrMean, Delta, Distance)

  }#for each plate

  #somehow Result is not always organized with metadata first
  Result <- Result %>% select(setdiff(colnames(Result), c("category_name_r", "n_features", "n_pc", "d")), everything())

  ##############  5. Write results into collection

  message("5. Write results into collection")


  htpp_cat_mah <- mongo(collection="htpp_cat_mah", url=mongoUrl, verbose=getOption("verbose"))

  newLine <- list(Result, VarianceExplainedList)

  return(newLine)
}#end function

