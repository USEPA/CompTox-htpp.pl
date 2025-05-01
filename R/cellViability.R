#' Create mongo collection for cell viability by well (cv_well) from well treated collection
#'
#' @param file_path character string: file path to the top level directory of cell viability Harmony files for an HTPP dataset (i.e., the directory above plate-level directories)
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing cv_well collection and reinsert; FALSE by default
#' @param use_db  boolean: Determines whether mongoDB will be used or not; default is TRUE
#' @param json_collection_path character: Full file path to where JSON collections will be stored
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
#' @export generate_cvWell
#'
generate_cvWell <- function(file_path, mongoUrl="", rerun=FALSE, use_db=T, json_collection_path=""){
  if(use_db==T){
  #------------------------------------------------------------------------------------#
  # 1. Read in cell-level data and normalize it --> cv_well (~ 1min/plate)
  #------------------------------------------------------------------------------------#


  List <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))$find()%>% filter(assay=="CV")
  List <- unique(List$plate_id)


  if(rerun == TRUE){
    mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))$drop()
  }


  for(PlateID in List){
    message(paste("*****************", PlateID, "***************************"))
    #check if there is a folder with the respective data
    FolderList = list.files(file_path) #note that the number of folder is more than the total number of plates in List

    #find the folder that matches the PlateID
    FolderName = grep(PlateID, FolderList, value=T)
    if(length(FolderName)==0){
      print("There is no folder for this plate")
    }else{
      #check if the raw file exists
      if(file.exists(paste0(file_path, FolderName, "/Evaluation2/Objects_Population - Selected Nuclei.txt"))){
        InputPath = paste0(file_path, FolderName, "/Evaluation2/")
      }else if(file.exists(paste0(file_path, FolderName, "/Evaluation1/Objects_Population - Selected Nuclei.txt"))) {
        InputPath = paste0(file_path, FolderName, "/Evaluation1/")
      }  #if the file exists
      try(CVanalysis(InputPath=InputPath, PlateID=PlateID, SType = "vehicle control",  mongoUrl=mongoUrl,
                     minNucleiArea=30, maxNucleiArea=1000, minRoundness=0.5) )
    }#if there was a folder for this plateID

  }#for each PlateID

  #------------------------------------------------------------------------------------#
  # 2. Sanity checks
  #------------------------------------------------------------------------------------#



  message(mongo(collection="cv_image_metadata", url=mongoUrl, verbose=getOption("verbose"))$count())

  #How many wells have a good qc flag?
  SampleKey = mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(assay="CV"))

  message(mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))$count())}else{
    #------------------------------------------------------------------------------------#
    # 1. Read in cell-level data and normalize it --> cv_well (~ 1min/plate)
    #------------------------------------------------------------------------------------#


    if(file.exists(paste(json_collection_path,"htpp_well_trt.JSON",sep="/"))){
      htpp_well_trt  <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_trt.JSON",sep="/")))
      htpp_well_trt[htpp_well_trt == "NA"] <- NA
      List<-htpp_well_trt[assay=="CV"]
      if(nrow(List)<1){
        stop("The there are no cell viability documents in htpp_well_trt. These are needed to generate htpp_cv_well.  Please ensure htpp_well_trt is filled correctly before proceeding.")
      }
    } else {
      stop("htpp_well_trt collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }
    List <- unique(List[plate_id])


    if(rerun == TRUE){
      cv_wellJSON <- toJSON(data.table())
      write(cv_wellJSON, file=paste(json_collection_path,"cv_well.JSON",sep="/"))
    }


    for(PlateID in List){
      message(paste("*****************", PlateID, "***************************"))
      #check if there is a folder with the respective data
      FolderList = list.files(file_path) #note that the number of folder is more than the total number of plates in List

      #find the folder that matches the PlateID
      FolderName = grep(PlateID, FolderList, value=T)
      if(length(FolderName)==0){
        print("There is no folder for this plate")
      }else{
        #check if the raw file exists
        if(file.exists(paste0(file_path, FolderName, "/Evaluation2/Objects_Population - Selected Nuclei.txt"))){
          InputPath = paste0(file_path, FolderName, "/Evaluation2/")
        }else if(file.exists(paste0(file_path, FolderName, "/Evaluation1/Objects_Population - Selected Nuclei.txt"))) {
          InputPath = paste0(file_path, FolderName, "/Evaluation1/")
        }  #if the file exists
        try(CVanalysis(InputPath=InputPath, PlateID=PlateID, SType = "vehicle control",  mongoUrl="",
                       minNucleiArea=30, maxNucleiArea=1000, minRoundness=0.5, json_collection_path=json_collection_path, use_db=F))
      }#if there was a folder for this plateID

    }#for each PlateID

    #------------------------------------------------------------------------------------#
    # 2. Sanity checks
    #------------------------------------------------------------------------------------#

    if(file.exists(paste(json_collection_path,"cv_image_metadata.JSON",sep="/"))){
      cv_image_metadata  <-data.table(fromJSON(txt=paste(json_collection_path,"cv_image_metadata.JSON",sep="/")))
      cv_image_metadata[cv_image_metadata == "NA"] <- NA
      message(nrow(cv_image_metadata))
    } else {
      stop("cv_image_metadata failed to generate.  Check cell viability documents in  htpp_well_trt")
    }

    if(file.exists(paste(json_collection_path,"cv_well.JSON",sep="/"))){
      cv_well <-data.table(fromJSON(txt=paste(json_collection_path,"cv_well.JSON",sep="/")))
      cv_well[cv_well == "NA"] <- NA
      message(nrow(cv_well))
    } else {
      stop("cv_well failed to generate.  Check cell viability documents in  htpp_well_trt")
    }
  }
}


#' Creates and populates cell viability bmc (cv_bmc) collection in mongo
#'
#' @param cell_viability boolean: if cell_viability = TRUE, retrieve both `relative_cell_count` and `percent_responder_pi` data from cv_tcpl,
#' otherwise, only retrieve `relative_cell_count` data
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun TRUE will drop existing collection and reinsert; FALSE by default
#' @param use_db boolean: Determines whether mongoDB will be used or not; default is TRUE
#' @param json_collection_path character: Full file path to where JSON collections will be stored
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
#' @export generate_cvBMC
#'

generate_cvBMC <- function(cell_viability, mongoUrl="", rerun=FALSE, use_db=TRUE, json_collection_path=""){




  if(use_db==TRUE){
    #---------------------------------------------------------------------------------------------#
    # 1. Check what is in cv_bmc and delete if rerun == TRUE
    #---------------------------------------------------------------------------------------------#

    cv_bmc <- mongo(collection="cv_bmc", url=mongoUrl, verbose=getOption("verbose"))

    if(rerun == TRUE){
    cv_bmc$drop()
  }

    if(cell_viability == TRUE){
      #---------------------------------------------------------------------------------------------#
      # 2. Get the CC and PI data and modify the tables with BMR and BMC info
      #---------------------------------------------------------------------------------------------#

      CCData <- mongo(collection = "cv_tcpl", url = mongoUrl, verbose = getOption("verbose"))$find(query=mongoQuery(endpoint = "rel_cell_count"))
      CCData <- CCData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytostatic_bmr = "EC50",
                      cytostatic_bmc = ifelse(top<0 & bmd<max_conc, signif(bmd,3), NA),
                      cytostatic_flag = ifelse(!is.na(cytostatic_bmc) & bmd<min_conc, T, F),
                      cytostatic_flag = ifelse(n_conc<4, NA, cytostatic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)

      PIData <- mongo(collection="cv_tcpl", url = mongoURL, verbose=getOption("verbose"))$find(query=mongoQuery(endpoint = "percent_responder_pi"))
      PIData <- PIData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytotoxic_bmr = "3 nMad",
                      cytotoxic_bmc = ifelse(top>0 & bmd<max_conc, signif(bmd,3), NA),
                      cytotoxic_flag = ifelse(!is.na(cytotoxic_bmc) & bmd<min_conc, T, F),
                      cytotoxic_flag = ifelse(n_conc<4, NA, cytotoxic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)

      #---------------------------------------------------------------------------------------------#
      # 3. Get the treatment info and summarise the CC data to retain only the no-observed-effect dose level
      #---------------------------------------------------------------------------------------------#

      Trt <- mongo(collection = "htpp_well_trt", url = mongoUrl, verbose = getOption("verbose"))$find()

      Trt <- Trt %>% filter(qc_flag == "OK" & stype != "vehicle control") %>% select(pg_id, stype, chem_id, cell_type, dose_level, conc, conc_unit, trt_name) %>% distinct() %>%
        arrange(pg_id, chem_id, dose_level)

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(CCData) %>% dplyr::filter(is.na(cytostatic_bmc) | conc < cytostatic_bmc | (cytostatic_flag & dose_level == 1))

      Summary_CC <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytostatic_bmr, cytostatic_bmc, cytostatic_flag) %>%
        dplyr::summarise(cytostatic_noec_dose_level = max(dose_level)) %>%
        dplyr:: mutate(cytostatic_noec_dose_level = ifelse(cytostatic_flag, 0, cytostatic_noec_dose_level)) %>% data.table()

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(PIData) %>% dplyr::filter(is.na(cytotoxic_bmc) | conc < cytotoxic_bmc | (cytotoxic_flag & dose_level == 1))

      Summary_PI <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytotoxic_bmr, cytotoxic_bmc, cytotoxic_flag) %>%
        dplyr::summarise(cytotoxic_noec_dose_level = max(dose_level)) %>%
        dplyr::mutate(cytotoxic_noec_dose_level = ifelse(cytotoxic_flag, 0, cytotoxic_noec_dose_level))

      #---------------------------------------------------------------------------------------------#
      # 4. Define overall CV NOEC from CC and PI and write to collection -- Do by cell type
      #---------------------------------------------------------------------------------------------#
      for(cell in unique(Summary_CC$cell_type)){

        Summary = dplyr::full_join(Summary_CC[cell_type == cell,], Summary_PI[cell_type == cell,]) %>% group_by(stype, pg_id, chem_id) %>%
          dplyr::mutate(cv_bmc = min(cytostatic_bmc, cytotoxic_bmc, na.rm=T),
                        cv_flag = any(cytostatic_flag, cytotoxic_flag),
                        cv_noec_dose_level = min(cytostatic_noec_dose_level, cytotoxic_noec_dose_level)) %>%
          dplyr::mutate(cv_bmc = ifelse(cv_bmc==Inf, NA, cv_bmc)) %>%
          ungroup() %>% arrange(stype, pg_id, chem_id)

        message(paste0(capture.output(table(Summary$stype, Summary$cv_flag, useNA="ifany")), collapse = "\n")) #check output for problems

        message(paste0(capture.output(table(Summary$stype, Summary$cv_noec_dose_level, useNA="ifany")), collapse = "\n")) #check output for problems

        cat("Adding", dim(Summary)[1], "documents to cv_bmc for", cell, "cells\n")
        for(i in 1:dim(Summary)[1]){
          cv_bmc$insert(Summary[i,], auto_unbox=TRUE, na = "null")
        }

      } #for each cell_type

    }else{

      #---------------------------------------------------------------------------------------------#
      # 2. Get the CC data and modify the tables with BMR and BMC info
      #---------------------------------------------------------------------------------------------#

      CCData <- mongo(collection = "cv_tcpl", url = mongoUrl, verbose = getOption("verbose"))$find(query=mongoQuery(endpoint = "rel_cell_count"))
      CCData <- CCData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytostatic_bmr = "EC50",
                      cytostatic_bmc = ifelse(top<0 & bmd<max_conc, signif(bmd,3), NA),
                      cytostatic_flag = ifelse(!is.na(cytostatic_bmc) & bmd<min_conc, T, F),
                      cytostatic_flag = ifelse(n_conc<4, NA, cytostatic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)


      #---------------------------------------------------------------------------------------------#
      # 3. Get the treatment info and summarise the CC data to retain only the no-observed-effect dose level
      #---------------------------------------------------------------------------------------------#

      Trt <- mongo(collection = "htpp_well_trt", url = mongoUrl, verbose = getOption("verbose"))$find()

      Trt <- Trt %>% filter(qc_flag == "OK" & stype != "vehicle control") %>% select(pg_id, stype, chem_id, cell_type, dose_level, conc, conc_unit, trt_name) %>% distinct() %>%
        arrange(pg_id, chem_id, dose_level)

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(CCData) %>% dplyr::filter(is.na(cytostatic_bmc) | conc<cytostatic_bmc | (cytostatic_flag & dose_level==1))

      Summary_CC <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytostatic_bmr, cytostatic_bmc, cytostatic_flag) %>%
        dplyr::summarise(cytostatic_noec_dose_level = max(dose_level)) %>%
        dplyr:: mutate(cytostatic_noec_dose_level = ifelse(cytostatic_flag, 0, cytostatic_noec_dose_level)) %>% data.table()

      #---------------------------------------------------------------------------------------------#
      # 4. Define overall CV NOEC from CC and write to collection -- Do by cell type
      #---------------------------------------------------------------------------------------------#
      for(cell in unique(Summary_CC$cell_type)){

        Summary <- Summary_CC[cell_type == cell,] %>%
          dplyr::mutate(cv_bmc = cytostatic_bmc,
                        cv_flag = cytostatic_flag,
                        cv_noec_dose_level = cytostatic_noec_dose_level) %>%
          ungroup() %>%arrange(stype, pg_id, chem_id)

        message(paste0(capture.output(table(Summary$stype, Summary$cv_flag, useNA="ifany")), collapse = "\n")) #check output for problems

        message(paste0(capture.output(table(Summary$stype, Summary$cv_noec_dose_level, useNA="ifany")), collapse = "\n")) #check output for problems


        cat("Adding", dim(Summary)[1], "documents to cv_bmc for", cell, "cells\n")
        for(i in 1:dim(Summary)[1]){
          cv_bmc$insert(Summary[i,], auto_unbox=TRUE, na = "null")
        }

      } #for each cell type
    }

  cv_bmc <- mongo(collection="cv_bmc", url =mongoUrl, verbose=getOption("verbose"))
  bmcCount <- cv_bmc$count()
  cv_tcpl <- mongo(collection="cv_tcpl", url =mongoUrl, verbose=getOption("verbose"))
  tcplCount <-cv_tcpl$count()

  if(bmcCount != tcplCount){
    warning(paste("Expected", tcplCount, "documents in cv_bmc, based on cv_tcpl, instead there are", bmcCount, "documents."))
  }
  } else {
    #---------------------------------------------------------------------------------------------#
    # 1. Check what is in cv_bmc and delete if rerun == TRUE
    #---------------------------------------------------------------------------------------------#

    #load cv_bmc into a datatable if the json exists and you are not rerunning the function. Otherwise make an empty datatable.
    if(file.exists(paste(json_collection_path,"cv_bmc.JSON",sep="/")) & rerun==FALSE){
      cv_bmc  <-data.table(fromJSON(txt=paste(json_collection_path,"cv_bmc.JSON",sep="/")))
      cv_bmc[cv_bmc == "NA"] <- NA
    } else {
      cv_bmc <- data.table()
    }


    if(cell_viability == TRUE){
      #---------------------------------------------------------------------------------------------#
      # 2. Get the CC and PI data and modify the tables with BMR and BMC info
      #---------------------------------------------------------------------------------------------#

      if(file.exists(paste(json_collection_path,"cv_tcpl.JSON",sep="/"))){
        cv_tcpl <-data.table(fromJSON(txt=paste(json_collection_path,"cv_tcpl.JSON",sep="/")))
        cv_tcpl[cv_tcpl == "NA"] <- NA
        CCData<-cv_tcpl[endpoint == "rel_cell_count"]
      } else {
        stop("cv_tcpl.JSON required to generate cv_bmc.JSON.  Check that generate_cvTcpl ran correctly and that json_collection_path is accurate.")
      }
      CCData <- CCData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytostatic_bmr = "EC50",
                      cytostatic_bmc = ifelse(top<0 & bmd<max_conc, signif(bmd,3), NA),
                      cytostatic_flag = ifelse(!is.na(cytostatic_bmc) & bmd<min_conc, T, F),
                      cytostatic_flag = ifelse(n_conc<4, NA, cytostatic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)

      PIData<-cv_tcpl[endpoint == "percent_responder_pi"]
      PIData <- PIData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytotoxic_bmr = "3 nMad",
                      cytotoxic_bmc = ifelse(top>0 & bmd<max_conc, signif(bmd,3), NA),
                      cytotoxic_flag = ifelse(!is.na(cytotoxic_bmc) & bmd<min_conc, T, F),
                      cytotoxic_flag = ifelse(n_conc<4, NA, cytotoxic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)

      #---------------------------------------------------------------------------------------------#
      # 3. Get the treatment info and summarise the CC data to retain only the no-observed-effect dose level
      #---------------------------------------------------------------------------------------------#

      if(file.exists(paste(json_collection_path,"cv_tcpl.JSON",sep="/"))){
        Trt <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_trt.JSON",sep="/")))
        Trt[Trt == "NA"] <- NA
      } else {
        stop("htpp_well_trt.JSON required to generate cv_bmc.JSON.  Check that the file exists.")
      }

      Trt <- Trt %>% filter(qc_flag == "OK" & stype != "vehicle control") %>% select(pg_id, stype, chem_id, cell_type, dose_level, conc, conc_unit, trt_name) %>% distinct() %>%
        arrange(pg_id, chem_id, dose_level)

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(CCData) %>% dplyr::filter(is.na(cytostatic_bmc) | conc < cytostatic_bmc | (cytostatic_flag & dose_level == 1))

      Summary_CC <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytostatic_bmr, cytostatic_bmc, cytostatic_flag) %>%
        dplyr::summarise(cytostatic_noec_dose_level = max(dose_level)) %>%
        dplyr:: mutate(cytostatic_noec_dose_level = ifelse(cytostatic_flag, 0, cytostatic_noec_dose_level)) %>% data.table()

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(PIData) %>% dplyr::filter(is.na(cytotoxic_bmc) | conc < cytotoxic_bmc | (cytotoxic_flag & dose_level == 1))

      Summary_PI <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytotoxic_bmr, cytotoxic_bmc, cytotoxic_flag) %>%
        dplyr::summarise(cytotoxic_noec_dose_level = max(dose_level)) %>%
        dplyr::mutate(cytotoxic_noec_dose_level = ifelse(cytotoxic_flag, 0, cytotoxic_noec_dose_level))

      #---------------------------------------------------------------------------------------------#
      # 4. Define overall CV NOEC from CC and PI and write to collection -- Do by cell type
      #---------------------------------------------------------------------------------------------#
      for(cell in unique(Summary_CC$cell_type)){

        Summary = dplyr::full_join(Summary_CC[cell_type == cell,], Summary_PI[cell_type == cell,]) %>% group_by(stype, pg_id, chem_id) %>%
          dplyr::mutate(cv_bmc = min(cytostatic_bmc, cytotoxic_bmc, na.rm=T),
                        cv_flag = any(cytostatic_flag, cytotoxic_flag),
                        cv_noec_dose_level = min(cytostatic_noec_dose_level, cytotoxic_noec_dose_level)) %>%
          dplyr::mutate(cv_bmc = ifelse(cv_bmc==Inf, NA, cv_bmc)) %>%
          ungroup() %>% arrange(stype, pg_id, chem_id)

        message(paste0(capture.output(table(Summary$stype, Summary$cv_flag, useNA="ifany")), collapse = "\n")) #check output for problems

        message(paste0(capture.output(table(Summary$stype, Summary$cv_noec_dose_level, useNA="ifany")), collapse = "\n")) #check output for problems

        cat("Adding", dim(Summary)[1], "documents to cv_bmc for", cell, "cells\n")
        for(i in 1:dim(Summary)[1]){
          cv_bmc<-rbind(cv_bmc, Summary[i,])
        }

      } #for each cell_type

    }else{

      #---------------------------------------------------------------------------------------------#
      # 2. Get the CC data and modify the tables with BMR and BMC info
      #---------------------------------------------------------------------------------------------#

      if(file.exists(paste(json_collection_path,"cv_tcpl.JSON",sep="/"))){
        cv_tcpl <-data.table(fromJSON(txt=paste(json_collection_path,"cv_tcpl.JSON",sep="/")))
        cv_tcpl[cv_tcpl == "NA"] <- NA
        CCData<-cv_tcpl[endpoint == "rel_cell_count"]
      } else {
        stop("cv_tcpl.JSON required to generate cv_bmc.JSON.  Check that generate_cvTcpl ran correctly and that json_collection_path is accurate.")
      }

      CCData <- CCData %>% select(pg_id, chem_id, cell_type, stype, min_conc, max_conc, n_conc, top, bmd) %>%
        dplyr::mutate(cytostatic_bmr = "EC50",
                      cytostatic_bmc = ifelse(top<0 & bmd<max_conc, signif(bmd,3), NA),
                      cytostatic_flag = ifelse(!is.na(cytostatic_bmc) & bmd<min_conc, T, F),
                      cytostatic_flag = ifelse(n_conc<4, NA, cytostatic_flag)) %>%
        select(-min_conc, -max_conc, -n_conc, -top, -bmd)


      #---------------------------------------------------------------------------------------------#
      # 3. Get the treatment info and summarise the CC data to retain only the no-observed-effect dose level
      #---------------------------------------------------------------------------------------------#

      if(file.exists(paste(json_collection_path,"cv_tcpl.JSON",sep="/"))){
        Trt <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_trt.JSON",sep="/")))
        Trt[Trt == "NA"] <- NA
      } else {
        stop("htpp_well_trt.JSON required to generate cv_bmc.JSON.  Check that the file exists.")
      }

      Trt <- Trt %>% filter(qc_flag == "OK" & stype != "vehicle control") %>% select(pg_id, stype, chem_id, cell_type, dose_level, conc, conc_unit, trt_name) %>% distinct() %>%
        arrange(pg_id, chem_id, dose_level)

      #retain only dose-levels that are below the BMC
      Data <- Trt %>% dplyr::left_join(CCData) %>% dplyr::filter(is.na(cytostatic_bmc) | conc<cytostatic_bmc | (cytostatic_flag & dose_level==1))

      Summary_CC <- Data %>% dplyr::group_by(pg_id, stype, chem_id, cell_type, cytostatic_bmr, cytostatic_bmc, cytostatic_flag) %>%
        dplyr::summarise(cytostatic_noec_dose_level = max(dose_level)) %>%
        dplyr:: mutate(cytostatic_noec_dose_level = ifelse(cytostatic_flag, 0, cytostatic_noec_dose_level)) %>% data.table()

      #---------------------------------------------------------------------------------------------#
      # 4. Define overall CV NOEC from CC and write to collection -- Do by cell type
      #---------------------------------------------------------------------------------------------#
      for(cell in unique(Summary_CC$cell_type)){

        Summary <- Summary_CC[cell_type == cell,] %>%
          dplyr::mutate(cv_bmc = cytostatic_bmc,
                        cv_flag = cytostatic_flag,
                        cv_noec_dose_level = cytostatic_noec_dose_level) %>%
          ungroup() %>%arrange(stype, pg_id, chem_id)

        message(paste0(capture.output(table(Summary$stype, Summary$cv_flag, useNA="ifany")), collapse = "\n")) #check output for problems

        message(paste0(capture.output(table(Summary$stype, Summary$cv_noec_dose_level, useNA="ifany")), collapse = "\n")) #check output for problems


        cat("Adding", dim(Summary)[1], "documents to cv_bmc for", cell, "cells\n")
        for(i in 1:dim(Summary)[1]){
          cv_bmc<-rbind(cv_bmc, Summary[i,])
        }

      } #for each cell type
    }


    #convert cv_bmc to JSON and write it to disc
    bmcJSON <- toJSON(cv_bmc, na = "string", digits = 8)
    write(bmcJSON, file=paste(json_collection_path,"cv_bmc.JSON",sep="/"))

    #check that counts make sense
    bmcCount <- dim(cv_bmc)[1]
    tcplCount <- dim(cv_tcpl)[1]
    if(bmcCount != tcplCount){
      warning(paste("Expected", tcplCount, "documents in cv_bmc, based on cv_tcpl, instead there are", bmcCount, "documents."))
    }
}
}


#' Creates cell viability collection cv_tcpl based on well and chem data
#'
#' @param cell_viability boolean: if cell_viability = TRUE, function expects the cv_well collection to be populated,
#' otherwise use relative cell counts from htpp_well collection
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing cv_tcpl collection and reinsert; FALSE by default
#' @param use_db boolean: Determines whether mongoDB will be used or not; default is TRUE
#' @param json_collection_path character: Full file path to where JSON collections will be stored
#' @param minObjects numeric: minimum number of objects used for PI filtering; default is 50
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
#' @import tcplfit2
#' @export generate_cvTcpl
#'
generate_cvTcpl<-function(cell_viability, mongoUrl="", rerun=FALSE, use_db=TRUE, json_collection_path="", minObjects = 50){
  if(use_db==TRUE){  #setup
    htpp_well <- mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))
    cv_tcpl <- mongo(collection="cv_tcpl", url=mongoUrl, verbose=getOption("verbose"))
    cv_well <- mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))

    if(cell_viability==TRUE){
      message("cell_viability = TRUE: will generate cv_tcpl collection from cell viability data from the cv_well collection")
      if(cv_well$count()<1){
        stop("The cv_well collection is empty but is required to create cv_tcpl. Please ensure cv_well is created before proceeding.")
      }
    } else {
      message("cell_viability = FALSE: will generate cv_tcpl collection using relative cell counts from the htpp_well collection")
      if(htpp_well$count()<1){
        stop("The htpp_well collection is empty but is required to create cv_tcpl if cell_viability is FALSE. Please ensure htpp_well is created before proceeding.")
      }
    }

    if(rerun == TRUE){
      cv_tcpl$drop()
    }


    #Check number of plate groups and chemicals
    metaData <- data.table(mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))$find())


    for(PG in as.character(unique(metaData[ ,pg_id]))){
      message(paste("*************", PG, "*************"))

      if(cell_viability == TRUE){
        CVData_all <- data.table(mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))

        #Use htpp_well_trt to add in dtxsid to CVData_all
        htpp_well_trt <- data.table(mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(pg_id = as.character(PG))))
        htpp_well_trt <- htpp_well_trt[,c("sample_id", "dtxsid"), with = FALSE]
        CVData_all <- CVData_all %>% dplyr::inner_join(htpp_well_trt, by = c("sample_id", "dtxsid"))

        #Use htpp_chem to add chem_name to CVData_all
        htpp_chem <- data.table(mongo(collection="htpp_chem", url=mongoURL, verbose=getOption("verbose"))$find())
        htpp_chem <- htpp_chem[, c("dtxsid", "chem_name")]
        CVData_all <- CVData_all %>% dplyr::inner_join(htpp_chem, by = c("dtxsid"))

      }else{
        CVData_all <- data.table(mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
      }

      for(cell in unique(CVData_all[, cell_type])){

        cat("Generating CV data for", cell, "cells\n")

        CVData <- CVData_all[cell_type == cell]

        CompoundList <- unique(CVData$chem_id[which(CVData$stype != "vehicle control")])

        if(cell_viability == TRUE){

          ############## fit cell count ##############
          message(paste("******", "cell count"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype == "vehicle control" & rel_cell_count >= 50 & n_fields == nCV)

          Control <- Control %>% summarise(Median = median(rel_cell_count, na.rm = T),
                                           nMad = mad(rel_cell_count, constant = 1.4826, na.rm = T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct()

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$rel_cell_count,
                          bmed=Control$Median,
                          cutoff=2*Control$nMad, #only fit curve if it exceeds 2*nMad;
                          onesd=50/1.349, #onesd is 50% divided by 1.349 to have a BMR of 50%
                          approach = "CV", endpoint = "rel_cell_count")

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = TRUE))

              if(is.null(newLine) | class(newLine) == "try-error"){
                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine) == "try-error"){
                newLine <- row
              }
              rownames(newLine) <- ""

              if(length(cv_tcpl$find(query = mongoQuery(pg_id = PG, chem_id = Chem, cell_type = cell, stype = sample_type, approach = "CV", endpoint = "rel_cell_count"), fields = '{"_id" : 1}')) > 0){
                cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                cv_tcpl$remove(query = mongoQuery(pg_id = PG, chem_id = Chem, stype = sample_type, approach = "CV", endpoint = "rel_cell_count") )
              }
              cv_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")
            } #for each stype

          }#for each Chem

          ############## fit  PI response ##############
          message(paste("******", "PI response"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype=="vehicle control" & n_cells_keep >= minObjects & n_fields == nCV)

          Control <- Control%>% summarise(Median = median(percent_responder_pi, na.rm=T),
                                          nMad   = mad(percent_responder_pi, constant=1.4826, na.rm=T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem & n_cells_keep >= minObjects) %>% filter(stype == sample_type)

              Metadata <- Subset %>% mutate(min_conc = min(conc),
                                            max_conc = max(conc),
                                            n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct() %>%
                dplyr::mutate(ctr_median = Control$Median, ctr_nmad = Control$nMad)

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$percent_responder_pi,
                          bmed=Control$Median,
                          cutoff=5*Control$nMad, #only fit curve if it exeeds 5*nMad;
                          onesd=Control$nMad/1.349*3,  #calculate bmd at a BMR of 3 nMad
                          approach = "CV", endpoint = "percent_responder_pi")

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels=c("cnst", "hill",  "gnls"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = FALSE))

              if(is.null(newLine) | class(newLine)=="try-error"){
                newLine = try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine)=="try-error"){
                newLine=row
              }
              rownames(newLine) <- ""

              if(length(cv_tcpl$find(query = mongoQuery(pg_id = PG, chem_id = Chem, cell_type = cell, stype = sample_type, approach = "CV", endpoint = "percent_responder_pi"), fields = '{"_id" : 1}')) > 0){
                cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                cv_tcpl$remove(query = mongoQuery(pg_id = PG, chem_id = Chem, stype = sample_type, approach = "CV", endpoint = "percent_responder_pi") )
              }
              cv_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")

            }#for each stype

          }#for each Chem

        }else{

          ############## fit cell count ##############
          message(paste("******", "cell count"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype == "vehicle control" & rel_cell_count >= 50 & n_fields == nCV)

          Control <- Control %>% summarise(Median = median(rel_cell_count, na.rm = T),
                                           nMad = mad(rel_cell_count, constant = 1.4826, na.rm = T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct()

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$rel_cell_count,
                          bmed=Control$Median,
                          cutoff=2*Control$nMad, #only fit curve if it exceeds 2*nMad;
                          onesd=50/1.349, #onesd is 50% divided by 1.349 to have a BMR of 50%
                          approach = "HTPP", endpoint = "rel_cell_count") #note `approach` is "HTPP" since these data come from HTPP cell painting plates

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = TRUE))

              if(is.null(newLine) | class(newLine) == "try-error"){
                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine) == "try-error"){
                newLine <- row
              }
              rownames(newLine) <- ""

              if(length(cv_tcpl$find(query = mongoQuery(pg_id = PG, chem_id = Chem, cell_type = cell, stype = sample_type, approach = "HTPP", endpoint = "rel_cell_count"), fields = '{"_id" : 1}')) > 0){
                cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                cv_tcpl$remove(query = mongoQuery(pg_id = PG, chem_id = Chem, stype = sample_type, approach = "HTPP", endpoint = "rel_cell_count") )
              }
              cv_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")
            } #for each stype

          }#for each Chem

        }
        htpp_chem <- mongo(collection="htpp_chem", url=mongoUrl, verbose=getOption("verbose"))
        chemCount <- htpp_chem$count(query = mongoQuery(stype = c("test sample", "viability positive control", "reference chemical")))
        cvCount<-cv_tcpl$count(query = mongoQuery(cell_type = cell))
        if(cvCount != chemCount){
          warning(paste("Expected", (chemCount), "documents in cell viability collection, based on htpp_chem, for", cell, "instead there are", cvCount, "documents."))
        }
      }#for each cell type
    }#for each plate group
  }else{
    if(file.exists(paste(json_collection_path,"htpp_well.JSON",sep="/"))){
      htpp_well <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well.JSON",sep="/")))
      htpp_well[htpp_well == "NA"] <- NA
    } else {
      message("htpp_well.JSON not found, check json_collection_path parameter.")
    }

    if(rerun == TRUE){
      cv_tcpl <- toJSON(data.table())
      write(cv_tcpl, paste(json_collection_path,"cv_tcpl.JSON", sep = "/"))
    }



    #Check number of plate groups and chemicals
    #load in htpp_well_trt as metaData (skipping unneeded conversion step from use-db version)
    if(file.exists(paste(json_collection_path,"htpp_well_trt.JSON",sep="/"))){
      metaData <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_trt.JSON",sep="/")))
      metaData[metaData == "NA"] <- NA
    } else {
      message("htpp_well_trt.JSON not found, check json_collection_path parameter.")
    }

    #load in htpp_chem
    if(file.exists(paste(json_collection_path,"htpp_chem.JSON",sep="/"))){
      htpp_chem <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_chem.JSON",sep="/")))
      htpp_chem[htpp_chem == "NA"] <- NA
    } else {
      message("htpp_chem.JSON not found, check json_collection_path parameter.")
    }

    #load in tcpl collection
    cv_tcpl <-data.table(fromJSON(txt=paste(json_collection_path,"cv_tcpl.JSON",sep="/")))
    cv_tcpl[cv_tcpl == "NA"] <- NA

    if(cell_viability==TRUE){
      message("cell_viability = TRUE: will generate cv_tcpl collection from cell viability data from the cv_well collection")
      if(file.exists(paste(json_collection_path,"cv_well.JSON",sep="/"))){
        cv_well <-data.table(fromJSON(txt=paste(json_collection_path,"cv_well.JSON",sep="/")))
        cv_well[cv_well == "NA"] <- NA
      } else {
        stop("The cv_well collection is empty but is required to create cv_tcpl. Please ensure cv_well is created before proceeding.")
      }
    } else {
      message("cell_viability = FALSE: will generate cv_tcpl collection using relative cell counts from the htpp_well collection")
    }

    for(PG in as.character(unique(metaData[ ,pg_id]))){
      message(paste("*************", PG, "*************"))

      if(cell_viability == TRUE){
        CVData_all <- cv_well[pg_id == PG,]

        #Use htpp_well_trt to add in dtxsid to CVData_all
        htpp_well_trt <- htpp_well_trt[,c("sample_id", "dtxsid"), with = FALSE]
        CVData_all <- CVData_all %>% dplyr::inner_join(htpp_well_trt, by = c("sample_id", "dtxsid"))

        #Use htpp_chem to add chem_name to CVData_all
        htpp_chem <- htpp_chem[, c("dtxsid", "chem_name")]
        CVData_all <- CVData_all %>% dplyr::inner_join(htpp_chem, by = c("dtxsid"))

      }else{
        CVData_all <- htpp_well[pg_id == PG,]
      }

      for(cell in unique(CVData_all[, cell_type])){

        cat("Generating CV data for", cell, "cells\n")

        CVData <- CVData_all[cell_type == cell]

        CompoundList <- unique(CVData$chem_id[which(CVData$stype != "vehicle control")])

        if(cell_viability == TRUE){

          ############## fit cell count ##############
          message(paste("******", "cell count"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype == "vehicle control" & rel_cell_count >= 50 & n_fields == nCV)

          Control <- Control %>% summarise(Median = median(rel_cell_count, na.rm = T),
                                           nMad = mad(rel_cell_count, constant = 1.4826, na.rm = T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct()

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$rel_cell_count,
                          bmed=Control$Median,
                          cutoff=2*Control$nMad, #only fit curve if it exceeds 2*nMad;
                          onesd=50/1.349, #onesd is 50% divided by 1.349 to have a BMR of 50%
                          approach = "CV", endpoint = "rel_cell_count")

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = TRUE))

              if(is.null(newLine) | class(newLine) == "try-error"){
                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine) == "try-error"){
                newLine <- row
              }
              rownames(newLine) <- ""

              if(cv_tcpl[pg_id = PG &  chem_id == Chem & cell_type == cell & stype == sample_type & approach == "CV" & endpoint == "rel_cell_count"] %>%
                 nrow() > 0) {
                cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                cv_tcpl<-cv_tcpl[!(pg_id == PG & chem_id == Chem & stype == sample_type & approach == "CV" & endpoint == "rel_cell_count")]
              }
              cv_tcpl<-rbind(cv_tcpl, newLine)
            } #for each stype

          }#for each Chem

          ############## fit  PI response ##############
          message(paste("******", "PI response"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype=="vehicle control" & n_cells_keep >= minObjects & n_fields == nCV)

          Control <- Control%>% summarise(Median = median(percent_responder_pi, na.rm=T),
                                          nMad   = mad(percent_responder_pi, constant=1.4826, na.rm=T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem & n_cells_keep >= minObjects) %>% filter(stype == sample_type)

              Metadata <- Subset %>% mutate(min_conc = min(conc),
                                            max_conc = max(conc),
                                            n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct() %>%
                dplyr::mutate(ctr_median = Control$Median, ctr_nmad = Control$nMad)

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$percent_responder_pi,
                          bmed=Control$Median,
                          cutoff=5*Control$nMad, #only fit curve if it exeeds 5*nMad;
                          onesd=Control$nMad/1.349*3,  #calculate bmd at a BMR of 3 nMad
                          approach = "CV", endpoint = "percent_responder_pi")

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels=c("cnst", "hill",  "gnls"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = FALSE))

              if(is.null(newLine) | class(newLine)=="try-error"){
                newLine = try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine)=="try-error"){
                newLine=row
              }
              rownames(newLine) <- ""

              if(length(cv_tcpl$find(query = mongoQuery(pg_id = PG, chem_id = Chem, cell_type = cell, stype = sample_type, approach = "CV", endpoint = "percent_responder_pi"), fields = '{"_id" : 1}')) > 0){
                cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                cv_tcpl$remove(query = mongoQuery(pg_id = PG, chem_id = Chem, stype = sample_type, approach = "CV", endpoint = "percent_responder_pi") )
              }
              cv_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")

            }#for each stype

          }#for each Chem

        }else{

          ############## fit cell count ##############
          message(paste("******", "cell count"))

          nCV<-as.integer(max(CVData$n_fields, na.rm = TRUE))
          Control <- CVData %>% filter(stype == "vehicle control" & rel_cell_count >= 50 & n_fields == nCV)

          Control <- Control %>% summarise(Median = median(rel_cell_count, na.rm = T),
                                           nMad = mad(rel_cell_count, constant = 1.4826, na.rm = T))

          for(Chem in CompoundList){
            message(Chem)
            #Subset <- CVData %>% filter(chem_id == Chem)
            if(Chem == "DMSO"){
              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
            }

            for(sample_type in unique(CVData[chem_id == Chem & stype != "vehicle control", stype])){

              Subset <- CVData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(assay, pg_id, cell_type, stype, chem_id, dtxsid, chem_name, min_conc,  max_conc, n_conc) %>% distinct()

              row <- list(assay = Metadata$assay,
                          pg_id = Metadata$pg_id,
                          cell_type = Metadata$cell_type,
                          stype = Metadata$stype,
                          chem_id = Metadata$chem_id,
                          dtxsid = Metadata$dtxsid,
                          chem_name = Metadata$chem_name,
                          min_conc = Metadata$min_conc,
                          max_conc = Metadata$max_conc,
                          n_conc = Metadata$n_conc,
                          ctr_median = Control$Median,
                          ctr_nmad = Control$nMad,
                          conc=Subset$conc,
                          resp=Subset$rel_cell_count,
                          bmed=Control$Median,
                          cutoff=2*Control$nMad, #only fit curve if it exceeds 2*nMad;
                          onesd=50/1.349, #onesd is 50% divided by 1.349 to have a BMR of 50%
                          approach = "HTPP", endpoint = "rel_cell_count") #note `approach` is "HTPP" since these data come from HTPP cell painting plates

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill"),
                                          conthits = F, aicc = F, force.fit = FALSE, bidirectional = TRUE))

              if(is.null(newLine) | class(newLine) == "try-error"){
                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
              }
              if(class(newLine) == "try-error"){
                newLine <- row
              }
              rownames(newLine) <- ""

              if(rerun == FALSE){
                if((cv_tcpl[pg_id == PG & chem_id == Chem & cell_type == cell & stype == sample_type & approach == "HTPP" & endpoint == "rel_cell_count"] %>% nrow()) > 0){
                  cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
                  cv_tcpl<-cv_tcpl[!(pg_id == PG & chem_id == Chem & stype == sample_type & approach == "HTPP" & endpoint == "rel_cell_count")]
                }
                cv_tcpl<-rbind(cv_tcpl,newLine)
              }else{
                cv_tcpl<-rbind(cv_tcpl,newLine)
              }
            } #for each stype

          }#for each Chem

        }
        cv_tcplJSON<-toJSON(cv_tcpl, na = "string", digits = 8)
        write(cv_tcplJSON, file=paste(json_collection_path,"cv_tcpl.JSON",sep="/"))

        chemCount <- htpp_chem[stype %in% c("test sample", "viability positive control", "reference chemical")] %>% nrow()
        cvCount<-cv_tcpl[cell_type==cell]%>% nrow()
        if(cvCount != chemCount){
          warning(paste("Expected", (chemCount), "documents in cell viability collection, based on htpp_chem, for", cell, "instead there are", cvCount, "documents."))
        }
      }#for each cell type
    }#for each plate group
  }
}
