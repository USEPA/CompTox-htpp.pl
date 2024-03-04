#' Create mongo collection for cell viability by well (cv_well) from well treated collection
#'
#' @param file_path character string: file path to the top level directory of cell viability Harmony files for an HTPP dataset (i.e., the directory above plate-level directories)
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing cv_well collection and reinsert; FALSE by default
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
generate_cvWell <- function(file_path, mongoUrl, rerun=FALSE){
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

  message(mongo(collection="cv_well", url=mongoUrl, verbose=getOption("verbose"))$count())

}


#' Creates and populates cell viability bmc (cv_bmc) collection in mongo
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun TRUE will drop existing collection and reinsert; FALSE by default
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

generate_cvBMC <- function(mongoUrl, rerun=FALSE){

  #---------------------------------------------------------------------------------------------#
  # 1. Check what is in cv_bmc and delete if rerun == TRUE
  #---------------------------------------------------------------------------------------------#

  cv_bmc <- mongo(collection="cv_bmc", url=mongoUrl, verbose=getOption("verbose"))


  if(rerun == TRUE){
    cv_bmc$drop()
  }

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

  cv_bmc <- mongo(collection="cv_bmc", url =mongoUrl, verbose=getOption("verbose"))
  bmcCount <- cv_bmc$count()
  cv_tcpl <- mongo(collection="cv_tcpl", url =mongoUrl, verbose=getOption("verbose"))
  tcplCount <-cv_tcpl$count()

  if(bmcCount != tcplCount){
    warning(paste("Expected", tcplCount, "documents in cv_bmc, based on cv_tcpl, instead there are", bmcCount, "documents."))
  }

}


#' Creates cell viability collection cv_tcpl based on well and chem data
#'
#' @param cell_viability boolean: if cell_viability = TRUE, the CVData_all object in line 218 should pull out data from the cv_well collection, otherwise it will pull data from the htpp_well collection
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing cv_tcpl collection and reinsert; FALSE by default
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
#'
#' @export generate_cvTcpl
#'
generate_cvTcpl<-function(cell_viability, mongoUrl, rerun=FALSE){
  #setup
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

    CVData_all <- data.table(mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))

    for(cell in unique(CVData_all[, cell_type])){

      cat("Generating CV data for", cell, "cells\n")

      CVData <- CVData_all[cell_type == cell]

      CompoundList <- unique(CVData$chem_id[which(CVData$stype != "vehicle control")])

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

          if(cell_viability == TRUE){
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
                        cutoff=5*Control$nMad, #only fit curve if it exceeds 5*nMad;
                        onesd=(Control$nMad/1.349)*3, #to calculate the BMC for a BMR of 3 nMad
                        approach = "CV", endpoint = "rel_cell_count")

            newLine <- NULL
            newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill", "gnls"),
                                        conthits = F, aicc = F, force.fit = FALSE, bidirectional = FALSE))

            if(is.null(newLine) | class(newLine) == "try-error"){
              newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = FALSE, fitmodels=c("cnst")))
            }
            if(class(newLine) == "try-error"){
              newLine <- row
            }
            rownames(newLine) <- ""
          } else{
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
          }

          if(length(cv_tcpl$find(query = mongoQuery(pg_id = PG, chem_id = Chem, cell_type = cell, stype = sample_type, approach = "CV", endpoint = "rel_cell_count"), fields = '{"_id" : 1}')) > 0){
            cat(PG, " ", Chem, " ", cell, " ", "old results deleted!\n")
            cv_tcpl$remove(query = mongoQuery(pg_id = PG, chem_id = Chem, stype = sample_type, approach = "CV", endpoint = "rel_cell_count") )
          }
          cv_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")
        } #for each stype

      }#for each Chem
      htpp_chem <- mongo(collection="htpp_chem", url=mongoUrl, verbose=getOption("verbose"))
      chemCount <- htpp_chem$count(query = mongoQuery(stype = c("test sample", "viability positive control", "reference chemical")))
      cvCount<-cv_tcpl$count(query = mongoQuery(cell_type = cell))
      if(cvCount != chemCount){
        warning(paste("Expected", (chemCount), "documents in cell viability collection, based on htpp_chem, for", cell, "instead there are", cvCount, "documents."))
      }
    }#for each cell type
  }#for each plate group
}
