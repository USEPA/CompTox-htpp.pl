#' Global Mahalanobis curve fitting, adds global Mahalanobis fit data to htpp_tcpl collection
#'
#' Retrieve data from the htpp_global_mah MongoDB collection for one plate group at a time
#' Identify baseline levels using the lowest two concentrations of test chemicals
#' Only include wells with relative cell counts > 50.
#' Calculate Mean and nMad.
#' For each chemical, subset the data based on cell viability flags as well as whether 'n_cells_keep' > 'minObjects' parameter
#' Run run the tcplfit2::concRespCore on the subset
#'
#' Typical parameters are:
#' cutoff = 1 * nMad of controls
#' onesd = nMad / 1.349  (to calculate the BMC for a BMR of 1 nMad)
#' fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5")
#' bidirectional = FALSE
#'
#' Insert results into the htpp_tcpl MongoDB collection for 'approach' == "global"
#'
#' @param minObjects numeric: The minimum number of objects used to filter the dataset for analysis
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun  boolean: rerun = TRUE will drop existing htpp_tcpl collection for global mah values (htpp_tcpl$remove(query=mongoQuery(approach="global", endpoint="global")) and reinsert; FALSE by default
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
#' @import tcplfit2
#'
#'
#' @export curveFit_htppGlobalMah
#'
curveFit_htppGlobalMah <- function(minObjects, mongoUrl="", rerun=FALSE, use_db=TRUE, json_collection_path=""){
 if(use_db==TRUE){
  #---------------------------------------------------------------------------------------------#
  # 1. Setup
  #---------------------------------------------------------------------------------------------#
  htpp_global_mah <- mongo(collection = "htpp_global_mah", url = mongoUrl, verbose = getOption("verbose"))

  if(htpp_global_mah$count() <1 ){
    stop("The htpp_global_mah collection is empty but is required to fill htpp_tcpl. Please ensure htpp_global_mah is created before proceeding.")
  }

  #Check collection and delete any global mahalanobis distance curve fits to htpp_tcpl if rerun == TRUE
  htpp_tcpl <- mongo(collection = "htpp_tcpl", url = mongoUrl, verbose = getOption("verbose"))



  if(rerun == TRUE){
    htpp_tcpl$remove(query=mongoQuery(approach="global", endpoint="global") )
  }

  #---------------------------------------------------------------------------------------------#
  # 2. tcplfit for Global Mahalanobis distances (~ 12 min) - FOR EACH CELL TYPE
  #---------------------------------------------------------------------------------------------#

  #select plate groups
  PlateGroupList <- as.character(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "pg_id"))


  #run tcpl curve fitting for each cell type
  for(PG in PlateGroupList){
    print(paste("************", PG, "************"))

    ## get distance information
    DistanceData_all <- data.table(mongo(collection="htpp_global_mah", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(pg_id = PG)))
    DistanceData_all <- DistanceData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

    ## get CV data
    CV_BMC_all <- data.table(mongo(collection = "cv_bmc", url = mongoUrl, verbose = getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
    CV_BMC_all <- CV_BMC_all %>% select(pg_id, stype, cell_type, chem_id, cv_noec_dose_level, cv_flag)

    for(cell in unique(DistanceData_all[, cell_type])){

      #print message
      cat("Beginning tcpl curve fits of global mahalanobis distance data for", cell, "cells and", PG, "plate group\n")

      #filter by cell type
      DistanceData <- as.data.frame(DistanceData_all[cell_type == cell, ])
      CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])

      # attach
      DistanceData <- DistanceData %>% left_join(CV_BMC) %>%
        #modify the information for null chemicals, as they can not be properly processed otherwise
        dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                      cv_flag = ifelse(stype == "null", F, cv_flag))

      #define baseline
      Control <- DistanceData %>% filter((stype == "test sample") & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep > minObjects)
      Control <- Control %>% dplyr::summarise(Mean = mean(d, na.rm = T),
                                              SD = sd(d, na.rm = T))

      CompoundList <- unique(DistanceData$chem_id)

      for(Chem in CompoundList){

        print(Chem)

        #add additional filter to handle DMSO as test chemical
        if(Chem == "DMSO"){
          Subset <- DistanceData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
        }

        #add logic to handle test chemicals that are also ref chems
        for(sample_type in unique(DistanceData$stype[DistanceData$chem_id == Chem & DistanceData$stype != "vehicle control"])){

          Subset <- DistanceData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

          if(is.na(Subset$cv_flag[1])){ #there were too little concentrations to run CV curve fit
            Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
              dplyr::summarise(min_conc = min(conc),
                               max_conc = NA,
                               n_conc = length(unique(conc)))

            Subset <- Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

          }else if(Subset$cv_flag[1]){

            Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
              dplyr::summarize(min_conc = min(conc),
                               max_conc = NA,
                               n_conc = 0) %>% distinct()

            Subset = Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

          }else{

            Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

            Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                 max_conc = max(conc),
                                                 n_conc = length(unique(conc))) %>%
              select(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc) %>% distinct()
          }

          row <- list(pg_id = Metadata$pg_id,
                      stype = Metadata$stype,
                      cell_type = Metadata$cell_type,
                      chem_id = Metadata$chem_id,
                      min_conc = Metadata$min_conc,
                      max_conc = Metadata$max_conc,
                      n_conc = Metadata$n_conc,
                      ctr_mean = Control$Mean,
                      ctr_sd = Control$SD,
                      conc = Subset$conc,
                      resp = Subset$d,
                      bmed = Control$Mean,
                      cutoff = Control$SD, #1 nMad was used for the HitDef manuscript
                      onesd = Control$SD/1.349, #we want a BMR of 1*nMad; (1 nMad was also used for the HitDef manuscript)
                      approach = "global",
                      endpoint = "global")

          newLine <- NULL
          newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
                                      conthits = T, aicc = F, force.fit = FALSE, bidirectional = FALSE))

          if(is.null(newLine) | class(newLine)=="try-error"){

            newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
          }else if(class(newLine) == "try-error"){

            newLine <- row
          }else{

            rownames(newLine) <- ""
          }

          if(length(htpp_tcpl$find(query = mongoQuery(pg_id = PG, cell_type = cell, chem_id = Chem, stype = sample_type, approach = "global", endpoint = "global"), fields = '{"_id" : 1}')) > 0){

            cat(PG, " ", cell, " ", Chem, " old results deleted!\n")

            htpp_tcpl$remove(query = mongoQuery(pg_id = PG, cell_type = cell, chem_id = Chem, stype = sample_type, approach = "global", endpoint = "global"))
          }

          htpp_tcpl$insert(newLine, auto_unbox=TRUE, na = "null")

        } #for each stype

      } #for each Chem

      cat("A total of", htpp_tcpl$count(query=mongoQuery(pg_id = PG, cell_type = cell, approach = "global", endpoint = "global")), "documents of global mahalanobis distance tcpl fits added to htpp_tcpl for", cell, "cells and", PG, "plate group\n")

    } #for each cell type

  } #for each plate group

  #add index
  htpp_tcpl$index(add = '{"pg_id" : 1, "stype": 1, "cell_type": 1, "chem_id": 1, "approach": 1, "endpoint": 1}')

  #---------------------------------------------------------------------------------------------#
  # 3. Check if collection is complete
  #---------------------------------------------------------------------------------------------#

  gMahCount <- list()
  i <- 1
  for(PG in PlateGroupList){
    gMahCount[i] <- length(htpp_global_mah$distinct(key = "chem_id", query = mongoQuery(stype = c("test sample", "null", "reference chemical", "viability positive control"), pg_id = PG)))
    i <- i+1
  }
  gMahCount <- do.call(sum, gMahCount)
  tcplCount <- htpp_tcpl$count(query = mongoQuery(approach = "global", endpoint = "global"))

  if( gMahCount != tcplCount){
    warning(paste("Expected", gMahCount, "documents with global approach and endpoint in htpp_tcpl, based on htpp_global_mah.  Instead, there are", tcplCount, "documents."))
  }

  for(cell in unique(DistanceData_all[, cell_type])){
    cat(htpp_tcpl$count(query = mongoQuery(cell_type = cell, approach = "global", endpoint = "global")), "global mahalanobis tcpl fit documents added to htpp_tcpl for", cell, "cells\n")
  }
 } else {
   #---------------------------------------------------------------------------------------------#
   # 1. Setup
   #---------------------------------------------------------------------------------------------#
   if(file.exists(paste(json_collection_path,"htpp_global_mah.JSON",sep="/"))){
     htpp_global_mah <- data.table(fromJSON(txt=paste(json_collection_path,"htpp_global_mah.JSON",sep="/")))
     htpp_global_mah[htpp_global_mah == "NA"] <- NA
   } else {
     stop("htpp_global_mah.JSON not found, check json_collection_path parameter.")
   }

   #Check collection and delete any global mahalanobis distance curve fits to htpp_tcpl if rerun == TRUE

   if(file.exists(paste(json_collection_path,"htpp_tcpl.JSON",sep="/")) & rerun == TRUE){
     htpp_tcpl <- data.table(fromJSON(txt=paste(json_collection_path,"htpp_tcpl.JSON",sep="/")))
     htpp_tcpl[htpp_tcpl == "NA"] <- NA

     htpp_tcpl <- htpp_tcpl[approach != "global" & endpoint != "global",]
   } else if(file.exists(paste(json_collection_path,"htpp_tcpl.JSON",sep="/")) == FALSE){
     htpp_tcpl <- data.table()
   }else if(rerun == FALSE){
     htpp_tcpl <- data.table(fromJSON(txt=paste(json_collection_path,"htpp_tcpl.JSON",sep="/")))
     htpp_tcpl[htpp_tcpl == "NA"] <- NA
   }

   if(file.exists(paste(json_collection_path,"cv_bmc.JSON",sep="/"))){
     cv_bmc <- data.table(fromJSON(txt=paste(json_collection_path,"cv_bmc.JSON",sep="/")))
     cv_bmc[cv_bmc == "NA"] <- NA
   } else {
     message("cv_bmc.JSON not found, check json_collection_path parameter.")
   }


   #---------------------------------------------------------------------------------------------#
   # 2. tcplfit for Global Mahalanobis distances (~ 12 min) - FOR EACH CELL TYPE
   #---------------------------------------------------------------------------------------------#

   #select plate groups
   PlateGroupList <- htpp_global_mah[,pg_id]%>%unique()


   #run tcpl curve fitting for each cell type
   for(PG in PlateGroupList){
     print(paste0("************", PG, "************"))

     ## get distance information
     DistanceData_all <- htpp_global_mah[pg_id == PG,]
     DistanceData_all <- DistanceData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

     ## get CV data
     CV_BMC_all <- cv_bmc[pg_id == PG,]
     CV_BMC_all <- CV_BMC_all %>% select(pg_id, stype, cell_type, chem_id, cv_noec_dose_level, cv_flag)

     for(cell in unique(DistanceData_all[, cell_type])){

       #print message
       cat("Beginning tcpl curve fits of global mahalanobis distance data for", cell, "cells and", PG, "plate group\n")

       #filter by cell type
       DistanceData <- as.data.frame(DistanceData_all[cell_type == cell, ])
       CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])

       # attach
       DistanceData <- DistanceData %>% left_join(CV_BMC) %>%
         #modify the information for null chemicals, as they can not be properly processed otherwise
         dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                       cv_flag = ifelse(stype == "null", F, cv_flag))

       #define baseline
       Control <- DistanceData %>% filter((stype == "test sample") & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep > minObjects)
       Control <- Control %>% dplyr::summarise(Mean = mean(d, na.rm = T),
                                               SD = sd(d, na.rm = T))

       CompoundList <- unique(DistanceData$chem_id)

       for(Chem in CompoundList){

         print(Chem)

         #add additional filter to handle DMSO as test chemical
         if(Chem == "DMSO"){
           Subset <- DistanceData %>% filter(chem_id == Chem) %>% filter(stype == "test sample")
         }

         #add logic to handle test chemicals that are also ref chems
         for(sample_type in unique(DistanceData$stype[DistanceData$chem_id == Chem & DistanceData$stype != "vehicle control"])){

           Subset <- DistanceData %>% filter(chem_id == Chem) %>% filter(stype == sample_type)

           if(is.na(Subset$cv_flag[1])){ #there were too little concentrations to run CV curve fit
             Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
               dplyr::summarise(min_conc = min(conc),
                                max_conc = NA,
                                n_conc = length(unique(conc)))

             Subset <- Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

           }else if(Subset$cv_flag[1]){

             Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
               dplyr::summarize(min_conc = min(conc),
                                max_conc = NA,
                                n_conc = 0) %>% distinct()

             Subset = Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

           }else{

             Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

             Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                  max_conc = max(conc),
                                                  n_conc = length(unique(conc))) %>%
               select(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc) %>% distinct()
           }

           row <- list(pg_id = Metadata$pg_id,
                       stype = Metadata$stype,
                       cell_type = Metadata$cell_type,
                       chem_id = Metadata$chem_id,
                       min_conc = Metadata$min_conc,
                       max_conc = Metadata$max_conc,
                       n_conc = Metadata$n_conc,
                       ctr_mean = Control$Mean,
                       ctr_sd = Control$SD,
                       conc = Subset$conc,
                       resp = Subset$d,
                       bmed = Control$Mean,
                       cutoff = Control$SD, #1 nMad was used for the HitDef manuscript
                       onesd = Control$SD/1.349, #we want a BMR of 1*nMad; (1 nMad was also used for the HitDef manuscript)
                       approach = "global",
                       endpoint = "global")

           newLine <- NULL
           newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
                                       conthits = T, aicc = F, force.fit = FALSE, bidirectional = FALSE))

           if(is.null(newLine) | class(newLine)=="try-error"){

             newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst")))
           }else if(class(newLine) == "try-error"){

             newLine <- row
           }else{

             rownames(newLine) <- ""
           }

           if(nrow(htpp_tcpl) < 1){
             htpp_tcpl <- rbind(htpp_tcpl, newLine)
           }else if(htpp_tcpl[pg_id == PG & cell_type == cell & chem_id == Chem & stype == sample_type & approach == "global" & endpoint == "global"] %>% nrow() > 0){
             cat(PG, " ", cell, " ", Chem, " old results deleted!\n")
             htpp_tcpl <- htpp_tcpl[!c(pg_id == PG & cell_type == cell & chem_id == Chem & stype == sample_type & approach == "global" & endpoint == "global"),] #remove existing document
             htpp_tcpl <- rbind(htpp_tcpl, newLine) #add new document
           }else{
             htpp_tcpl <- rbind(htpp_tcpl, newLine)
           }

         } #for each stype

       } #for each Chem

       cat("A total of", htpp_tcpl[pg_id == PG & cell_type == cell & approach == "global" & endpoint == "global"] %>% nrow(), "documents of global mahalanobis distance tcpl fits added to htpp_tcpl for", cell, "cells and", PG, "plate group\n")

     } #for each cell type

   } #for each plate group

   tcplJSON<-toJSON(htpp_tcpl, na="string", digits = 8)
   write(tcplJSON, file=paste(json_collection_path,"htpp_tcpl.JSON",sep="/"))

   #---------------------------------------------------------------------------------------------#
   # 3. Check if collection is complete
   #---------------------------------------------------------------------------------------------#

   gMahCount <- list()
   i <- 1
   for(PG in PlateGroupList){
     gMahCount[i] <- htpp_global_mah[stype %in% c("test sample", "null", "reference chemical", "viability positive control") & pg_id==PG, chem_id] %>% unique() %>% length()
     i <- i+1
   }
   gMahCount <- do.call(sum, gMahCount)
   tcplCount <- htpp_tcpl[approach=="global" & endpoint=="global"] %>% nrow()

   if( gMahCount != tcplCount){
     warning(paste("Expected", gMahCount, "documents with global approach and endpoint in htpp_tcpl, based on htpp_global_mah.  Instead, there are", tcplCount, "documents."))
   }

   for(cell in unique(DistanceData_all[, cell_type])){
     cat(htpp_tcpl[cell_type == cell & approach == "global" & endpoint == "global"] %>% nrow(), "global mahalanobis tcpl fit documents added to htpp_tcpl for", cell, "cells\n")

   }
  }
}



#' Category Mahalanobis curve fitting, adds category Mahalanobis fit data to htpp_tcpl collection
#'
#' Retrieve data from the htpp_global_mah MongoDB collection for one plate group at a time
#' FOR EACH category:**
#' Subset data for category of interest
#' Identify baseline levels using the lowest two concentrations of test chemicals
#' Only include wells with relative cell counts > 50.
#' Calculate Mean and nMad.
#' For each chemical, subset the data based on cell viability flags as well as whether 'n_cells_keep' > 'minObjects' parameter
#' Run run the tcplfit2::concRespCore on the subset
#'
#' Typical parameters are:
#' cutoff = 1 * nMad of controls
#' onesd = nMad / 1.349  (to calculate the BMC for a BMR of 1 nMad)
#' fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5")
#' bidirectional = FALSE
#'
#' Insert results into the htpp_tcpl MongoDB collection for 'approach' == "category"
#'
#'
#'
#' @param minObjects numeric: The minimum number of objects used to filter the dataset for analysis
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing htpp_tcpl collection for global mah values (htpp_tcpl$remove(query=mongoQuery(approach="category")) and reinsert; FALSE by default
#' @param nThreads numeric: the number of threads to use for processing; default is 1
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
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import tcplfit2
#'
#'
#' @export curveFit_htppCatMah
#'
curveFit_htppCatMah <- function(minObjects, mongoUrl="", rerun=FALSE, nThreads=1, use_db=T, json_collection_path=""){

  if(use_db==T){
  htpp_cat_mah<-mongo(collection = "htpp_cat_mah", url = mongoUrl, verbose = getOption("verbose"))
  if(htpp_cat_mah$count()<1){
    stop("The htpp_cat_mah collection is empty but is required to fill in category information for htpp_tcpl. Please ensure htpp_cat_mah is created before proceeding.")
  }

  #Check collection and delete any category-level mahalanobis distance curve fits to htpp_tcpl if rerun == TRUE
  htpp_tcpl <- mongo(collection = "htpp_tcpl", url = mongoUrl, verbose = getOption("verbose"))

  Table1_all = data.table(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$find(fields='{"_id":0}'))

  if(rerun == TRUE){
    htpp_tcpl$remove(query=mongoQuery(approach="category"))
  }

  #---------------------------------------------------------------------------------------------#
  # 2. tcplfit for Category-level Mahalanobis distances (ca ~ 8 min/pg) - FOR EACH CELL TYPE
  #---------------------------------------------------------------------------------------------#


  #To-Do: Explore ways to parallelize this part of the code.

  #select plate groups
  PlateGroupList = as.character(unique(Table1_all[, pg_id])) #35 plate groups in this dataset

  my.cluster <- parallel::makeCluster(
    nThreads,
    type = "PSOCK"
  )

  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  #run tcpl curve fitting for each cell type
  for(PG in PlateGroupList){

    print(paste("************", PG, "************"))

    ## get distance information
    DistanceData_all <- data.table(mongo(collection = "htpp_cat_mah", url = mongoUrl, verbose = getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
    DistanceData_all <- DistanceData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

    ## get CV data
    CV_BMC_all <- data.table(mongo(collection = "cv_bmc", url = mongoUrl, verbose = getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
    CV_BMC_all <- CV_BMC_all %>% select(pg_id, stype, cell_type, chem_id, cv_noec_dose_level, cv_flag)

    for(cell in unique(DistanceData_all[, cell_type])){

      #keep time
      tic()

      #print message
      cat("Beginning tcpl curve fits of category-level mahalanobis distance data for", cell, "cells\n")

      #filter by cell type
      DistanceData <- as.data.frame(DistanceData_all[cell_type == cell, ])
      CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])

      # attach
      DistanceData <- DistanceData %>% left_join(CV_BMC) %>%
        #modify the information for null chemicals, as they can not be properly processed otherwise
        dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                      cv_flag = ifelse(stype == "null", F, cv_flag))

      #for each category
      for(Category in unique(DistanceData$category_name_r)){

        cat("Curve fitting for", Category, "\n")

        DistanceData_1Category <- DistanceData %>% filter(category_name_r == Category)


        #define baseline
        Control <- DistanceData_1Category %>% filter(stype == "test sample" & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep>minObjects)
        Control <- Control%>% dplyr::summarise(Mean = mean(d, na.rm = T), SD = sd(d, na.rm = T))

        #CompoundList = unique(DistanceData_1Category$chem_id)

        #for each chem
        #parallelize here?
        Chems<-unique(DistanceData_1Category$chem_id)

        concResp <- foreach(x=Chems, .combine = 'rbind', .packages = c('tidyr', 'dplyr', 'tcplfit2')) %dopar% {

          #add additional filter to handle DMSO as test chemical
          if(x == "DMSO"){
            Subset <- DistanceData_1Category %>% filter(chem_id == x) %>% filter(stype == "test sample")
          }

          #add logic to handle test chemicals that are also ref chems
          for(sample_type in unique(DistanceData_1Category$stype[DistanceData_1Category$chem_id == x & DistanceData_1Category$stype != "vehicle control"])){

            Subset <- DistanceData_1Category %>% filter(chem_id == x) %>% filter(stype == sample_type)

            if(is.na(Subset$cv_flag[1])){ #if there were too few concentrations to run CV curve fit
              Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                dplyr::summarise(min_conc = min(conc),
                                 max_conc = NA,
                                 n_conc = length(unique(conc)))
              Subset <- Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

            }else if(Subset$cv_flag[1]){
              Metadata <-  Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                dplyr::summarise(min_conc = min(conc),
                                 max_conc = NA,
                                 n_conc = 0) %>% distinct()
              Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

            }else{
              Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)
              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc) %>% distinct()
            }

            row<-list(pg_id = Metadata$pg_id,
                      stype = Metadata$stype,
                      cell_type = Metadata$cell_type,
                      chem_id = Metadata$chem_id,
                      min_conc = Metadata$min_conc,
                      max_conc = Metadata$max_conc,
                      n_conc = Metadata$n_conc,
                      ctr_mean = Control$Mean,
                      ctr_sd = Control$SD,
                      conc = Subset$conc,
                      resp = Subset$d,
                      bmed = Control$Mean,
                      cutoff = Control$SD, #1 nMad was used for the HitDef manuscript
                      onesd = Control$SD/1.349, #we want a BMR of 1*nMad; (1 nMad was also used for the HitDef manuscript)
                      approach = "category",
                      endpoint = Category)
            newLine <- NULL
            newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
                                        conthits = T, aicc = F, force.fit = FALSE, bidirectional = FALSE))

            if(is.null(newLine) | class(newLine) == "try-error"){

              newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels = c("cnst")))
            } else if(class(newLine) == "try-error"){

              newLine <- row
            }else{

              rownames(newLine) <- ""
            }

          } #for each stype
          newLine
        } #for each Chem


        concResp<-as_tibble(concResp)



        htpp_tcpl$insert(concResp, auto_unbox=TRUE, na = "null")


        cat("A total of", htpp_tcpl$count(query=mongoQuery(pg_id = PG, cell_type = cell, approach = "category", endpoint = Category)), "documents of category-level mahalanobis distance tcpl fits were added to htpp_tcpl for", cell, "cells and", Category, "category\n")

      } #for each category



      cat("A total of", htpp_tcpl$count(query=mongoQuery(pg_id = PG, cell_type = cell, approach = "category")), "documents of category-level mahalanobis distance tcpl fits were added to htpp_tcpl for", cell, "cells\n")

      toc()

    } #for each cell type

  } #for each plate group
  stopCluster(my.cluster)

  #---------------------------------------------------------------------------------------------#
  # 3. Check if collection is complete
  #---------------------------------------------------------------------------------------------#

  catMahCount <- list()
  i <- 1
  for(PG in PlateGroupList){
    catMahCount[i] <- length(htpp_cat_mah$distinct(key = "chem_id", query = mongoQuery(stype = c("test sample", "null", "reference chemical", "viability positive control"), pg_id = PG)))
    i <- i+1
  }
  catMahCount <- do.call(sum, catMahCount)
  tcplCount <- htpp_tcpl$count(query = mongoQuery(approach = "category"))/49

  if( catMahCount != tcplCount){
    warning(paste("Expected", catMahCount, "documents with category approach in htpp_tcpl, based on htpp_cat_mah.  Instead, there are", tcplCount, "documents."))
  }

  for(cell in unique(DistanceData_all[, cell_type])){
    cat(htpp_tcpl$count(query = mongoQuery(cell_type = cell, approach = "category")), "category-level mahalanobis tcpl fit documents added to htpp_tcpl for", cell, "cells\n")
  }

  } else {

    if(file.exists(paste(json_collection_path,"htpp_cat_mah.JSON",sep="/"))){
      htpp_cat_mah <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_cat_mah.JSON",sep="/")))
      htpp_cat_mah[htpp_cat_mah == "NA"] <- NA
    } else {
      stop("The htpp_cat_mah.JSON does not exist. Please check the json_collection_path and ensure htpp_cat_mah is created before proceeding.")
    }

    #Check collection and delete any category-level mahalanobis distance curve fits to htpp_tcpl if rerun == TRUE
    if(file.exists(paste(json_collection_path,"htpp_tcpl.JSON",sep="/")) & rerun == TRUE){
      htpp_tcpl <- data.table(fromJSON(txt=paste(json_collection_path,"htpp_tcpl.JSON",sep="/")))
      htpp_tcpl[htpp_tcpl == "NA"] <- NA
      htpp_tcpl <- htpp_tcpl[approach != "category",] #if rerun == TRUE remove category fits from collection
    }else if(file.exists(paste(json_collection_path,"htpp_tcpl.JSON",sep="/")) & rerun == FALSE){
      htpp_tcpl <- data.table(fromJSON(txt=paste(json_collection_path,"htpp_tcpl.JSON",sep="/")))
      htpp_tcpl[htpp_tcpl == "NA"] <- NA
    }else{
      htpp_tcpl <- data.table()
    }

    if(file.exists(paste(json_collection_path,"htpp_well_norm.JSON",sep="/"))){
      Table1_all <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_norm.JSON",sep="/")))
      Table1_all[Table1_all == "NA"] <- NA
    } else {
      stop("The htpp_well_norm.JSON file does not exist. Please check the json_collection_path and ensure htpp_well_norm is created before proceeding.")
    }

    if(file.exists(paste(json_collection_path,"cv_bmc.JSON",sep="/"))){
      cv_bmc  <-data.table(fromJSON(txt=paste(json_collection_path,"cv_bmc.JSON",sep="/")))
      cv_bmc[cv_bmc == "NA"] <- NA
    } else {
      stop("The cv_bmc.JSON file does not exist. Please check the json_collection_path and ensure cv_bmc is created before proceeding.")

    }


    #---------------------------------------------------------------------------------------------#
    # 2. tcplfit for Category-level Mahalanobis distances (ca ~ 8 min/pg) - FOR EACH CELL TYPE
    #---------------------------------------------------------------------------------------------#

    #select plate groups
    PlateGroupList = as.character(unique(Table1_all[, pg_id])) #35 plate groups in this dataset

    my.cluster <- parallel::makeCluster(
      nThreads,
      type = "PSOCK"
    )

    #register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)

    #run tcpl curve fitting for each cell type
    for(PG in PlateGroupList){

      print(paste("************", PG, "************"))

      ## get distance information
      DistanceData_all <- htpp_cat_mah[pg_id==PG]

      DistanceData_all <- DistanceData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

      ## get CV data
      CV_BMC_all <- cv_bmc[pg_id==PG]
      CV_BMC_all <- CV_BMC_all %>% select(pg_id, stype, cell_type, chem_id, cv_noec_dose_level, cv_flag)

      for(cell in unique(DistanceData_all[, cell_type])){

        #keep time
        tic()

        #print message
        cat("Beginning tcpl curve fits of category-level mahalanobis distance data for", cell, "cells\n")

        #filter by cell type
        DistanceData <- as.data.frame(DistanceData_all[cell_type == cell, ])
        CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])

        # attach
        DistanceData <- DistanceData %>% left_join(CV_BMC) %>%
          #modify the information for null chemicals, as they can not be properly processed otherwise
          dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                        cv_flag = ifelse(stype == "null", F, cv_flag))

        #for each category
        for(Category in unique(DistanceData$category_name_r)){

          cat("Curve fitting for", Category, "\n")

          DistanceData_1Category <- DistanceData %>% filter(category_name_r == Category)


          #define baseline
          Control <- DistanceData_1Category %>% filter(stype == "test sample" & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep>minObjects)
          Control <- Control%>% dplyr::summarise(Mean = mean(d, na.rm = T), SD = sd(d, na.rm = T))

          #CompoundList = unique(DistanceData_1Category$chem_id)

          #for each chem
          #parallelize here?
          Chems<-unique(DistanceData_1Category$chem_id)

          concResp <- foreach(x=Chems, .combine = 'rbind', .packages = c('tidyr', 'dplyr', 'tcplfit2')) %dopar% {

            #add additional filter to handle DMSO as test chemical
            if(x == "DMSO"){
              Subset <- DistanceData_1Category %>% filter(chem_id == x) %>% filter(stype == "test sample")
            }

            #add logic to handle test chemicals that are also ref chems
            for(sample_type in unique(DistanceData_1Category$stype[DistanceData_1Category$chem_id == x & DistanceData_1Category$stype != "vehicle control"])){

              Subset <- DistanceData_1Category %>% filter(chem_id == x) %>% filter(stype == sample_type)

              if(is.na(Subset$cv_flag[1])){ #if there were too few concentrations to run CV curve fit
                Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                  dplyr::summarise(min_conc = min(conc),
                                   max_conc = NA,
                                   n_conc = length(unique(conc)))
                Subset <- Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

              }else if(Subset$cv_flag[1]){
                Metadata <-  Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                  dplyr::summarise(min_conc = min(conc),
                                   max_conc = NA,
                                   n_conc = 0) %>% distinct()
                Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

              }else{
                Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)
                Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                     max_conc = max(conc),
                                                     n_conc = length(unique(conc))) %>%
                  select(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc) %>% distinct()
              }

              row<-list(pg_id = Metadata$pg_id,
                        stype = Metadata$stype,
                        cell_type = Metadata$cell_type,
                        chem_id = Metadata$chem_id,
                        min_conc = Metadata$min_conc,
                        max_conc = Metadata$max_conc,
                        n_conc = Metadata$n_conc,
                        ctr_mean = Control$Mean,
                        ctr_sd = Control$SD,
                        conc = Subset$conc,
                        resp = Subset$d,
                        bmed = Control$Mean,
                        cutoff = Control$SD, #1 nMad was used for the HitDef manuscript
                        onesd = Control$SD/1.349, #we want a BMR of 1*nMad; (1 nMad was also used for the HitDef manuscript)
                        approach = "category",
                        endpoint = Category)
              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
                                          conthits = T, aicc = F, force.fit = FALSE, bidirectional = FALSE))

              if(is.null(newLine) | class(newLine) == "try-error"){

                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels = c("cnst")))
              } else if(class(newLine) == "try-error"){

                newLine <- row
              }else{

                rownames(newLine) <- ""
              }
            } #for each stype
            newLine
          } #for each Chem


          concResp<-as_tibble(concResp)

          htpp_tcpl<-rbind(htpp_tcpl, concResp)


          cat("A total of", nrow(htpp_tcpl[pg_id == PG & cell_type == cell & approach == "category" & endpoint == Category]), "documents of category-level mahalanobis distance tcpl fits were added to htpp_tcpl for", cell, "cells and", Category, "category\n")

        } #for each category



        cat("A total of", nrow(htpp_tcpl[pg_id == PG & cell_type == cell & approach == "category"]), "documents of category-level mahalanobis distance tcpl fits were added to htpp_tcpl for", cell, "cells\n")

        toc()

      } #for each cell type

    } #for each plate group
    stopCluster(my.cluster)

    #---------------------------------------------------------------------------------------------#
    # 3. Check if collection is complete
    #---------------------------------------------------------------------------------------------#

    catMahCount <- list()
    i <- 1
    for(PG in PlateGroupList){
      catMahCount[i] <-htpp_cat_mah[stype %in% c("test sample", "null", "reference chemical", "viability positive control") & pg_id == PG] %>% pull(chem_id) %>% unique() %>% length() #find the records with a relevant stype in the current plate group, pull out the chem_id column, and see how many unique chemical ids there are.
      i <- i+1
    }
    catMahCount <- do.call(sum, catMahCount)
    tcplCount <- nrow(htpp_tcpl[approach == "category"])/49

    if( catMahCount != tcplCount){
      warning(paste("Expected", catMahCount, "documents with category approach in htpp_tcpl, based on htpp_cat_mah.  Instead, there are", tcplCount, "documents."))
    }

    for(cell in unique(DistanceData_all[, cell_type])){
      cat(nrow(htpp_tcpl[cell_type == cell & approach == "category", ]), "category-level mahalanobis tcpl fit documents added to htpp_tcpl for", cell, "cells\n")
    }
    tcplJSON<-toJSON(htpp_tcpl, na="string", digits = 8)
    write(tcplJSON, file=paste(json_collection_path,"htpp_tcpl.JSON",sep="/"))

}
}

#' Feature-level curve fitting, adds feature-level fit data to htpp_tcpl collection
#'
#' Performs concentration-response modeling on the HTPP feature-level data.
#'
#' Retrieve the Level 5 data from the htpp_well_norm MongoDB collection for one plate group at a time
#'
#' FOR EACH feature:**
#' Subset data for feature of interest
#' Identify baseline levels using the lowest two concentrations of test chemicals
#' Only include wells with relative cell counts > 50.
#' Calculate Mean and nMad.
#' For each chemical, subset the data based on cell viability flags as well as whether 'n_cells_keep' > 'minObjects' parameter
#' Run run the tcplfit2::concRespCore on the subset.
#'
#' Typical parameters are:
#' cutoff = 1 * nMad of controls
#' onesd = nMad / 1.349  (to calculate the BMC for a BMR of 1 nMad)
#' fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5")
#' bidirectional = TRUE
#'
#' Insert results into the htpp_tcpl MongoDB collection for 'approach' == "feature"
#'
#' @param minObjects numeric: The minimum number of objects used to filter the dataset for analysis
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing htpp_tcpl collection for feature (htpp_tcpl$remove(query=mongoQuery(approach="feature")) and reinsert; FALSE by default
#' @param nThreads - numeric: the number of threads to use for processing; default is 1
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
#' @import tcplfit2
#' @import foreach
#' @import parallel
#' @import doParallel
#'
#'
#' @export curveFit_htppFeature
#'
curveFit_htppFeature <- function(minObjects, mongoUrl="", rerun=FALSE, nThreads=1, use_db=T, json_collection_path=""){
  if(use_db==T){
    #---------------------------------------------------------------------------------------------#
  # 1. Connect and check collection and delete if rerun == TRUE
  #---------------------------------------------------------------------------------------------#
  htpp_well_norm <- mongo(collection = "htpp_well_norm", url=mongoUrl, verbose = getOption("verbose"))
  if(htpp_well_norm$count()<1){
    stop("The htpp_well_norm collection is empty but is required to fill in feature data for http_tcpl. Please ensure htpp_well_norm is created before proceeding.")
  }


  htpp_tcpl <- mongo(collection = "htpp_tcpl", url=mongoUrl, verbose = getOption("verbose"))

  if(rerun == TRUE){
    htpp_tcpl$remove(query=mongoQuery(approach="feature"))
  }

  #---------------------------------------------------------------------------------------------#
  # 2. tcplfit for Feature-level data (ca 4 h / plategroup) - FOR EACH CELL TYPE
  #---------------------------------------------------------------------------------------------#


  #select plate groups
  PlateGroupList <- as.character(unique(data.table(mongo(collection = "htpp_well_trt", url = mongoUrl, verbose = getOption("verbose"))$find())[, pg_id])) #35 plate groups in this study

  #get featurelist
  FeatureList_all <- mongo(collection = "htpp_feature", url = mongoUrl, verbose = getOption("verbose"))$find()
  FeatureList_all <- FeatureList_all %>% select(any_of(c("feature_id", "feature_name_mongo", "feature_name_r", "category_name_r", "cell_type"))) %>% filter(!is.na(feature_id))

  my.cluster <- parallel::makeCluster(
    nThreads,
    type = "PSOCK"
  )

  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  #run tcpl fir on feature level data for each cell type
  for(PG in PlateGroupList){

    print(paste("************", PG, "************"))

    ## get Level5 information
    WellData_all <- data.table(mongo(collection = "htpp_well_norm", url = mongoUrl, verbose = getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
    WellData_all <- WellData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

    ## get CV data
    CV_BMC_all <- data.table(mongo(collection = "cv_bmc", url = mongoUrl, verbose = getOption("verbose"))$find(query = mongoQuery(pg_id = PG)))
    CV_BMC_all <- CV_BMC_all %>% select(pg_id, stype, cell_type, chem_id, cv_noec_dose_level, cv_flag)

    for(cell in unique(WellData_all[, cell_type])){ #TO-DO: This could probably be made into a function to run mclapply() and then use a for loop on the list object for the DB insert

      #keep time
      tic()

      #filter by cell type
      WellData <- as.data.frame(WellData_all[cell_type == cell, ])
      CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])
      FeatureList <- FeatureList_all %>%  dplyr::filter(cell_type == cell)

      # attach
      WellData <- WellData %>% left_join(CV_BMC) %>%
        #modify the information for null chemicals, as they can not be properly processed otherwise
        dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                      cv_flag = ifelse(stype == "null", F, cv_flag))

      #find which columns are NOT feature data
      MetadataColNames <- setdiff(colnames(WellData), FeatureList$feature_name_mongo)

      for(Feature in unique(FeatureList$feature_name_mongo)){

        WellData_1Feature <- WellData %>% select(any_of(MetadataColNames), any_of(Feature)) %>% dplyr::rename(value = Feature)

        #define baseline
        Control <- WellData_1Feature %>% filter(stype == "test sample" & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep > minObjects)
        CompoundList <- unique(WellData_1Feature$chem_id)
        print(CompoundList)

        Control <- Control %>% dplyr::summarise(Mean = mean(value, na.rm = T),
                                                SD = sd(value, na.rm = T))

        concResp <- foreach(x = CompoundList, .combine = 'rbind', .packages = c('tidyr', 'dplyr', 'tcplfit2')) %dopar% {

          #add additional filter to handle DMSO as test chemical
          if(x == "DMSO"){
            Subset <- WellData_1Feature %>% filter(chem_id == x) %>% filter(stype == "test sample")
          }

          #add logic to handle test chemicals that are also ref chems
          for(sample_type in unique(WellData_1Feature$stype[WellData_1Feature$chem_id == x & WellData_1Feature$stype != "vehicle control"])){



            Subset <- WellData_1Feature %>% filter(chem_id == x) %>% filter(stype == sample_type)

            if(is.na(Subset$cv_flag[1])){ #there were too little concentrations to run CV curve fit

              Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                dplyr::summarise(min_conc = min(conc),
                                 max_conc = NA,
                                 n_conc = length(unique(conc)))

              Subset <- Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

            }else if(Subset$cv_flag[1]){

              Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                dplyr::summarise(min_conc = min(conc),
                                 max_conc = NA,
                                 n_conc = 0) %>% distinct()

              Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

            }else{

              Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

              Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                   max_conc = max(conc),
                                                   n_conc = length(unique(conc))) %>%
                select(any_of(c("pg_id", "stype", "cell_type", "chem_id", "min_conc", "max_conc", "n_conc"))) %>% distinct()
            }


            row<-list(pg_id = Metadata$pg_id,
                      stype = Metadata$stype,
                      cell_type = Metadata$cell_type,
                      chem_id = Metadata$chem_id,
                      min_conc = Metadata$min_conc,
                      max_conc = Metadata$max_conc,
                      n_conc = Metadata$n_conc,
                      ctr_mean = Control$Mean,
                      ctr_sd = Control$SD,
                      conc = Subset$conc,
                      resp = Subset$value,
                      bmed = Control$Mean,
                      cutoff = Control$SD,
                      onesd = Control$SD/1.349,
                      approach = "feature",
                      endpoint = Feature)

            newLine <- NULL
            newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"), conthits = T, aicc = F, force.fit = FALSE, bidirectional = TRUE))

            if(is.null(newLine) | class(newLine) == "try-error"){

              newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels = c("cnst")))
            } else if(class(newLine) == "try-error"){

              newLine <- row
            }else{

              rownames(newLine) <- ""
            }

          } #for each stype
          newLine
        } #for each Chem

        concResp<-as_tibble(concResp)


        htpp_tcpl$insert(concResp, auto_unbox=TRUE, na = "null")

        cat("A total of", htpp_tcpl$count(query=mongoQuery(pg_id = PG, cell_type = cell, approach="feature", endpoint=Feature)), "documents of feature level tcpl fits were added to htpp_tcpl for", cell, "cells and", Feature, "feature\n")

      } #for each feature



      cat("A total of", htpp_tcpl$count(query=mongoQuery(pg_id = PG, cell_type = cell, approach="feature")), "documents of feature level tcpl fits were added to htpp_tcpl for", cell, "cells\n")

      toc()

    } #for each cell type

  } #for each plate group
  stopCluster(my.cluster)

  #---------------------------------------------------------------------------------------------#
  # 3. Check if collection is complete
  #---------------------------------------------------------------------------------------------#

  featureCount <- list()
  i <- 1
  for(PG in PlateGroupList){
    featureCount[i] <- length(htpp_well_norm$distinct(key = "chem_id", query = mongoQuery(stype = c("test sample", "null", "reference chemical", "viability positive control"), pg_id = PG)))
    i <- i+1
  }
  featureCount <- do.call(sum, featureCount)

  if(htpp_tcpl$count(query = mongoQuery(approach = "feature")) != featureCount*1300){
    warning(paste("Expected", featureCount*1300, "feature documents in htpp_tcpl, based on htpp_well_norm, instead there are", htpp_tcpl$count(query = mongoQuery(approach = "feature")), "documents."))
  }

  for(cell in unique(WellData_all[, cell_type])){
    cat(htpp_tcpl$count(query = mongoQuery(cell_type = cell, approach = "feature")), "feature-level tcpl fit documents added to htpp_tcpl for", cell, "cells\n")
  }
  } else {
    #---------------------------------------------------------------------------------------------#
    # 1. Connect and check collection and delete if rerun == TRUE
    #---------------------------------------------------------------------------------------------#
    if(file.exists(paste(json_collection_path,"htpp_well_norm.JSON",sep="/"))){
      htpp_well_norm  <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_norm.JSON",sep="/")))
      htpp_well_norm[htpp_well_norm == "NA"] <- NA
    } else {
      stop("htpp_well_norm collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }


    if(file.exists(paste(json_collection_path,"htpp_tcpl.JSON",sep="/"))){
      htpp_tcpl  <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_tcpl.JSON",sep="/")))
      htpp_tcpl[htpp_tcpl == "NA"] <- NA
    } else {
      stop("htpp_tcpl collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }

    if(rerun == TRUE){
      htpp_tcpl <- htpp_tcpl[approach != "feature"]
    }

    if(file.exists(paste(json_collection_path,"htpp_well_trt.JSON",sep="/"))){
      htpp_well_trt  <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_well_trt.JSON",sep="/")))
      htpp_well_trt[htpp_well_trt == "NA"] <- NA
    } else {
      stop("htpp_well_trt collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }

    if(file.exists(paste(json_collection_path,"cv_bmc.JSON",sep="/"))){
      cv_bmc  <-data.table(fromJSON(txt=paste(json_collection_path,"cv_bmc.JSON",sep="/")))
      cv_bmc[cv_bmc == "NA"] <- NA
    } else {
      stop("cv_bmc collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }

    #---------------------------------------------------------------------------------------------#
    # 2. tcplfit for Feature-level data (ca 4 h / plategroup) - FOR EACH CELL TYPE
    #---------------------------------------------------------------------------------------------#

    #To-Do: Explore ways to parallelize this part of the code.

    #select plate groups
    PlateGroupList <- as.character(unique(htpp_well_trt[,pg_id])) #35 plate groups in this study

    #get featurelist
    if(file.exists(paste(json_collection_path,"htpp_feature.JSON",sep="/"))){
      htpp_feature  <-data.table(fromJSON(txt=paste(json_collection_path,"htpp_feature.JSON",sep="/")))
      htpp_feature[htpp_feature == "NA"] <- NA
    } else {
      stop("htpp_feature collection is needed for this function, but is empty.  Check your json_collection_path parameter.")
    }
    FeatureList <- htpp_feature %>% select(any_of(c("feature_id", "feature_name_mongo", "feature_name_r", "category_name_r"))) %>% filter(!is.na(feature_id))

    my.cluster <- parallel::makeCluster(
      nThreads,
      type = "PSOCK"
    )

    #register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)

    #run tcpl fir on feature level data for each cell type
    for(PG in PlateGroupList){

      print(paste("************", PG, "************"))

      ## get Level5 information
      WellData_all <- htpp_well_norm[pg_id==PG,]
      WellData_all <- WellData_all %>% filter(stype %in% c("reference chemical", "test sample", "viability positive control", "null"))

      ## get CV data
      CV_BMC_all <- cv_bmc[pg_id==PG,]
      CV_BMC_all <- CV_BMC_all %>% select(any_of(c("pg_id", "stype", "cell_type", "chem_id", "cv_noec_dose_level", "cv_flag")))

      for(cell in unique(WellData_all[, cell_type])){

        #keep time
        tic()

        #filter by cell type
        WellData <- as.data.frame(WellData_all[cell_type == cell, ])
        CV_BMC <- as.data.frame(CV_BMC_all[cell_type == cell, ])

        # attach
        WellData <- WellData %>% left_join(CV_BMC) %>%
          #modify the information for null chemicals, as they can not be properly processed otherwise
          dplyr::mutate(cv_noec_dose_level = ifelse(stype == "null", 8, cv_noec_dose_level),
                        cv_flag = ifelse(stype == "null", F, cv_flag))

        #find which columns are NOT feature data
        MetadataColNames <- setdiff(colnames(WellData), FeatureList$feature_name_mongo)

        for(Feature in unique(FeatureList$feature_name_mongo)){

          WellData_1Feature <- WellData %>% select(one_of(MetadataColNames), one_of(Feature)) %>% dplyr::rename(value = Feature)

          #define baseline
          Control <- WellData_1Feature %>% filter(stype == "test sample" & dose_level %in% c(1,2) & rel_cell_count > 50 & n_cells_keep > minObjects)
          CompoundList <- unique(WellData_1Feature$chem_id)
          print(CompoundList)

          Control <- Control %>% dplyr::summarise(Mean = mean(value, na.rm = T),
                                                  SD = sd(value, na.rm = T))

          concResp <- foreach(x = CompoundList, .combine = 'rbind', .packages = c('tidyr', 'dplyr', 'tcplfit2')) %dopar% {

            #add additional filter to handle DMSO as test chemical
            if(x == "DMSO"){
              Subset <- WellData_1Feature %>% filter(chem_id == x) %>% filter(stype == "test sample")
            }

            #add logic to handle test chemicals that are also ref chems
            for(sample_type in unique(WellData_1Feature$stype[WellData_1Feature$chem_id == x & WellData_1Feature$stype != "vehicle control"])){



              Subset <- WellData_1Feature %>% filter(chem_id == x) %>% filter(stype == sample_type)

              if(is.na(Subset$cv_flag[1])){ #there were too little concentrations to run CV curve fit

                Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                  dplyr::summarise(min_conc = min(conc),
                                   max_conc = NA,
                                   n_conc = length(unique(conc)))

                Subset = Subset %>% filter(n_cells_keep > minObjects) #keep all concentrations

              }else if(Subset$cv_flag[1]){

                Metadata <- Subset %>% group_by(pg_id, stype, cell_type, chem_id) %>%
                  dplyr::summarise(min_conc = min(conc),
                                   max_conc = NA,
                                   n_conc = 0) %>% distinct()

                Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

              }else{

                Subset <- Subset %>% filter(dose_level <= (cv_noec_dose_level + 1) & n_cells_keep > minObjects)

                Metadata <- Subset %>% dplyr::mutate(min_conc = min(conc),
                                                     max_conc = max(conc),
                                                     n_conc = length(unique(conc))) %>%
                  select(any_of(c("pg_id", "stype", "cell_type", "chem_id", "min_conc", "max_conc", "n_conc"))) %>% distinct()
              }


              row<-list(pg_id = Metadata$pg_id,
                        stype = Metadata$stype,
                        cell_type = Metadata$cell_type,
                        chem_id = Metadata$chem_id,
                        min_conc = Metadata$min_conc,
                        max_conc = Metadata$max_conc,
                        n_conc = Metadata$n_conc,
                        ctr_mean = Control$Mean,
                        ctr_sd = Control$SD,
                        conc = Subset$conc,
                        resp = Subset$value,
                        bmed = Control$Mean,
                        cutoff = Control$SD, #1 nMad was used for the HitDef manuscript
                        onesd = Control$SD/1.349, #we want a BMR of 1*nMad; (1 nMad was also used for the HitDef manuscript)
                        approach = "feature",
                        endpoint = Feature)

              newLine <- NULL
              newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
                                          conthits = T, aicc = F, force.fit = FALSE, bidirectional = TRUE))

              if(is.null(newLine) | class(newLine) == "try-error"){

                newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels = c("cnst")))
              } else if(class(newLine) == "try-error"){

                newLine <- row
              }else{

                rownames(newLine) <- ""
              }

            } #for each stype
            newLine
          } #for each Chem

          concResp<-as_tibble(concResp)

          htpp_tcpl<-rbind(htpp_tcpl, concResp)

          cat("A total of", nrow(htpp_tcpl[pg_id == PG & cell_type == cell & approach=="feature" & endpoint==Feature,]), "documents of feature level tcpl fits were added to htpp_tcpl for", cell, "cells and", Feature, "feature\n")

        } #for each feature



        cat("A total of", nrow(htpp_tcpl[pg_id == PG & cell_type == cell & approach=="feature",]), "documents of feature level tcpl fits were added to htpp_tcpl for", cell, "cells\n")

        toc()

      } #for each cell type

    } #for each plate group
    stopCluster(my.cluster)

    tcplJSON<-toJSON(htpp_tcpl, na="string", digits = 8)
    write(tcplJSON, file=paste(json_collection_path,"htpp_tcpl.JSON",sep="/"))


    #---------------------------------------------------------------------------------------------#
    # 3. Check if collection is complete
    #---------------------------------------------------------------------------------------------#

    featureCount <- list()
    i <- 1
    for(PG in PlateGroupList){
      featureCount[i] <- htpp_well_norm[stype %in% c("test sample", "null", "reference chemical", "viability positive control") & pg_id == PG,] %>% select(chem_id) %>% unique() %>% nrow()
      i <- i+1
    }
    featureCount <- do.call(sum, featureCount)

    if(nrow(htpp_tcpl[approach=="feature",]) != featureCount*1300){
      warning(paste("Expected", featureCount*1300, "feature documents in htpp_tcpl, based on htpp_well_norm, instead there are", nrow(htpp_tcpl[approach=="feature"]), "documents."))
    }

    for(cell in unique(WellData_all[, cell_type])){
      cat(nrow(htpp_tcpl[approach=="feature" & cell_type==cell,]), "feature-level tcpl fit documents added to htpp_tcpl for", cell, "cells\n")
    }

}
}

