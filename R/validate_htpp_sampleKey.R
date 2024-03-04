#' Reformatting the sample key into a more machine-readable form for later processing and checking it for errors
#'
#' @param SampleKey File name, path to file, of the sample key file (in csv format) being used
#' @param skipped_tests The names of QC tests you want to skip
#' @param max_dose_level The maximum dose level used
#' @param dataFrame Boolean, if TRUE, will treat SampleTable as data.frame/data.table input and not read in file; default is FALSE
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
#' @import foreach
#'
#' @return A reformatted sample key if there were no errors, or a list of errors if there were
#'
#' @export validate_htpp_sampleKey
#'
#' @examples
#' sample_key <- data.table::fread(
#' file = system.file("extdata", "example_sampleKey.csv", package = "htpp.pl"),
#' sep = ",")
#'
#' validated_sample_key <- validate_htpp_sampleKey(SampleKey = sample_key,
#' max_dose_level = 8, skipped_tests = "QCV_0", dataFrame = TRUE)
#'

validate_htpp_sampleKey<-function(SampleKey, skipped_tests=c(), max_dose_level=8, dataFrame = FALSE){

  if(dataFrame == FALSE){
    SampleTable <- data.table(read_csv(SampleKey))
  }else{
    SampleTable <- data.table(SampleKey)
  }

  required_cols = c("replicate_num", "cell_type", "culture_id", "pg_id", "doseplate_id", "stype", "chem_id", "dtxsid", "casrn", "chem_name", "dose_level", "conc", "conc_unit", "sample_id", "plate_id", "well_id", "trt_name", "assay", "qc_flag", "qc_flag_description")

  if(all(required_cols %in% colnames(SampleTable)) == FALSE){
    stop("SampleKey is missing the following required columns: ", paste(setdiff(required_cols, colnames(SampleTable)), collapse = ", "))
  }

  NAcheck<-unique(SampleTable[is.na(chem_id), .(chem_name, stype)])
  if(nrow(NAcheck)>0){
    warning("Warning, NAs in Chemical ID field. Check that vehicle contorls, viability positive controls, and reference chemicals have a defined chem_id.")
    message(paste(capture.output(head(NAcheck, 5)), collapse="\n"))
  }

  SampleTable$qc_flag <- as.character(SampleTable$qc_flag)

  Table <- SampleTable %>%
    dplyr::mutate(plate_id = plate_id,
           well_id = well_id,
           sample_id = paste0(plate_id, "_", well_id),
           replicate_num = replicate_num,
           cell_type = cell_type,
           culture_id = culture_id,
           doseplate_id = doseplate_id,
           stype = ifelse(stype=="viability postive control", "viability positive control", stype),
           pg_id = as.character(pg_id),
           chem_id = ifelse(is.na(chem_id),chem_name, chem_id),
           chem_name = chem_name,
           dtxsid = dtxsid,
           casrn = casrn,
           conc = signif(conc, 3),
           conc_unit = conc_unit,
           dose_level = as.integer(dose_level),
           qc_flag = qc_flag,
           qc_flag_description = ifelse(qc_flag_description == "NA", NA, qc_flag_description),
           trt_name =trt_name,
           assay = paste(assay)) %>%
    select(replicate_num, cell_type, culture_id, pg_id, doseplate_id, stype, chem_id, dtxsid, casrn, chem_name, dose_level, conc, conc_unit,
           sample_id, plate_id, well_id, trt_name, assay, qc_flag, qc_flag_description)
  PlateGroups <- unique(SampleTable[, c("plate_id", "pg_id")])

  tests <- fromJSON(txt=system.file("extdata", "sampleKeyTests.json", package = "htpp.pl"))
  #check for sample_id duplication

  extra_cols = c()
  required_cols = c("replicate_num", "cell_type", "culture_id", "pg_id", "doseplate_id", "stype", "chem_id", "dtxsid", "casrn", "chem_name", "dose_level", "conc", "conc_unit", "sample_id", "plate_id", "well_id", "trt_name", "assay", "qc_flag", "qc_flag_description")

  dbCols <- union(c("replicate_num", "cell_type", "culture_id", "pg_id", "doseplate_id", "stype", "chem_id", "dtxsid", "casrn", "chem_name", "dose_level", "conc", "conc_unit", "sample_id", "plate_id", "well_id", "trt_name", "assay", "qc_flag", "qc_flag_description"), extra_cols)
  results <- list()


  if (!"CIT_0" %in% skipped_tests){
    message(tests$CIT_0)
    if(!all(required_cols %in% colnames(Table))) {
      newelem <- paste("Missing the following columns: ", paste(setdiff(dbCols, colnames(Table)), collapse=", "), "(CIT_0)")
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$CIT_0))

  if (!"UCT_0" %in% skipped_tests){
    message(tests$UCT_0)
    if(!all(colnames(Table) %in% dbCols)) {
      newelem <- paste("Table contains unexpected columns: ", paste(setdiff(colnames(Table), dbCols), collapse=", "), "(UCT_0)")
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$UCT_0))

  if (!"NAS_0" %in% skipped_tests){
    message(tests$NAS_0)
    if (sum(is.na(Table[,sample_id])) != 0){
      newelem <- "Table contains NA sample_ids (NAS_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAS_0))

  if (!"DST_0" %in% skipped_tests){
    message(tests$DST_0)
    if (sum(duplicated(Table[,sample_id])) != 0){
      newelem <- "Table contains duplicated sample_ids (DST_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$DST_0))

  if (!"SIF_0" %in% skipped_tests){
    message(tests$SIF_0)
    if (all(grepl("^TC[0-9]{7,8}_[A-Z][0-9]{2}$", Table[,sample_id])) != TRUE){
      newelem <- "Table contains sample_ids that don't follow the ^TC[0-9]{8}_[A-Z][0-9]{2} pattern (SIF_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$SIF_0))

  if (!"SPR_0" %in% skipped_tests && 'plate_id' %in% colnames(Table) && 'well_id' %in% colnames(Table)){
    message(tests$SPR_0)
    if (all(Table[,sample_id] == paste(Table[,plate_id], Table[,well_id], sep="_")) != TRUE){
      newelem <- "Table contains sample_ids that don't correspond to the plate_id_well_id (SPR)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$SPR_0))

  if (!"NAP_0" %in% skipped_tests && 'plate_id' %in% colnames(Table)){
    message(tests$NAP_0)
    if (sum(is.na(Table[,plate_id])) != 0){
      newelem <- "Table contains plate_ids that are NA (NAP_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAP_0))

  if (!"PIF_0" %in% skipped_tests && 'plate_id' %in% colnames(Table)){
    message(tests$PIF_0)
    if (all(grepl("^TC[0-9]{7,8}$", Table[,plate_id])) != TRUE){
      newelem <- "Table contains plate_ids that don't follow the TC[0-9]{8} format (PIF_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$PIF_0))

  if (!"WIF_0" %in% skipped_tests && 'well_id' %in% colnames(Table)){
    message(tests$WIF_0)
    if (all(grepl("^[A-Z][0-9]{2}$", Table[,well_id])) != TRUE){
      newelem <- "Table contains well_id that don't follow the TC[0-9]{2} format (WIF_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$WIF_0))

  if (!"WPC_0" %in% skipped_tests && 'plate_id' %in% colnames(Table) && 'well_id' %in% colnames(Table)){
    message(tests$WPC_0)
    if (all(table(Table[,well_id]) <= length(unique(Table[,plate_id]))) != TRUE){
      newelem <- "less well_ids found than plate_ids (WPC_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$WPC_0))

  if (!"NAT_0" %in% skipped_tests){
    message(tests$NAT_0)
    if (sum(is.na(Table[,trt_name])) != 0){
      newelem <- "Table contains trt_name that are NA (NAT_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAT_0))

  if (!"TRR_0" %in% skipped_tests){
    message(tests$TRR_0)
    if (nrow(Table[stype=="QC sample"])!=0 && length(unique(table(Table[stype=="QC sample", trt_name]))) != 1){
      newelem <- "Table contains QC samples for which there isn't a fixed number of replicates per trt_name (TRR_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$TRR_0))

  if (!"TRS_0" %in% skipped_tests){
    message(tests$TRS_0)

    # Every trt_name should correspond to exactly one stype
    all_trt_name <- unique(Table[,trt_name])
    all_trt_stype <- foreach(tn = all_trt_name, .combine='c') %do% {
      length(unique(Table[trt_name == tn, stype]))
    }
    names(all_trt_stype) <- all_trt_name

    if (all(all_trt_stype != 1)){
      newelem <- "Each trt_name doesn't correspond to exactly one stype (TRS_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$TRS_0))

  if (!"TRP_0" %in% skipped_tests && 'plate_id' %in% colnames(Table)){
    # For test samples, there should be 1 replicate of each trt_name on each plate
    message(tests$TRP_0)
    for(pid in unique(Table[, plate_id])) {
      pid_test_wells <- Table[(stype == "test sample") & (plate_id == pid), ]
      if (any(table(pid_test_wells[, trt_name]) != 1)){
        newelem <- "Some test samples have more than 1 replicate of each trt_name on each plate (TRP_0)"
        results <- c(results, newelem)
        break
      }
    }
  }
  else
    message(paste("skipping", tests$TRP_0))

  if (!"QCN_0" %in% skipped_tests){
    # qc_flag - No NAs, Should be in set {OK, CELL_VIABILITY}
    message(tests$QCN_0)
    if (sum(is.na(Table[,qc_flag])) != 0){
      newelem <- "qc_flag should not be NA (QCN_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$QCN_0))

  if (!"QCV_0" %in% skipped_tests){
    message(tests$QCV_0)
    if (all(Table[,qc_flag] %in% c("OK", "CELL_VIABILITY", "DOSEPLATE_FAIL", "DISPENSE_FAIL", "SINGLE_REP")) != TRUE){
      newelem <- "qc_flag must be one of those values: OK, CELL_VIABILITY, DOSEPLATE_FAIL, DISPENSE_FAIL or SINGLE_REP (QCV_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$QCV_0))

  if (!"NAS_1" %in% skipped_tests && ('stype' %in% colnames(Table))){
    message(tests$NAS_1)
    if (sum(is.na(Table[,stype])) != 0){
      newelem <- "stype can't be NA (NAS_1)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAS_1))

  if (!"STV_0" %in% skipped_tests && ('stype' %in% colnames(Table))){
    message(tests$STV_0)
    if (all(Table[,stype] %in% c("test sample", "viability positive control", "vehicle control", "reference chemical")) != TRUE){
      newelem <- "stype must be one of those values: test sample, viability positive control, vehicle control, reference chemical (STV_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$STV_0))

  if (!"NAR_0" %in% skipped_tests && ('rna_src' %in% colnames(Table))){
    message(tests$NAR_0)
    if (sum(is.na(Table[,rna_src]))!=0){
      newelem <- "rna_src can't be NA NAR_0"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAR_0))

  if (!"CON_0" %in% skipped_tests && ('conc' %in% colnames(Table))){
    message(tests$CON_0)
    if (!is.numeric(Table[,conc])){
      newelem <- "conc should be numeric (CON_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$CON_0))

  if (!"COC_0" %in% skipped_tests && 'conc' %in% colnames(Table) && 'conc_unit' %in% colnames(Table)){
    message(tests$COC_0)
    # conc_unit - NA status should match conc
    if (all(is.na(Table[,conc_unit]) == is.na(Table[,conc])) != TRUE){
      newelem <- "conc_unit - NA status should match conc (COC_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$COC_0))

  if (!"DMC_0" %in% skipped_tests && 'dose_level' %in% colnames(Table) && 'conc' %in% colnames(Table)){

    message(tests$DMC_0)
    # dose_level - NA status should match conc, should be integer in range 0:8, should be 0 for vehicle controls and BL DMSOs only
    if (all(is.na(Table[,conc]) == is.na(Table[,dose_level])) != TRUE){
      newelem <- "dose_level - NA status should match conc (DMC_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$DMC_0))

  if (!"DLI_0" %in% skipped_tests && 'dose_level' %in% colnames(Table)){
    message(tests$DLI_0)
    if (!is.integer(Table[,dose_level])){
      newelem <- "dose_level should be integer (DLI_0)"
      results <- c(results, newelem)
    }
    #message(Table[,dose_level])
  }
  else
    message(paste("skipping", tests$DLI_0))

  if (!"DLR_0" %in% skipped_tests && 'dose_level' %in% colnames(Table)){
    message(paste0(tests$DLR_0, max_dose_level, "]"))
    if (all(Table[!is.na(dose_level), dose_level] %in% 0:max_dose_level) != TRUE){
      newelem <- paste0("dose_level should be in range [0:", max_dose_level, "] (DLR_0)")
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$DLR_0))

  if (!"DL0_0" %in% skipped_tests && 'stype' %in% colnames(Table) && 'chem_id' %in% colnames(Table) && 'dose_level' %in% colnames(Table)){
    message(tests$DL0_0)
    if (all(Table[(((stype == "vehicle control") | (chem_id == "DMSO")) & !is.na(dose_level)), dose_level] == 0) != TRUE){
      newelem <- "dose_level should be 0 for vehicle controls and BL DMSOs only (DL0_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$DL0_0))

  if (!"DCC_0" %in% skipped_tests && 'stype' %in% colnames(Table) && 'chem_id' %in% colnames(Table) && 'dose_level' %in% colnames(Table)){
    message(tests$DCC_0)
    dose_warning=FALSE
    dose_conc_warning=FALSE
    conc_increasing_warning=FALSE
    # For each chem_id, each dose level should have the same conc and be monotonic increasing
    for(chem in unique(Table[stype %in% c("test sample", "reference chemical"), chem_id])) {
      # Get the dose levels for this chem_id
      chem_levels <- sort(unique(Table[!is.na(dose_level) & (chem_id == chem), dose_level]))

      if (length(chem_levels) <= 0 && dose_warning==FALSE){
        newelem <- "For each chem_id, each dose level don't have the same conc or are monotonic increasing (DCC_0)"
        results <- c(results, newelem)
        dose_warning=TRUE
      }


      if (max(chem_levels) != length(chem_levels) && dose_warning==FALSE){
        newelem <- "For each chem_id, each dose level don't have the same conc or are monotonic increasing (DCC_0)"
        results <- c(results, newelem)
        dose_warning=TRUE
      }

      # There should be a single conc per dose level, and should be in increasing order:
      chem_concs <- foreach(dose = chem_levels, .combine='c') %do% {
        dose_conc <- unique(Table[(chem_id == chem) & (dose_level == dose), conc])
        if (length(dose_conc) != 1 && dose_conc_warning==FALSE){
          newelem <- "Some dose have more than 1 conc (DCC_0)"
          results <- c(results, newelem)
          dose_conc_warning=TRUE
        }
        return(dose_conc)
      }
      if (all(chem_concs == sort(chem_concs, decreasing = F)) != TRUE && conc_increasing_warning==FALSE){
        newelem <- "conc should be in increasing level (DCC_0)"
        results <- c(results, newelem)
        conc_increasing_warning=TRUE
      }
    }
  }
  else
    message(paste("skipping", tests$DCC_0))

  if (!"CPC_0" %in% skipped_tests && 'stype' %in% colnames(Table) && 'chem_id' %in% colnames(Table)){
    # For each test chemical - each chem_id should be on a single pg_id, and have no more than 8 * 4 samples
    chem_dose_rep_pb_detected = FALSE
    chem_id_pb_detected = FALSE
    chem_pgs_pb_detected = FALSE

    message(tests$CPC_0)

    for(chem in unique(Table[stype == "test sample", chem_id])) {
      chem_pgs <- unique(Table[chem_id == chem, pg_id])
      if (chem_pgs_pb_detected ==FALSE && length(chem_pgs) != 1){
        newelem <- "Some test chemical have chem_id that is not on a single pg_id (CPC_0)"
        results <- c(results, newelem)
        chem_pgs_pb_detected = TRUE
      }

      chem_wells <- Table[chem_id == chem, ]
      if (chem_id_pb_detected == FALSE && nrow(chem_wells) > (8*4)){
        newelem <- "Some test chemicals have more than 8 * 4 samples (CPC_0)"
        results <- c(results, newelem)
        chem_id_pb_detected = TRUE
      }

      # Every dose_level, replicate_num should be unique
      if ("replicate_num" %in% colnames(Table)){
        chem_dose_rep <- paste(chem_wells$dose_level, chem_wells$replicate_num, "_")
        if (chem_dose_rep_pb_detected == FALSE && sum(duplicated(chem_dose_rep)) != 0){
          newelem <- "For some dose level, replicate_num is not unique (CPC_0)"
          results <- c(results, newelem)
          chem_dose_rep_pb_detected = TRUE
        }
        if (chem_id_pb_detected && chem_dose_rep_pb_detected && chem_pgs_pb_detected)
          break
      }

    }
  }
  else
    message(paste("skipping", tests$CPC_0))

  if (!"PIC_0" %in% skipped_tests && 'pg_id' %in% colnames(Table)){
    message(tests$PIC_0)
    # pg_id - Character type2
    if (!is.character(Table[,pg_id])){
      newelem <- "pg_id should be of type Character (PIC_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$PIC_0))

  if (!"NAP_1" %in% skipped_tests && 'pg_id' %in% colnames(Table)){
    message(tests$NAP_1)

    if (sum(is.na(Table[,pg_id])) != 0){
      newelem <- "pg_id should not be NA (NAP_1)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAP_1))

  if (!"PGC_0" %in% skipped_tests && 'plate_id' %in% colnames(Table) && 'pg_id' %in% colnames(Table) && 'stype' %in% colnames(Table)){
    message(tests$PGC_0)
    # Make sure plates are balanced within plate_group, and all contain the same chemicals
    plt_chem_set <- foreach(plt = unique(Table[,plate_id]), .combine='c') %do% {
      plt_chems <- unique(Table[(plate_id == plt) & (stype == "test sample"), chem_id])
      return(paste(sort(plt_chems), collapse=","))
    }
    names(plt_chem_set) <- unique(Table[,plate_id])
    plates_balanced=FALSE
    plates_contain_chem=FALSE
    plates_contain_only=FALSE
    for(pg in unique(Table[,pg_id])) {
      pg_plates <- Table[pg_id == pg, plate_id]
      if (plates_balanced==FALSE && length(unique(table(pg_plates))) != 1){
        newelem <- "plates are not balanced within plate group (PGC_0)"
        results <- c(results, newelem)
        plates_balanced=TRUE
      }

      pg_plates <- unique(pg_plates)
      if (plates_contain_chem==FALSE && all(pg_plates %in% names(plt_chem_set)) != TRUE){
        newelem <- "some plates do not contain the same chemical (PGC_0)"
        results <- c(results, newelem)
        plates_contain_chem=TRUE
      }

      if (plates_contain_only==FALSE && length(unique(plt_chem_set[pg_plates])) != 1){
        newelem <- "some plates do not contain the same chemical (PGC_0)"
        results <- c(results, newelem)
        plates_contain_only=TRUE
      }
    }
  }
  else
    message(paste("skipping", tests$PGC_0))


  if (!"NAR_0" %in% skipped_tests && 'replicate_num' %in% colnames(Table)){
    message(tests$NAR_0)
    if (sum(is.na(Table[, replicate_num])) != 0){
      newelem <- "replicate_num should not be NA (NAR_0)"
      results <- c(results, newelem)
    }
  }
  else
    message(paste("skipping", tests$NAR_0))


  message("\n\n")
  message(paste("Check that all wells have metadata"))
  message(paste0(capture.output(table(Table$plate_id)), collapse = "\n"))

  message("\n\n")
  message(paste("\n Check that there are the right number of vehicle control, test sample, reference chemical and viability postive control wells"))
  message(paste0(capture.output(table(Table$stype)/4), collapse = "\n"))

  if(length(results)==0){
    return(Table)
  }else{
    return(results)
  }

}
