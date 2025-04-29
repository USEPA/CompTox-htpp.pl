
#' Adding chemical id and sample data to htpp_well_trt, htpp_chem
#'
#' @param SampleKey dataframe: The sample key with chemical info and metadata
#' @param rerun boolean: Whether you want to clear the mongo database as you go and refill it; false by default
#' @param replace boolean: Whether you want to replace existing records in the mongo database; false by default
#' @param mongoUrl characters string A mongoUrl with credentials to access the database
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
#' @export generate_htppWellTrt_htppChem
#'

generate_htppWellTrt_htppChem <- function(SampleKey, mongoUrl, rerun=FALSE, replace=TRUE){

  htpp_well_trt <- mongo(collection = "htpp_well_trt", url = mongoUrl, verbose = getOption("verbose"))
  htpp_well_trt$count()

  if(rerun == TRUE){
    htpp_well_trt$drop() #delete everything that was in it
  }
  keysmith <- htpp_well_trt$distinct("sample_id")


  if(replace == FALSE) {
    for(i in 1:dim(SampleKey)[1]){
      if(SampleKey[i,]$sample_id %in% keysmith){
        next
      }else{
        htpp_well_trt$insert(SampleKey[i,], auto_unbox=TRUE, na = "null")
      }
    }
  } else{
    #insert into collection
    for(i in 1:dim(SampleKey)[1]){
      htpp_well_trt$remove(query = mongoQuery(sample_id = SampleKey[i,]$sample_id)) #remove existing record
      htpp_well_trt$insert(SampleKey[i,], auto_unbox=TRUE, na = "null") #replace existing record
    }
  }




  #Use chem info that is part of the sample key
  ChemInfo <- SampleKey

  #columns needed for htpp_chem: stype, dtxsid, casrn, chem_name
  ChemInfo <- ChemInfo[, c("stype", "chem_name", "chem_id","dtxsid", "casrn")] #keep sample type for reference
  ChemInfo <- unique(ChemInfo)



  htpp_chem <- mongo(collection = "htpp_chem", url = mongoUrl, verbose = getOption("verbose"))

  if(rerun == TRUE){
    htpp_chem$drop()
  }

  ChemInfo <- dplyr::mutate(ChemInfo, chem_id=ifelse(is.na(chem_id), chem_name, chem_id))
  keysmith1 <- htpp_chem$distinct("chem_id")


  if(replace == FALSE) {
    for(i in 1:dim(ChemInfo)[1]){
      if(ChemInfo[i,]$chem_id %in% keysmith1){
        next
      }else{
        htpp_chem$insert(ChemInfo[i,], auto_unbox=TRUE)
      }
    }
  } else{
    #insert into collection
    for(i in 1:dim(ChemInfo)[1]){
      htpp_chem$remove(query = mongoQuery(chem_id = ChemInfo[i,]$chem_id)) #remove existing record
      htpp_chem$insert(ChemInfo[i,], auto_unbox=TRUE) #replace record
    }
  }

  if(htpp_well_trt$count() != length(SampleKey$sample_id)){
    message(paste("sample_id counts do not match.  Sample counts in SampleKey=", length(SampleKey$sample_id), "and sample counts in Mongo database=",htpp_well_trt$count()))
  }

  if(length(unique((ChemInfo$chem_id))) != htpp_chem$count()){
    message(paste("Chemical counts do not match.  Chemical counts in ChemInfo =", length(ChemInfo$chem_id), "and chemical counts in Mongo database =", htpp_chem$count()))
  }
}

#' Inserts feature and category data into mongo collection htpp_feature and htpp_category
#'
#' @param inputPath  character string: Can either be a truncated path, or a full path to a HTPP data file. If it is truncated, the function will rebuild a full path using file_path
#' @param PlateID character string: The PlateID for the plate being analyzed
#' @param mongoUrl character string: The database where the collections will be stored
#' @param file_path character string: The path to where the input file is located
#' @param rerun boolean:  Whether to delete and reinsert into both collections; false by default
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
#' @export generate_htppFeature_htppCategory

generate_htppFeature_htppCategory <- function(inputPath, PlateID, mongoUrl, file_path="", rerun=FALSE){
  if (nchar(file_path) > 0){
    Table1 <- suppressMessages(read_delim(paste(file_path, inputPath, sep = "/"), delim = "\t", col_names = T, skip = 9, n_max = 1000) )
  }else{
    Table1 <- suppressMessages(read_delim(inputPath, delim = "\t", col_names = T, skip = 9, n_max = 1000) )
  }
  Table1 <- as_tibble(Table1) %>% select(-starts_with("X1")) %>%
    dplyr::mutate(plate_id = PlateID,
                  well_id = Well3Digit(paste0(LETTERS[Row], Column)),
                  sample_id = paste0(plate_id, "_", well_id),
                  sample_id_field = paste0(sample_id, "_f", Field),
                  cell_keep = T,
                  Position_both_Centroid_Distance = 0)
  FeatureList <- tibble(feature_name_harmony = colnames(Table1)) %>%
    dplyr::mutate(feature_name_r = feature_name_harmony)


  #rename some R names
  FeatureList$feature_name_r[which(FeatureList$feature_name_harmony=="Field")] = "field"
  FeatureList$feature_name_r[which(FeatureList$feature_name_harmony=="Object No")] = "object_no"
  FeatureList$feature_name_r[which(FeatureList$feature_name_harmony=="X")] = "position_x"
  FeatureList$feature_name_r[which(FeatureList$feature_name_harmony=="Y")] = "position_y"
  FeatureList$feature_name_r[which(FeatureList$feature_name_harmony=="Bounding Box")] = "bounding_box"

  #remove the prefix
  FeatureList$feature_name_r <- str_replace(FeatureList$feature_name_r, "Cells Non-Border - ", "")

  #remplace empty spaces
  FeatureList$feature_name_r <- str_replace_all(FeatureList$feature_name_r, " ", "_")

  #remove the unit for the two features with um2
  i=which(grepl("Shape", FeatureList$feature_name_r) & grepl("_Area_", FeatureList$feature_name_r))
  FeatureList$feature_name_r[i] <-str_sub(FeatureList$feature_name_r[i], 1, -7)

  #remove the unit for the four features with um
  i=which(grepl("Shape", FeatureList$feature_name_r) & xor(grepl("_Width_", FeatureList$feature_name_r),  grepl("_Length", FeatureList$feature_name_r) ) )
  FeatureList$feature_name_r[i] <-str_sub(FeatureList$feature_name_r[i], 1, -6)

  #remove the unit for 7 features with um
  i=which(grepl("Position", FeatureList$feature_name_r) & (grepl("_Distance", FeatureList$feature_name_r) |  grepl("_Centroid", FeatureList$feature_name_r) ) )
  FeatureList$feature_name_r[i] <-str_sub(FeatureList$feature_name_r[i], 1, -6)

  FeatureList <- FeatureList %>% dplyr::mutate(column_type = NA)

  FeatureList$column_type[which(FeatureList$feature_name_r %in% c("plate_id", "well_id", "sample_id", "sample_id_field",
                                                                  "field", "object_no", "position_x", "position_y", "bounding_box", "cell_keep"))] <- "M"

  for(Channel in c("^Shape", "^Position", "^DNA", "^RNA", "^ER", "^AGP", "Mito")){
    B <- grepl(pattern=Channel, FeatureList$feature_name_r)
    FeatureList$column_type[B] <-"F"
  }

  FeatureList <- FeatureList %>% dplyr::mutate(category_name_r = NA)

  ###### Position
  FeatureList$category_name_r[(grepl(pattern="^Position_", FeatureList$feature_name_r) &
                                 !(grepl(pattern="?m", FeatureList$feature_name_r) | grepl(pattern="in_Image", FeatureList$feature_name_r)))] <- "Position"

  ###### Shape
  FeatureList$category_name_r[grepl(pattern="Shape", FeatureList$feature_name_r)] <- "Shape"

  ###### Mito
  #STAR
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Symmetry",              FeatureList$feature_name_r)] <- "Mito_Symmetry_Cells"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Threshold_Compactness", FeatureList$feature_name_r)] <- "Mito_Compactness_Cells"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Axial",                 FeatureList$feature_name_r)] <- "Mito_Axial_Cells"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Radial",                FeatureList$feature_name_r)] <- "Mito_Radial_Cells"
  #Intensity
  FeatureList$category_name_r[grepl(pattern="Mito_Ring_Intensity",      FeatureList$feature_name_r)] <- "Mito_Intensity_Ring"
  FeatureList$category_name_r[grepl(pattern="Mito_Cytoplasm_Intensity", FeatureList$feature_name_r)] <- "Mito_Intensity_Cytoplasm"
  #Texture
  FeatureList$category_name_r[grepl(pattern="Mito_Ring_Texture",      FeatureList$feature_name_r)] <- "Mito_Texture_Ring"
  FeatureList$category_name_r[grepl(pattern="Mito_Cytoplasm_Texture", FeatureList$feature_name_r)] <- "Mito_Texture_Cytoplasm"

  ###### AGP
  #STAR
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Symmetry",              FeatureList$feature_name_r)] <- "AGP_Symmetry_Cells"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Threshold_Compactness", FeatureList$feature_name_r)] <- "AGP_Compactness_Cells"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Axial",                 FeatureList$feature_name_r)] <- "AGP_Axial_Cells"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Radial",                FeatureList$feature_name_r)] <- "AGP_Radial_Cells"
  #Intensity
  FeatureList$category_name_r[grepl(pattern="AGP_Ring_Intensity",      FeatureList$feature_name_r)] <- "AGP_Intensity_Ring"
  FeatureList$category_name_r[grepl(pattern="AGP_Cytoplasm_Intensity", FeatureList$feature_name_r)] <- "AGP_Intensity_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="AGP_Membrane_Intensity",  FeatureList$feature_name_r)] <- "AGP_Intensity_Membrane"
  #Texture
  FeatureList$category_name_r[grepl(pattern="AGP_Ring_Texture",      FeatureList$feature_name_r)] <- "AGP_Texture_Ring"
  FeatureList$category_name_r[grepl(pattern="AGP_Cytoplasm_Texture", FeatureList$feature_name_r)] <- "AGP_Texture_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="AGP_Membrane_Texture",  FeatureList$feature_name_r)] <- "AGP_Texture_Membrane"

  ###### ER
  #STAR
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Symmetry",              FeatureList$feature_name_r)] <- "ER_Symmetry_Cells"
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Threshold_Compactness", FeatureList$feature_name_r)] <- "ER_Compactness_Cells"
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Axial",                 FeatureList$feature_name_r)] <- "ER_Axial_Cells"
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Radial",                FeatureList$feature_name_r)] <- "ER_Radial_Cells"
  #Intensity
  FeatureList$category_name_r[grepl(pattern="ER_Ring_Intensity",      FeatureList$feature_name_r)] <- "ER_Intensity_Ring"
  FeatureList$category_name_r[grepl(pattern="ER_Cytoplasm_Intensity", FeatureList$feature_name_r)] <- "ER_Intensity_Cytoplasm"
  #Texture
  FeatureList$category_name_r[grepl(pattern="ER_Ring_Texture",      FeatureList$feature_name_r)] <- "ER_Texture_Ring"
  FeatureList$category_name_r[grepl(pattern="ER_Cytoplasm_Texture", FeatureList$feature_name_r)] <- "ER_Texture_Cytoplasm"

  ###### RNA
  #STAR
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Morph_STAR_Symmetry",              FeatureList$feature_name_r)] <- "RNA_Symmetry_Nuclei"
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Morph_STAR_Threshold_Compactness", FeatureList$feature_name_r)] <- "RNA_Compactness_Nuclei"
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Morph_STAR_Axial",                 FeatureList$feature_name_r)] <- "RNA_Axial_Nuclei"
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Morph_STAR_Radial",                FeatureList$feature_name_r)] <- "RNA_Radial_Nuclei"
  #Intensity
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Intensity",      FeatureList$feature_name_r)] <- "RNA_Intensity_Nuclei"
  #Texture
  FeatureList$category_name_r[grepl(pattern="RNA_Nuclei_Texture",      FeatureList$feature_name_r)] <- "RNA_Texture_Nuclei"

  ###### DNA
  #STAR
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Morph_STAR_Symmetry",              FeatureList$feature_name_r)] <- "DNA_Symmetry_Nuclei"
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Morph_STAR_Threshold_Compactness", FeatureList$feature_name_r)] <- "DNA_Compactness_Nuclei"
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Morph_STAR_Axial",                 FeatureList$feature_name_r)] <- "DNA_Axial_Nuclei"
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Morph_STAR_Radial",                FeatureList$feature_name_r)] <- "DNA_Radial_Nuclei"
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Radial",                 FeatureList$feature_name_r)] <- "DNA_Radial_Cells"
  #Intensity
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Intensity",      FeatureList$feature_name_r)] <- "DNA_Intensity_Nuclei"
  #Texture
  FeatureList$category_name_r[grepl(pattern="DNA_Nuclei_Texture",      FeatureList$feature_name_r)] <- "DNA_Texture_Nuclei"

  ##### Profile
  #Mito
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Profile_1",  FeatureList$feature_name_r)] <- "Mito_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Profile_2",  FeatureList$feature_name_r)] <- "Mito_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Profile_3",  FeatureList$feature_name_r)] <- "Mito_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Profile_4",  FeatureList$feature_name_r)] <- "Mito_Profile_Nuclei"
  FeatureList$category_name_r[grepl(pattern="Mito_Cells_Morph_STAR_Profile_5",  FeatureList$feature_name_r)] <- "Mito_Profile_Nuclei"
  #AGP
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Profile_1",   FeatureList$feature_name_r)] <- "AGP_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Profile_2",   FeatureList$feature_name_r)] <- "AGP_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Profile_3",   FeatureList$feature_name_r)] <- "AGP_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Profile_4",   FeatureList$feature_name_r)] <- "AGP_Profile_Nuclei"
  FeatureList$category_name_r[grepl(pattern="AGP_Cells_Morph_STAR_Profile_5",   FeatureList$feature_name_r)] <- "AGP_Profile_Nuclei"
  #ER
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Profile_1",    FeatureList$feature_name_r)] <- "ER_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Profile_2",    FeatureList$feature_name_r)] <- "ER_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="ER_Cells_Morph_STAR_Profile_3",    FeatureList$feature_name_r)] <- "ER_Profile_Cytoplasm"
  #RNA
  FeatureList$category_name_r[grepl(pattern="RNA_Cells_Morph_STAR_Profile_4",   FeatureList$feature_name_r)] <- "RNA_Profile_Nuclei"
  FeatureList$category_name_r[grepl(pattern="RNA_Cells_Morph_STAR_Profile_5",   FeatureList$feature_name_r)] <- "RNA_Profile_Nuclei"
  FeatureList$category_name_r[(grepl(pattern="RNA_Nuclei_Morph_STAR_Profile",   FeatureList$feature_name_r) &
                                 !grepl(pattern="SER", FeatureList$feature_name_r))] <- "RNA_Profile_Nuclei"
  #(The SER endpoints are identical for Nuclei and Cell, therefore they need to be removed once)

  #DNA
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Profile_1",     FeatureList$feature_name_r)] <- "DNA_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Profile_2",     FeatureList$feature_name_r)] <- "DNA_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Profile_3",     FeatureList$feature_name_r)] <- "DNA_Profile_Cytoplasm"
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Profile_4",     FeatureList$feature_name_r)] <- "DNA_Profile_Nuclei"
  FeatureList$category_name_r[grepl(pattern="DNA_Cells_Morph_STAR_Profile_5",     FeatureList$feature_name_r)] <- "DNA_Profile_Nuclei"
  FeatureList$category_name_r[(grepl(pattern="DNA_Nuclei_Morph_STAR_Profile",   FeatureList$feature_name_r) &
                                 !grepl(pattern="SER", FeatureList$feature_name_r))] <- "DNA_Profile_Nuclei"

  #(The SER endpoints are identical for Nuclei and Cell, therefore they need to be removed once)
  dummy <- FeatureList %>% filter(column_type =="F" & !is.na(category_name_r)) %>%
    arrange(desc(column_type),category_name_r, feature_name_r) %>%
    dplyr::mutate(feature_id = row_number(), feature_name_mongo = paste0("f_", feature_id))

  FeatureList <- FeatureList %>% left_join(dummy) %>%
    dplyr::mutate(feature_name_mongo = ifelse(column_type == "M", feature_name_r, feature_name_mongo)) %>%
    arrange(desc(column_type), category_name_r,  feature_name_r)

  write_csv(FeatureList, file = "./htpp_feature.csv")
  range(FeatureList$feature_id, na.rm=T)

  htpp_feature <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))

  #delete and reinsert if rerun
  if(rerun == TRUE){
    htpp_feature$drop()
  }

  htpp_feature$insert(FeatureList, na = "null")
  if(htpp_feature$count() != 1410){
    warning(paste("Incorrect number of feature documents.  Expected 1410 features, instead there are", htpp_feature$count(), "features in collection."))
  }

  List <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))$find()

  List2 <- List %>% filter(!is.na(feature_id)) %>%
    group_by(category_name_r) %>%
    dplyr::summarise(n_features = n_distinct(feature_id), .groups = 'drop') %>%
    dplyr::mutate(category_name_description = category_name_r,
                  channel = str_split_fixed(string = category_name_r, pattern = "_", n = 3)[,1],
                  module = str_split_fixed(string = category_name_r, pattern = "_", n = 3)[,2],
                  compartment = str_split_fixed(string = category_name_r, pattern = "_", n = 3)[,3])

  htpp_category <- mongo(collection="htpp_category", url=mongoUrl, verbose=getOption("verbose"))

  if(rerun == TRUE){
    htpp_category$drop()
  }

  htpp_category$insert(List2, na = "null")

  if(htpp_category$count() != 49){
    warning(paste("Incorrect number of category documents.  Expected 49 categories, instead there are", htpp_category$count(), "categories."))
  }

}


#' Creates htppWell collections (htpp_well_raw, htpp_well, htpp_image_metadata)
#'
#' @param file_path character string: file path to the top level directory of Harmony files for an HTPP dataset (i.e., the directory above plate-level directories)
#' @param mongoUrl character string: The URL of the mongo database holding the collection, with user credentials to access it
#' @param Cell_Type character string or list of strings: the cell type or types being used
#' @param CellArea.Limit dictionary: A dictionary of cells and their corresponding cell area limits, of the form c("celltype" = list(c(lower,upper)))
#' @param NucleiArea.Limit dictionary: A dictionary of cells and their corresponding nuclei area limits, of the form c("celltype" = list(c(lower,upper)))
#' @param SType character string: The sample type used for normalization. Set to "vehicle control" by default
#' @param n_max numeric: The maximum dimensions of the table
#' @param rerun boolean: Whether to drop and replace the collections in the dataframe before inserting
#' @param replace boolean: Whether you want to replace existing records in the mongo database; false by default
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
#' @export generate_htppWell

generate_htppWell <- function(file_path, mongoUrl, Cell_Type, CellArea.Limit, NucleiArea.Limit, SType = "vehicle control", n_max=2000*384, rerun=FALSE, replace=FALSE){
  #Grab all plateIDs:
  htpp_well_trt <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))

  List <- data.table(htpp_well_trt$find())

  List <- unique(List[, c("plate_id", "pg_id", "cell_type")])

  #------------------------------------------------------------------------------------#
  # 2. process every plate (~ 10-12 min/plate)
  #------------------------------------------------------------------------------------#

  #check documents in existing collections
  mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$count()
  mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))$count()
  mongo(collection="htpp_image_metadata", url=mongoUrl, verbose=getOption("verbose"))$count()

  #delete collections is rerun == TRUE
  if(rerun == TRUE){
    mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$drop()
    mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))$drop()
    mongo(collection="htpp_image_metadata", url=mongoUrl, verbose=getOption("verbose"))$drop()
  }

  ##add index
  mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$index(add = '{"plate_id" : 1, "sample_id": 1, "trt_name":1, "pg_id":1, "stype":1, "chem_id":1}')

  file_path<-paste0(file_path, "/")

  for(PlateID in unique(List[, plate_id])){
    #Need to add input file variable for Raw2Level4 function
    FolderList <- list.files(file_path) #note that the number of folder is more than the total number of plates in PlateGroups

    #find the folder that matches the PlateID
    FolderName <- grep(PlateID, FolderList, value=T)

    #Define cell line
    cell_type <- unique(List[plate_id == PlateID, cell_type])
    if(!(cell_type %in% Cell_Type)){
      stop(paste0("You input for cell type, ", Cell_Type, ", does not match the ", cell_type," in htpp_well_trt for plate ",  PlateID))
    }

    if(!(cell_type)%in% names(NucleiArea.Limit)){
      stop(paste("NucleiArea.Limit does not include", cell_type, ", this is required"))
    }

    if(!(cell_type)%in%names(CellArea.Limit)){
      stop(paste("CellArea.Limit does not include", cell_type, ", this is required"))
    }

    if(length(FolderName)==0){
      message("There is no folder for this plate:")
      message(PlateID)
    }else if(file.exists(paste0(file_path, FolderName, "/Evaluation2/Objects_Population - Cells Non-Border.txt"))){  #check if the raw file exist
      InputPath = paste0(file_path, FolderName, "/Evaluation2/")

      if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") & replace == TRUE){
        message(PlateID)
        message(InputPath)
        message("Plate ", PlateID, " already existed and REPLACE = TRUE. Will rerun this plate\n")

        mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))
        mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))
        mongo(collection="htpp_image_metadata", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))

        tic()
        try(Raw2Level4(InputPath = InputPath, PlateID = PlateID, Cell_Type = cell_type, SType = SType, n_max = n_max, mongoUrl = mongoUrl, CellArea.Limit = CellArea.Limit[cell_type], NucleiArea.Limit = NucleiArea.Limit[cell_type]))
        toc()
      }else if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") == FALSE){
        message(PlateID)
        message(InputPath)
        tic()
        try(Raw2Level4(InputPath = InputPath, PlateID = PlateID, Cell_Type = cell_type, SType = SType, n_max = n_max, mongoUrl = mongoUrl, CellArea.Limit = CellArea.Limit[cell_type], NucleiArea.Limit = NucleiArea.Limit[cell_type]))
        toc()
      }else if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") & replace == FALSE){
        message(PlateID)
        message(InputPath)
        message("Plate ", PlateID, " already existed and REPLACE = FALSE. Will skip this plate\n")
        next
      }
    }else if(file.exists(paste0(file_path, FolderName, "/Evaluation1/Objects_Population - Cells Non-Border.txt"))) {
      InputPath = paste0(file_path, FolderName, "/Evaluation1/")

      if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") & replace == TRUE){
        message(PlateID)
        message(InputPath)
        message("Plate ", PlateID, " already existed and REPLACE = TRUE. Will rerun this plate\n")

        mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))
        mongo(collection="htpp_well_raw", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))
        mongo(collection="htpp_image_metadata", url=mongoUrl, verbose=getOption("verbose"))$remove(query=mongoQuery(plate_id = PlateID))

        tic()
        try(Raw2Level4(InputPath = InputPath, PlateID = PlateID, Cell_Type = cell_type, SType = SType, n_max = n_max, mongoUrl = mongoUrl, CellArea.Limit = CellArea.Limit[cell_type], NucleiArea.Limit = NucleiArea.Limit[cell_type]))
        toc()
      }else if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") == FALSE){
        message(PlateID)
        message(InputPath)
        tic()
        try(Raw2Level4(InputPath = InputPath, PlateID = PlateID, Cell_Type = cell_type, SType = SType, n_max = n_max, mongoUrl = mongoUrl, CellArea.Limit = CellArea.Limit[cell_type], NucleiArea.Limit = NucleiArea.Limit[cell_type]))
        toc()
      }else if(PlateID %in% mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$distinct(key = "plate_id") & replace == FALSE){
        message(PlateID)
        message(InputPath)
        message("Plate ", PlateID, " already existed and REPLACE = FALSE. Will skip this plate\n")
        next
      }
    }
  }
}


#' Create http_well_norm, a collection of normalized well data for all plate groups
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing collection and reinsert; have FALSE by default
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
#' @export generate_htppWellNorm
#'
generate_htppWellNorm<-function(mongoUrl, rerun=FALSE){


  #------------------------------------------------------------------------------------#
  # 1. find which plate groups are in the database
  #------------------------------------------------------------------------------------#

  List4 <- mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$find(fields='{"pg_id":1, "sample_id":1, "_id":0}')

  if (length(List4) < 1){
    stop("The htpp_well collection is empty but is required to create htpp_well_norm. Please ensure htpp_well is created before proceeding.")
  }


  #------------------------------------------------------------------------------------#
  # 2. process every plate group (~ 4 min/plate group)
  #------------------------------------------------------------------------------------#


  mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$index(add = '{"plate_id" : 1, "sample_id": 1, "trt_name":1, "pg_id":1, "stype":1, "chem_id":1}')


  #delete collections is rerun == TRUE
  if(rerun == TRUE){
    mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$drop()
  }

  for(PlateGroup in unique(List4$pg_id)){
    message(paste("**************", PlateGroup, "**************"))
    tic()
    Level5(PlateGroup = as.character(PlateGroup), SType = "vehicle control", mongoUrl=mongoUrl)
    toc()
  }

  #check documents in existing collections
  normCount<-mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$count()
  wellCount<-mongo(collection="htpp_well", url=mongoUrl, verbose=getOption("verbose"))$count()

  if(normCount != wellCount){
    warning(paste("Expected", wellCount, "documents in well_norm, based on htpp_well, instead there are", normCount, "documents."))
  }


}

#' Create plots for signal strength and add NULL chemicals to htpp_well_norm
#'
#' @param n_lowest_conc integer: The number of the lowest concentrations in a concentration series for modeling Null chemical data; Default is 2 (i.e., dose_level 1 and 2)
#' @param n_cv_active_conc integer: The number of cell viability active concentrations to be excluded; default is 6 (i.e., exclude chemicals where dose_level >= 6 are cell viability actives)
#' @param rel_cellCount integer: The relative cell count threshold for excluding well data for Null chemical sampling; default is 50 (i.e., exclude wells with rel_cell_count < 50)
#' @param plot_file_path character string: file path where signal strength plots will be created
#' @param study_name character string: the name of the study
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param ConcList numeric vector: vector of 8 test concentrations to be used for the NULL chemicals. c(100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001) by default.
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
#' @export generate_htppNullChems
#'
generate_htppNullChems <- function(n_lowest_conc = 2, n_cv_active_conc = 6, rel_cellCount = 50, plot_file_path, study_name, mongoUrl, ConcList=c(100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001), rerun=FALSE){


  #---------------------------------------------------------------------------------------------#
  # 1. Check to see if null chemicals exist in htpp_well_norm and delete if rerun == TRUE
  #---------------------------------------------------------------------------------------------#

  htpp_well_norm <-  mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))
  if(htpp_well_norm$count() < 1){
    stop("htpp_well_norm collection is needed for well data for this function, but is empty.  Check your mongo database.")
  }

  cv_bmc <-  mongo(collection="cv_bmc", url=mongoUrl, verbose=getOption("verbose"))
  if(cv_bmc$count() < 1){
    stop("cv_bmc collection is needed for biomarker data for this function, but is empty.  Check your mongo database.")
  }

  if(rerun == TRUE){
    htpp_well_norm$remove(query=mongoQuery(stype="null")) #if putting in new null chemicals
  }

  #------------------------------------------------------------------------------------#
  # 2. Identify test samples that were not CV active at x concentration - FOR EACH CELL TYPE
  #------------------------------------------------------------------------------------#

  CVResult_all <- data.table(mongo(collection="cv_bmc", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(stype="test sample")))

  #FOR EACH CELL TYPE
  null_chemicals <- list() #this will store the list of null chems for each cell type
  i <- 1

  for(cell in unique(CVResult_all[, cell_type])){


    CVResult <- CVResult_all[cell_type == cell]
    cat("There are", length(unique(CVResult[, chem_id])), "chemicals in", cell, "cells\n")

    ## only include chemicals that have at least 6 non CV active concentrations
    CVResult <- CVResult %>% filter(cv_noec_dose_level>=6)
    cat("There are", length(unique(CVResult[, chem_id])), "chemicals in", cell, "cells that have at least 6 non CV active concentrations\n")

    #Pull in CPData
    CPData <- data.table(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(stype="test sample", cell_type = cell, dose_level=c(1:n_lowest_conc))))

    #retain only test chemicals that were identified in step 1 as eligible
    CPData <- CPData %>% filter(chem_id %in% CVResult$chem_id)

    #filter out wells with low cell count
    CPData <- CPData %>% filter(rel_cell_count > 50)

    message(paste0(capture.output(table(CPData$pg_id, CPData$replicate_num)), collapse = "\n"))

    ##Generate matrix of Euclidean norm values for each chem_id across all plates
    htpp_signalStrength <- lapply(unique(CPData[, chem_id]), function(x){
      dat <- CPData[chem_id == x, ]
      fCols<-grep(pattern="^f_", colnames(dat), value=TRUE)
      idCols<-colnames(dat[, !c("n_fields", "n_cells_total", "n_cells_keep", "rel_cell_count")])
      idCols<-idCols[!(idCols %in% fCols)]
      dat_long <- melt.data.table(dat[, !c("n_fields", "n_cells_total", "n_cells_keep", "rel_cell_count"), with = FALSE],
                                  id.vars = idCols, value.name = "norm_signal", variable.name = "htpp_feature")
      dat_long$norm_signal<- as.numeric(dat_long$norm_signal, na.rm=TRUE)
      euclidean_norm <- data.table(chem_id = x, norm = Euclidean_norm_vec(dat_long[, norm_signal]))
      return(euclidean_norm)
    })




    htpp_signalStrength <- do.call(rbind, htpp_signalStrength)




    ##What is the outer fence of the signal strength data
    message(paste("the outer fences of the signal strength data for", cell, "are", outerFences(input_vector = htpp_signalStrength[, norm]), ".  "))

    ##Plot distribution of points
    p <- ggplot(data = htpp_signalStrength, aes(x = norm)) +
      geom_boxplot() +
      scale_x_log10() +
      geom_vline(xintercept = outerFences(input_vector = htpp_signalStrength[, norm])[2], linetype = "dashed") +
      geom_vline(xintercept = outerFences(input_vector = htpp_signalStrength[, norm])[1], linetype = "dashed") +
      theme_bw() +
      ggtitle("HTPP Feature Signal Strength", subtitle = "Dashed line is the Tukey's Outer Fence")

    png(filename = paste0(plot_file_path, "/", study_name, "_", cell, "_",  "signalStrengthBoxplot.png"), res = 200, width = 6, height = 3, units = "in")
    print(p)
    dev.off()

    ##Filter chemicals based on Tukey's outer fence of signal strength
    null_chems <- htpp_signalStrength[norm <= outerFences(input_vector = htpp_signalStrength[, norm])[2],]
    null_chems <- null_chems[norm >= outerFences(input_vector = htpp_signalStrength[, norm])[1],]
    null_chems <- data.table(chem_id = null_chems[, chem_id], cell_type = cell)

    cat("For", cell, "cells there are", length(null_chems[, chem_id]), "null chemicals\n")

    null_chemicals[[i]] <- null_chems
    i <- i+1

  } #for each cell type

  #Retrieve all null chemicals into data.table
  null_chemicals <- data.table(do.call(rbind, null_chemicals))

  #------------------------------------------------------------------------------------#
  # 3. Retrieve CP data from db & prepare for null chemical generation - FOR EACH CELL TYPE AND PLATE GROUP
  #------------------------------------------------------------------------------------#

  CPData_all = data.table(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$find(query=mongoQuery(stype="test sample", dose_level=c(1,2))))#11377 obs

  for(cell in unique(CPData_all[, cell_type])){
    for(PG in unique(CPData_all[cell_type == cell, pg_id])){

      CPData <- CPData_all[cell_type == cell & pg_id == PG]

      #retain only test samples that were identified in step 1 as eligible
      CPData <- CPData %>% filter(chem_id %in% null_chemicals[cell_type == cell, chem_id])

      #filter out wells with low cell count
      CPData <- CPData %>% filter(rel_cell_count>50)

      cat("For", cell, "cells", "and plate group", PG, "there should be a total of", floor(min(CPData[, .N, by = c("plate_id", "replicate_num")]$N)/length(ConcList)), "null chemicals\n")

      ## Evaluate how may chemicals are eligible for NULL modeling on each plate group.
      ## Multiply the minimum number of eligible chemicals on any plate group by 2, representing the lowest 2 dose levels. (note this has already been done here)
      ## Divide that number by 8, representing an eight point concentration series.
      ## Example: We have a minimum of 82 chemicals across the 4 plates. 82 / 8 = 10.25. Round down to 10 for the number of null chemicals
      ## This may not be a useful metric for datasets without CV-specific data. See below for alternative calculation.

      null_chems <- floor(min(CPData[, .N, by = c("plate_id", "replicate_num")]$N)/length(ConcList))

      ## delete unneeded meta information
      colnames(CPData)[ which(!grepl("^f_", colnames(CPData))) ]
      CPData = CPData %>% dplyr::mutate(stype=NA, chem_id="", dose_level=NA, trt_name="", conc=NA)

      #------------------------------------------------------------------------------------#
      # 4. Generate null chemicals
      #------------------------------------------------------------------------------------#

      #########  Sample 10 chemicals from each plate
      set.seed(126) #adding seed for reproducibility
      NullSet <- CPData %>%
        group_by(plate_id) %>% dplyr::sample_n(null_chems*length(ConcList), replace=F) %>%
        dplyr::mutate(chem_id = paste0("lowest2conc_1", base::sample(x=rep(letters[1:null_chems], length(ConcList)), replace=F))) %>%
        group_by(plate_id, chem_id) %>%
        dplyr::mutate(conc = base::sample(x=ConcList, replace=F),
                      dose_level = rank(conc),
                      stype = "null",
                      trt_name = paste(pg_id, chem_id, dose_level, sep="_")) %>%
        dplyr::arrange(pg_id, plate_id, chem_id, dose_level) %>% ungroup()

      message(paste0(capture.output(table(NullSet$plate_id, NullSet$conc)), collapse = "\n"))
      #------------------------------------------------------------------------------------#
      # 5. Write data to collection
      #------------------------------------------------------------------------------------#

      tic()
      for(i in 1:dim(NullSet)[1]){
        #print(i)
        newDocument = NullSet[i,]
        htpp_well_norm$insert(newDocument, auto_unbox=TRUE, na = "null")
      }
      toc()

      cat("A total of", dim(NullSet)[1], "null chemical documents were added to htpp_well_norm for", cell, "cells\n\n\n")
    } #for each plate group
  } #for each cell type

  #------------------------------------------------------------------------------------#
  # 6. Check to see that null chemicals were added
  #------------------------------------------------------------------------------------#
  message(paste("There are a total of", htpp_well_norm$count(query = mongoQuery(stype = "null")), "null chemical wells in htpp_well_norm."))
  message(paste("There are", htpp_well_norm$count(query = mongoQuery(stype = c("reference chemical", "vehicle control", "test sample", "viability positive control"))), "wells treated with chemicals that are either reference chemicals, vehicle controls, test samples or viability positive controls in htpp_well_norm."))
  message(paste("There are", htpp_well_norm$count(), "total wells after normalization.  There should be", htpp_well_norm$count(query = mongoQuery(stype = "null"))+htpp_well_norm$count(query = mongoQuery(stype = c("reference chemical", "vehicle control", "test sample", "viability positive control"))), "total wells, based on the sum of the prior two totals."))
}

#' create htpp_profile collection
#'
#' @param n_cells  numeric: Minimum threshold for the number of cells to keep for filtering
#' @param relative_cellCount numeric: Minimum threshold of the count of relative number of cells to use for filtering
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param rerun boolean: rerun = TRUE will drop existing collection and reinsert; FALSE by default
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
#' @export generate_htppProfile
#'
generate_htppProfile <- function(n_cells, relative_cellCount, mongoUrl, rerun=FALSE){
  #------------------------------------------------------------------------------------#
  # 1. Setup
  #------------------------------------------------------------------------------------#




  htpp_profile <- mongo(collection="htpp_profile", url=mongoUrl, verbose=getOption("verbose"))
  htpp_well_norm <- mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))
  htpp_well_trt <- mongo(collection="htpp_well_trt", url=mongoUrl, verbose=getOption("verbose"))



  if(htpp_well_norm$count()<1){
    stop("The htpp_well_norm collection is empty but is required to create htpp_profile. Please ensure htpp_well_norm is created before proceeding.")
  }

  if (htpp_well_trt$count()<1){
    stop("The htpp_well_trt collection is empty but is required to create htpp_profile. Please ensure htpp_well_norm is created before proceeding.")
  }

  message(htpp_profile$count())

  #delete documents in collection if rerun == TRUE
  if(rerun == TRUE){
    htpp_profile$drop()
  }


  #------------------------------------------------------------------------------------#
  # 2. Insert into htpp_profile by cell type
  #------------------------------------------------------------------------------------#
  fieldsOfTrt<-htpp_well_trt$find(fields='{"pg_id":1, "sample_id":1, "_id":0}')
  PGList = as.character(unique(fieldsOfTrt$pg_id))


  for(PG in PGList){

    if (PG == PGList[1]){
      FeatureList <- mongo(collection="htpp_feature", url=mongoUrl)$find() %>% filter(!is.na(feature_id))
      mongo(collection="htpp_profile",  url=mongoUrl, verbose=getOption("verbose"))$index(add = '{  "pg_id":1, "stype":1, "cell_type":1, "chem_id":1, "n_wells":1,
                                                                                           "dose_level":1, "conc":1, "trt_name":1 }')
    }

    message(paste("*****", PG, "******"))
    tic()
    Input <- data.table(mongo(collection="htpp_well_norm", url=mongoUrl)$find(query=mongoQuery(pg_id=PG, stype=c("reference chemical", "test sample"))))
    toc()

    for(cell in unique(Input$cell_type)){
      message(paste("Generating htpp_profile input for", cell, "cells", sep = " "))

      Data <- Input[cell_type == cell, ] %>% dplyr::mutate(use.me = (n_cells_keep>n_cells & rel_cell_count>relative_cellCount)) %>%
        group_by(pg_id, stype, cell_type, chem_id, dose_level, conc, conc_unit, trt_name)

      N <- Data %>% dplyr::summarise(n_wells = sum(use.me, na.rm=T))

      tic()
      Profiles <- Data %>% filter(use.me) %>%
        summarise_at(.vars = FeatureList$feature_name_mongo, .funs = "median", na.rm=T) %>%
        dplyr::mutate_at(.vars = FeatureList$feature_name_mongo, .funs = "round", 3)
      toc()

      Output <- N %>% dplyr::left_join(Profiles)

      htpp_profile <- mongo(collection="htpp_profile", url=mongoUrl, verbose=getOption("verbose"))

      tic()
      #if there are already a documents for this plate group and cell type, delete it
      if(length(htpp_profile$find(query=mongoQuery(pg_id = PG, cell_type = cell), fields='{"_id" : 1}')) > 0){
        cat("Entry for plate group", PG, "and", cell, "cells already existed!!!\n")
        htpp_profile$remove(query=mongoQuery(pg_id = PG, cell_type = cell))
      }

      for(iRow in 1:dim(Output)[1]){
        newDocument = Output[iRow,]
        htpp_profile$insert(newDocument, auto_unbox=TRUE, na = "null")
      }#for each well
      toc()

      cat("Inserted", dim(Output)[1], "documents into htpp_profile for", cell, "cells\n")
    } #for each cell type

    if(PG=="1"){
      mongo(collection="htpp_profile",  url=mongoUrl, verbose=getOption("verbose"))$index(add = '{  "pg_id":1, "stype":1, "cell_type":1, "chem_id":1, "n_wells":1,
                                                                                           "dose_level":1, "conc":1, "trt_name":1 }')
    }

    rm(Input, N, Profiles, Output, newDocument)
  }


  return(paste("The htpp_profile collection now contains", htpp_profile$count(), "documents after insert."))
}

#' Calculate, plot and record global Mahalanobis distances from mongo data
#'
#' @param coverVariance numeric: The value of variance explained used to determine the number of eigen features used in analysis
#' @param minObjects numeric: The minimum number of objects used to filter the dataset for analysis
#' @param plot_file_path character string: file path where variance explained plots will be created
#' @param study_name character string: the name of the experiment used to title the plots
#' @param mongoUrl URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
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
#' @import ggplot2
#'
#'
#' @export generate_htppGlobalMah
#'
generate_htppGlobalMah <- function(coverVariance, minObjects, plot_file_path, study_name, mongoUrl, rerun=FALSE){
  #---------------------------------------------------------------------------------------------#
  # 1. Connect and check all collections and delete if rerun == TRUE
  #---------------------------------------------------------------------------------------------#
  htpp_well_norm<-mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))

  if(htpp_well_norm$count()<1){
    stop("The htpp_well_norm collection is empty but is required to create htpp_global_mah. Please ensure htpp_well_norm is created before proceeding.")
  } else if(htpp_well_norm$count(query=mongoQuery(stype="null"))<1){
    stop("There are no NULL chemicals in the htpp_well_norm collection.  Please check that generate_htppNullChems() ran correctly")
  }



  #global mahalanobis collection
  htpp_global_mah <- mongo(collection="htpp_global_mah", url=mongoUrl, verbose=getOption("verbose"))


  if(rerun == TRUE){
    htpp_global_mah$drop()
  }

  #---------------------------------------------------------------------------------------------#
  # 2. Global Mahalanobis Distance Calculation - FOR EACH CELL TYPE
  #---------------------------------------------------------------------------------------------#

  Table1_all <- data.table(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$find(fields='{"_id":0}'))

  for(cell in unique(Table1_all[, cell_type])){

    Table1 <- as.data.frame(Table1_all[cell_type == cell])

    tic()
    Output <- globalMahalanobisDistances(Table1 = Table1, coverVariance = coverVariance, minObjects = minObjects, SType = "vehicle control",  url = mongoUrl)
    toc()

    CumProportion <- Output$CumProportion

    #How many PC are needed to explain  at least x% of variance?
    PC90 <- length(which(CumProportion<0.90))+1
    PC95 <- length(which(CumProportion<0.95))+1
    PC99 <- length(which(CumProportion<0.99))+1

    png(paste0(plot_file_path, "/", study_name, "_", cell,"_", "htpp_global_mah_Proportion_of_variance.png"), width=8, height=6, units="in", res=144)
    plot(x=1:1300, y=CumProportion, col="gray50", pch=19, cex=0.5, type="p",
         ylim=c(0,1), xlab="# of components", ylab="Proportion of variance retained", main=paste(study_name, "Principal components of telohaec_apcra data for",cell, "cells", sep = " "))
    #horizontal part
    segments(x0=30, y0=0.90, x1 = PC90, col="blue", lty='dashed')
    segments(x0=30, y0=0.95, x1 = PC95, col="blue", lty='solid', lwd=2)
    segments(x0=30, y0=0.99, x1 = PC99, col="blue", lty='dotted')
    #vertical part
    segments(x0=PC90, y0=0.1, y1 = 0.90, col="blue", lty='dashed')
    segments(x0=PC95, y0=0.1, y1 = 0.95, col="blue", lty='solid', lwd=2)
    segments(x0=PC99, y0=0.1, y1 = 0.99, col="blue", lty='dotted')

    text(x=c(PC90, PC95, PC99), y=0.05, labels=c(PC90, PC95, PC99), srt=90)
    text(x=0, y=c(0.9, 0.95, 0.99),  labels=paste0(c(90,95,99), "%"), cex=0.7)
    dev.off()

    #Check if all data is there
    normCount<-mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$count(query=mongoQuery(cell_type=cell))
    gMahCount<-htpp_global_mah$count(query=mongoQuery(cell_type=cell))
    if(gMahCount != normCount){
      warning(paste("Expected", normCount, "documents in htpp_global_mah, based on htpp_well_norm, for", cell, ", instead there are", gMahCount, "documents."))
    } else {
      message(paste("htpp_global_mah collection was created successfully, contains", gMahCount, "records for", cell, " This matches the size of the htpp_well_norm collection for this cell type"))
    }
  } #for all cell types
}

#' Create htpp collection htpp_cat_mah for category Mahalanobis distances
#'
#' @param coverVariance numeric: The value of variance explained used to determine the number of eigen features used in analysis
#' @param minObjects numeric: The minimum number of objects used to filter the dataset for analysis
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param varianceExplainedPath character string: the path where the function will write variance explained metadata
#' @param nThreads numeric:  the number of threads to use for processing; default is 1
#' @param rerun boolean: rerun = TRUE will drop existing htpp_cat_mah collection and reinsert; have FALSE by default
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
#' @import parallel
#' @import doParallel
#'
#'
#' @export
#'
generate_htppCatMah <- function(coverVariance, minObjects, mongoUrl, varianceExplainedPath, nThreads=1, rerun=FALSE){
  #---------------------------------------------------------------------------------------------#
  # 1. Connect and check collection and delete if rerun == TRUE
  #---------------------------------------------------------------------------------------------#

  htpp_well_norm<-mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))
  if(htpp_well_norm$count()<1){
    stop("The htpp_well_norm collection is empty but is required to create htpp_cat_mah. Please ensure htpp_well_norm is created before proceeding.")
  } else if(htpp_well_norm$count(query = mongoQuery(stype = "null"))<1){
    stop("The htpp_well_norm collection has no records with 'stype':'null' which are required for category mahalanobis.  Please check that htpp_well_norm has the correct recordsS")
  }

  #category mahalanobis collection
  htpp_cat_mah <- mongo(collection="htpp_cat_mah", url=mongoUrl, verbose=getOption("verbose"))

  if(rerun == TRUE){
    htpp_cat_mah$drop()
  }

  #add index
  mongo(collection = "htpp_cat_mah", url = mongoUrl, verbose = getOption("verbose"))$index(add = '{"pg_id" : 1, "stype": 1, "chem_id": 1, "category_name_r ": 1}')

  #---------------------------------------------------------------------------------------------#
  # 2. Category-level Mahalanobis Distance Calculation - FOR EACH CELL TYPE
  #---------------------------------------------------------------------------------------------#

  #load all well data
  Table1_all <- data.table(mongo(collection="htpp_well_norm", url=mongoUrl, verbose=getOption("verbose"))$find(fields='{"_id":0}'))


  #grab category information
  CategoryList <- mongo(collection  ="htpp_feature", url = mongoUrl, verbose = getOption("verbose"))$find() %>%
    select(category_name_r) %>% distinct() %>% filter(!is.na(category_name_r))


  FeatureList <- mongo(collection="htpp_feature", url=mongoUrl, verbose=getOption("verbose"))$find() %>% as_tibble()


  #CategoryList = CategoryList %>% mutate(n=row_number()) %>% filter(n<=40& n>20) %>% select(-n)
  #CategoryList = CategoryList %>% mutate(n=row_number()) %>% filter(n>40) %>% select(-n)

  my.cluster <- parallel::makeCluster(
    nThreads,
    type = "PSOCK"
  )

  registerDoParallel(cl <- my.cluster)

  clusterExport(cl, c("categoryMahalanobisDistances", 'mongoQuery'))

  CatNames<-CategoryList$category_name_r
  CatSplit<-split(CatNames, cut(seq_along(CatNames),3,labels = FALSE))


  #define combine function for the foreach loops
  comb <- function(x, ...){
    lapply(seq_along(x),FUN = function(i){ c(x[[i]], lapply(list(...), function(y) y[[i]]))})
  }

  for(cell in unique(Table1_all[, cell_type])){

    Table1 <- as.data.frame(Table1_all[cell_type == cell])
    varianceExplained1<-c()
    varianceExplained2<-c()
    varianceExplained3<-c()



    tic()
    ResultList<-foreach(x=CatSplit$`1`, .packages=c('tidyr', 'dplyr', 'stringr', 'mongolite'), .combine = comb, .init=list(list(), list())) %dopar% {
      Output <- try(categoryMahalanobisDistances(Level5=Table1, FeatureList = FeatureList, CategoryName = x,
                                                 coverVariance = coverVariance, minObjects = minObjects, SType = "vehicle control",  mongoUrl = mongoUrl))

      list(Output[[1]], Output[[2]])
    }

    newLine <- do.call(rbind, ResultList[[1]])
    varianceExplained1 <- do.call(rbind, ResultList[[2]])

    #bulk insert into collection
    htpp_cat_mah$insert(newLine, auto_unbox=TRUE, na = "null")

    ResultList<-foreach(x=CatSplit$`2`, .packages=c('tidyr', 'dplyr', 'stringr', 'mongolite'), .combine = comb, .init=list(list(), list())) %dopar% {
      Output <- try(categoryMahalanobisDistances(Level5=Table1, FeatureList = FeatureList, CategoryName = x,
                                                 coverVariance = coverVariance, minObjects = minObjects, SType = "vehicle control",  mongoUrl = mongoUrl))

      list(Output[[1]], Output[[2]])
    }

    newLine <- do.call(rbind, ResultList[[1]])
    varianceExplained2 <- do.call(rbind, ResultList[[2]])

    #bulk insert into collection
    htpp_cat_mah$insert(newLine, auto_unbox=TRUE, na = "null")

    ResultList<-foreach(x=CatSplit$`3`, .packages=c('tidyr', 'dplyr', 'stringr', 'mongolite'), .combine = comb, .init=list(list(), list())) %dopar% {
      Output <- try(categoryMahalanobisDistances(Level5=Table1, FeatureList = FeatureList, CategoryName = x,
                                                 coverVariance = coverVariance, minObjects = minObjects, SType = "vehicle control",  mongoUrl = mongoUrl))

      list(Output[[1]], Output[[2]])
    }

    newLine <- do.call(rbind, ResultList[[1]])
    varianceExplained3 <- do.call(rbind, ResultList[[2]])

    #bulk insert into collection
    htpp_cat_mah$insert(newLine, auto_unbox=TRUE, na = "null")

    #combine all varianceExplained outputs
    varianceExplained <- rbind(varianceExplained1, varianceExplained2, varianceExplained3)
    write_csv(varianceExplained, paste0(varianceExplainedPath, "/", cell, "_varianceExplained.csv") )

    toc()


  }
  stopCluster(cl)

  #Check if all data is there
  if(htpp_cat_mah$count()/49 != htpp_well_norm$count()){
    warning(paste("Expected", (49*htpp_well_norm$count()), "documents in htpp_cat_mah, based on 49 * htpp_well_norm, instead there are", htpp_cat_mah$count(), "documents."))
  }
}


#' Create htpp_bmc collection based on htpp_tcpl and adds global mahalanobis distances into htpp_bmc
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param hitCall numeric (between 0-1): Hitcall threshold from tcplfit2 to use for filtering good BMD values; default is 0 for no hitcall filtering
#' @param bmc_max numeric: The maximum bmc value if bmd > highest tested conc; default is NA
#' @param bmc_min Defines the denominator for calculating the minimum bmc value for cases where the bmc is less that the lowest tested conc (i.e., minimum tested conc/bmc_min); default is 10^0.5
#' @param rerun rerun = TRUE will drop existing entries in htpp_bmc for approach = "global" and endpoint = "global", and reinsert; FALSE by default
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
#' @export generate_htppBmc_globalMah
#'
generate_htppBmc_globalMah <- function(mongoUrl, hitCall=0, bmc_max=NA, bmc_min=10^0.5, rerun=FALSE){
  htpp_tcpl <- mongo(collection = "htpp_tcpl", url = mongoUrl, verbose = getOption("verbose"))
  if (htpp_tcpl$count(query = mongoQuery(approach = "global", endpoint = "global"))<1){
    stop("The there are no global documents in htpp_tcpl. These are needed to generate htpp_bmc.  Please ensure htpp_tcpl is filled correctly before proceeding.")
  }

  ## 1) Write into htpp_bmc
  htpp_bmc <- mongo(collection = "htpp_bmc", url=mongoUrl, verbose=getOption("verbose"))
  htpp_bmc$count(query = mongoQuery(approach = "global", endpoint = "global"))

  #if rerun == TRUE remove global fits from collection
  if(rerun == TRUE){
    htpp_bmc$remove(query = mongoQuery(approach = "global", endpoint = "global"))
  }

  htpp_bmc$index(add = '{"pg_id" : 1, "stype": 1, "cell_type": 1, "chem_id": 1, "approach": 1, "endpoint": 1}')

  CPData <- data.table(htpp_tcpl$find(query = mongoQuery(approach = "global", endpoint = "global")))


  CPData <- CPData %>%
    dplyr::mutate(bmc = ifelse(hitcall > hitCall, bmd, NA),
                  bmc = ifelse(bmc > max_conc, bmc_max, bmc),
                  bmc = ifelse(bmc < min_conc/bmc_min, min_conc/bmc_min, bmc),
                  bmc = signif(bmc, 3),
                  cp_flag = ifelse(!is.na(bmc) & bmc < min_conc, T, F),
                  cp_flag = ifelse(n_conc < 4, NA, cp_flag))

  message(paste0(capture.output(table(CPData$stype, !is.na(CPData$bmc))), collapse = "\n")) #look across all cell types



  for(i in 1:dim(CPData)[1]){
    htpp_bmc$insert(CPData[i,], auto_unbox=TRUE, na = "null")
  }
  message(paste("Inserted", dim(CPData)[1], "bmc documents to htpp_bmc for all chemicals and cell types", sep = " "))

  #check collection dimension after inserting documents
  bmcCount<-htpp_bmc$count(query = mongoQuery(approach = "global", endpoint = "global"))
  tcplGlobCount<-htpp_tcpl$count(query = mongoQuery(approach = "global", endpoint = "global"))
  if(bmcCount != tcplGlobCount){
    warning(paste("Expected", tcplGlobCount, "documents in htpp_bmc based on htpp_tcpl.  Instead, there are", bmcCount, "documents."))
  }
}


#' Print out Category-level Null Chemical Maximum Hitcall Probabilities
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param null_max_conc integer: Maximum concentration of Null chemicals; default is 100 (uM)
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
#' @export nullProbs_catMah
#'
nullProbs_catMah <- function(mongoUrl, null_max_conc = 100){

  ## 1.a) load data
  htpp_tcpl <- mongo(collection = "htpp_tcpl", url = mongoUrl, verbose = getOption("verbose"))
  if(htpp_tcpl$count(query=mongoQuery(approach = "category")) <1 ){
    stop("The htpp_tcpl collection has no category documents.  These are required to fill in htpp_bmc. Please ensure htpp_tcpl is created correctly before proceeding.")
  }

  CPData <- data.table(htpp_tcpl$find(query = mongoQuery(approach = "category")))

  # Investigate the null chemicals
  NullChem <- CPData %>% filter(stype == "null") %>%
    dplyr::mutate(Hitcall = ifelse(bmd > null_max_conc | is.na(bmd), 0, hitcall))

  ZF_Null <- NullChem %>% dplyr::group_by(pg_id, cell_type, chem_id) %>% dplyr::summarise(maxHitcall = max(Hitcall))
  hist(ZF_Null$maxHitcall)

  #print out hitcall probabilities
  message(paste0(capture.output(quantile(ZF_Null$maxHitcall, probs=c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99))), collapse = "\n"))

}

#' Add category Mahalanobis distance information to htpp_bmc collection
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param hitCall numeric (between 0-1): Hitcall threshold from tcplfit2 to use for filtering good BMD values; default is 0.95
#' @param bmc_max numeric: The maximum bmc value if bmd > highest tested conc; default is NA
#' @param bmc_min numeric: Defines the denominator for calculating the minimum bmc value for cases where the bmc is less that the lowest tested conc (i.e., minimum tested conc/bmc_min); default is 10^0.5
#' @param rerun boolean: rerun = TRUE will drop existing entries in htpp_bmc for approach = "global" and endpoint = "global", and reinsert; FALSE by default
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
#' @export
#'
generate_htppBmc_catMah <- function(mongoUrl, hitCall=0.95, bmc_max=NA, bmc_min=10^0.5, rerun=FALSE){
  #---------------------------------------------------------------------------------------------#
  # 1.Apply the thresholds & move data from htpp_tcpl --> htpp_bmc
  #---------------------------------------------------------------------------------------------#

  ## 1.a) load data
  htpp_tcpl <- mongo(collection = "htpp_tcpl", url = mongoUrl, verbose = getOption("verbose"))
  if(htpp_tcpl$count(query=mongoQuery(approach = "category")) <1 ){
    stop("The htpp_tcpl collection has no category documents.  These are required to fill in htpp_bmc. Please ensure htpp_tcpl is created correctly before proceeding.")
  }

  htpp_bmc <- mongo(collection = "htpp_bmc", url = mongoUrl, verbose = getOption("verbose"))

  #if rerun == TRUE delete global documents in htpp_pac
  if(rerun == TRUE){
    htpp_bmc$remove(query = mongoQuery(approach = "category"))
  }

  # index
  htpp_bmc$index(add = '{"pg_id" : 1, "stype": 1, "cell_type": 1, "chem_id": 1, "approach": 1, "endpoint": 1}')


  CPData <- data.table(htpp_tcpl$find(query = mongoQuery(approach = "category")))


  # Apply `thresholds`; identify valid BMCs
  CPData <- CPData %>%
    dplyr::mutate(bmc = ifelse(hitcall > hitCall, bmd, NA),
                  bmc = ifelse(bmc > max_conc, bmc_max, bmc),
                  bmc = ifelse(bmc < min_conc/bmc_min, min_conc/bmc_min, bmc),
                  bmc = signif(bmc, 3),
                  cp_flag = ifelse(!is.na(bmc) & bmc < min_conc, T, F),
                  cp_flag = ifelse(n_conc < 4, NA, cp_flag))


  # Write into htpp_bmc


  for(i in 1:dim(CPData)[1]){
    htpp_bmc$insert(CPData[i,], auto_unbox = TRUE, na = "null")
  }

  #count the number of inserted documents
  if (htpp_bmc$count(query = mongoQuery(approach = "category")) != htpp_tcpl$count(query = mongoQuery(approach = "category"))){
    warning(paste("Expected", htpp_tcpl$count(query = mongoQuery(approach = "category")), "documents added to htpp_bmc, based on category information from htpp_tcpl, instead there are", htpp_bmc$count(query = mongoQuery(approach = "category"))))
  }



}

#' generate htpp_pac from htpp_bmc and add global mahalanobis distance records to htpp_pac
#'
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param hit_n_conc numeric: Number of test concentrations needed during curve fitting to determine if a PAC is a hit; default is 4
#' @param rerun Boolean: TRUE will drop existing entries in htpp_bmc for approach = "global" and endpoint = "global", and reinsert; default is FALSE
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
#' @export generate_htppPac_globalMah
#'
generate_htppPac_globalMah <- function(mongoUrl, hit_n_conc=4, rerun=FALSE){

  htpp_bmc <- mongo(collection = "htpp_bmc", url = mongoUrl, verbose = getOption("verbose"))
  CPData <- data.table(htpp_bmc$find(query = mongoQuery(approach = "global", endpoint = "global")))
  if (length(CPData) <1){
    stop("htpp_bmc collection has no records with approach:'global' and endpoint:'global'.  These are needed to generate htpp_pac.  Please check that htpp_bmc was generated correctly.")
  }
  PAC <- CPData %>% select(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc,
                           approach, top_over_cutoff, hitcall, bmc, cp_flag) %>%
    dplyr::rename(pac = bmc) %>%
    dplyr::mutate(hit = ifelse(n_conc < hit_n_conc, NA, !is.na(pac)))



  Trt <- mongo(collection = "htpp_well_norm", url = mongoUrl, verbose = getOption("verbose"))$find(fields = '{ "_id":0, "pg_id":1, "stype" :1, "cell_type":1, "chem_id":1, "dose_level":1, "conc":1 }')

  Trt <- Trt %>% filter(stype != "vehicle control")  %>% distinct() %>%
    dplyr::arrange(pg_id, chem_id, dose_level)


  Data <- Trt %>% left_join(PAC) %>% filter(conc > pac | (cp_flag & dose_level == 1))

  LOEC <- Data %>% group_by(pg_id, stype, cell_type, chem_id, cp_flag) %>%
    dplyr::summarise(cp_loec_dose_level = min(dose_level)) %>% ungroup() %>%
    dplyr::mutate(cp_loec_dose_level = ifelse(cp_flag, 0, cp_loec_dose_level))

  PAC <- PAC %>% left_join(LOEC) %>%
    dplyr::mutate_at(.vars = c("min_conc", "max_conc", "pac"), .funs = "log10") %>%
    dplyr::mutate_at(.vars = c("min_conc", "max_conc", "top_over_cutoff", "hitcall", "pac"), .funs = "round", 3)

  message(paste0(capture.output(table(PAC$stype, PAC$cp_loec_dose_level, PAC$cell_type, useNA = "ifany")), collapse = "\n"))

  htpp_pac <- mongo(collection = "htpp_pac", url = mongoUrl, verbose = getOption("verbose"))
  htpp_pac$count(query = mongoQuery(approach = "global"))


  if(rerun == TRUE){
    htpp_pac$remove(query = mongoQuery(approach = "global"))
  }

  htpp_pac$index(add = '{"pg_id" : 1, "stype": 1, "cell_type":1, "chem_id": 1, "approach": 1}')


  for(i in 1:dim(PAC)[1]){
    htpp_pac$insert(PAC[i,], auto_unbox=TRUE, na = "null")
  }

  message(paste0(capture.output(table(PAC$stype, PAC$hit, useNA = "ifany")), collapse = "\n"))

  bmcCount<-htpp_bmc$count(query = mongoQuery(approach = "global", endpoint = "global"))
  pacCount<-htpp_pac$count(query = mongoQuery(approach = "global"))
  if(bmcCount != pacCount){
    warning(paste("Expected", bmcCount, "documents in htpp_pac based on htpp_bmc.  Instead, there are", pacCount, "documents."))
  }


}


#' add category mahalanobis records to htpp_pac
#' @param mongoUrl character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl
#' @param hit_n_conc numeric: Number of test concentrations needed during curve fitting to determine if a PAC is a hit; default is 4
#' @param rerun boolean: rerun = TRUE will drop existing entries in htpp_bmc for approach = "category", and reinsert; FALSE by default
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
#' @export generate_htppPac_catMah
#'
generate_htppPac_catMah <- function(mongoUrl, hit_n_conc=4,rerun=FALSE){

  htpp_bmc<-mongo(collection = "htpp_bmc", url = mongoUrl, verbose = getOption("verbose"))
  if(htpp_bmc$count(query = mongoQuery(approach = "category"))<1){
    stop("The htpp_bmc collection needs to first be created and contain documents for approach = 'category'. Please ensure those records were input correctly before proceeding.")
  }

  CPData<-htpp_bmc$find(query = mongoQuery(approach = "category"))


  PAC <- CPData %>% dplyr::group_by(pg_id, stype, cell_type, chem_id, min_conc, max_conc, n_conc, approach) %>%
    dplyr::summarise(top_over_cutoff = maxJN(top_over_cutoff),
                     hitcall = maxJN(hitcall),
                     pac = minJN(bmc),
                     cp_flag = any(cp_flag),
                     n_cat = sum(!is.na(bmc))) %>% ungroup() %>%
    dplyr::mutate(hit = ifelse(n_conc < 4, NA, !is.na(pac)),
                  n_cat = ifelse(n_conc < 4, NA, as.integer(n_cat)),
                  hitcall = ifelse(n_conc < 4, NA, hitcall))


  # Get the treatment info to define the LOEC

  Trt <- mongo(collection = "htpp_well_norm", url = mongoUrl, verbose = getOption("verbose"))$find(fields = '{ "_id":0, "pg_id":1, "stype" :1, "cell_type": 1, "chem_id":1, "dose_level":1, "conc":1 }')

  Trt <- Trt %>% dplyr::filter(stype != "vehicle control")  %>% dplyr::distinct() %>%
    dplyr::arrange(pg_id, cell_type, chem_id, dose_level)

  # retain only dose-levels that are above the PAC

  Data <- Trt %>% left_join(PAC) %>% filter(conc > pac | (cp_flag & dose_level == 1))

  LOEC <- Data %>% group_by(pg_id, stype, cell_type, chem_id, cp_flag) %>%
    dplyr::summarise(cp_loec_dose_level = min(dose_level)) %>% ungroup() %>%
    dplyr::mutate(cp_loec_dose_level = ifelse(cp_flag, 0, cp_loec_dose_level))

  PAC <- PAC %>% left_join(LOEC) %>%
    dplyr::mutate_at(.vars = c("min_conc", "max_conc", "pac"), .funs="log10") %>%
    dplyr::mutate_at(.vars = c("min_conc", "max_conc", "top_over_cutoff", "hitcall", "pac"), .funs = "round", 3)

  message(paste0(capture.output(table(PAC$stype, PAC$cp_loec_dose_level, PAC$cell_type, useNA = "ifany")), collapse = "\n"))


  ## 3.b) Write into htpp_pac
  htpp_pac <- mongo(collection = "htpp_pac", url = mongoUrl, verbose = getOption("verbose"))

  #if rerun == TRUE delete category documents in htpp_pac
  if(rerun == TRUE){
    htpp_pac$remove(query = mongoQuery(approach = "category"))
  }

  #insert documents into htpp_pac
  for(i in 1:dim(PAC)[1]){
    htpp_pac$insert(PAC[i,], auto_unbox = TRUE, na = "null")
  }

  if (htpp_pac$count(query = mongoQuery(approach = "category")) != htpp_bmc$count(query = mongoQuery(approach = "category"))/49){
    warning(paste("Expected", htpp_bmc$count(query = mongoQuery(approach = "category"))/49, "documents added to htpp_pac based on htpp_bmc.  Instead, there are", htpp_pac$count(query = mongoQuery(approach = "category")), "documents."))
  }
}

