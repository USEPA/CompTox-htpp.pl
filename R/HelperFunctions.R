#' Convert 2 digit values to 3 digit in a vector by adding a preceding 0
#'
#' @param Vector The vector being modified
#'
#' @return The vector with all 3 digit values
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
#' @export Well3Digit
#'
#'
#' @examples Well3Digit(c("100","10","30"))
Well3Digit<-function(Vector){
  i=which( str_length(Vector)==2 )
  if(length(i)>0){
    Vector[i]=paste(str_sub(Vector[i], 1, 1),
                    "0",
                    str_sub(Vector[i], 2, 2), sep="" )
  }
  return(Vector)
}

#' calculate the minimum, but return NA is no number is present (pmin doesn't work with summarise)
#'
#' @param x a vector one wants the minimum of
#'
#' @return The minimum, or NA if there are no nonzero numbers in X
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
#' @export minJN
#'
#' @examples
#' minJN(c(4,1,7,10,9))
#' minJN(c(NA, NA, NA, NA, NA, NA, NA, NA, 0))
minJN <- function(x){
  if(sum(!is.na(x))==0){
    y = NA
  }else{
    y=min(x, na.rm=T)
  }
  return(y)
}

#' calculate the maximum, but return NA is no number is present (pmax doesn't work with summarise)
#'
#' @param x a vector one wants the maximum of
#'
#' @return The maximum, or NA if there are no nonzero numbers in X
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
#' @export maxJN
#'
#' @examples
#' maxJN(c(4,1,7,10,9))
#' maxJN(c(NA, NA, NA, NA, NA, NA, NA, NA, 0))
maxJN <- function(x){
  if(sum(!is.na(x))==0){
    y = NA
  }else{
    y=max(x, na.rm=T)
  }
  return(y)
}

######################### by Logan Everett (updated for use with the SCDCD server) ##################################
# Just format a mongo URL based on mongoServer, username, password, and database
# This currently does no validation to check for problems
# It also does not currently support anonymous access (no username/password)


#' Combines credentials into a MongoURL to access a database
#'
#' @param host The database host
#' @param user The MongoDB username
#' @param passwd The password for that database
#' @param db The name of the database
#' @param authSource The authentication source
#' @param authMechanism The authentication mechanism
#'
#' @return a MongoDB URL granting read access to the database
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
#' @export mongoURL
#'
#' @examples
#' mongoURL(host="test", user="readonly", passwd="passwd", db="test_db")
mongoURL <- function(host, user, passwd, db, authSource="admin", authMechanism="SCRAM-SHA-256") {
  passwd <- URLencode(passwd, reserved = TRUE)  # Protect any special characters in passwd
  if(is.null(authSource) & is.null(authMechanism)){
    paste0("mongodb://", user, ":", passwd, "@", host, "/", db)
  }else {
    paste0("mongodb://", user, ":", passwd, "@", host, "/", db, "?", "authSource=", authSource, "&", "authMechanism=", authMechanism)
  }
}

#' Build a mongo query from an arbitrary set of params or a list
# Takes a list of named values/vectors and converts to a JSON-format mongo query. Each entry with length > 1 is interpreted
# as a set of potential matching values and converted to the $in keyword for mongo
#'
#' @param ... (any) = Any set of named parameters OR a single named list of all arguments
#'
#' @return JSON query ready to pass to mongolite functions
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
#' @export mongoQuery
#'
#' @examples
#' mongoQuery(approach = "global", endpoint = "global")
mongoQuery <- function(...) {
  # Get all args as a list:
  args <- list(...)
  # If singular list argument, elevate that up to args
  if((length(args)==1) && (class(args) == "list") && is.null(names(args))) {
    args <- args[[1]]
  }
  # Drop members that are NULL or have length 0
  for(key in names(args)) {
    if(is.null(args[[key]]) || (length(args[[key]])==0)) {
      args[[key]] <- NULL
    }
  }
  # If args is length 0, return '{}'
  if(length(args)==0) {
    null_query <- '{}'
    class(null_query) <- "json"
    return(null_query)
  }
  # Loop through args and convert any vector args to $in construct
  for(key in names(args)) {
    if((length(args[[key]]) > 1) && (class(args[[key]]) != "list")) {
      args[[key]] <- list('$in'=args[[key]])
    }
  }
  # Convert to JSON and return
  return(toJSON(args, auto_unbox = T))
}

#' Calculate the Euclidean norm of a vector
#'
#' @param vect Numeric vector: the vector you will take the euclidean norm of
#'
#' @return Numeric: the Euclidean norm of the vector
#' @export Euclidean_norm_vec
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
#' @examples Euclidean_norm_vec(c(1,5,3,4,12))
Euclidean_norm_vec <- function(vect){
  sqrt(sum(vect^2))
}

#' Tukey Outer fence function
#'
#' @param input_vector vector: the data whose outer fences you want to calculate
#' @param iqr.factor numeric: the interquartile range factor.  3 by default.
#'
#' @return upper and lower fences
#' @export outerFences
#'
#' @examples
#' outerFences(input_vector= c(113.86844, 108.47126,
#' 125.22345, 115.17092, 123.61978, 112.45098, 128.88594,
#' 100.98654, 103.95449, 109.41060, 114.59485, 107.66969,
#' 101.75302, 100.46038, 103.86191, 110.25791,
#' 113.22980, 101.13067, 112.93010, 105.88312, 98.46702,
#' 92.14626, 97.33307, 117.61402, 110.81441,
#' 106.10610, 110.93701, 115.87897, 116.07416, 104.51375,
#' 106.85793, 117.94644, 111.48693, 122.34972,
#' 103.06462, 127.57024, 120.70566, 98.62746, 110.22989,
#' 150.98643), iqr.factor=3)
#'
outerFences <- function(input_vector, iqr.factor = 3) {
  input_vector <- as.numeric(input_vector, na.rm=TRUE)
  x_iqr <- IQR(input_vector, na.rm=TRUE)
  x_q1 <- stats::quantile(x=input_vector, probs=0.25, na.rm=TRUE)
  x_q3 <- stats::quantile(x=input_vector, probs=0.75, na.rm=TRUE)
  x_lower <- x_q1 - (iqr.factor * x_iqr)
  x_upper <- x_q3 + (iqr.factor * x_iqr)
  return(c(x_lower, x_upper))
}
