# Proteome Coverage Functions
#
# Proteome coverage measures the non-redundant length of proteoform sequences
# In measuring the proteome coverage, the proteoform sequences need to be merged
# into contigs. Then the total length of the contigs is calculated.
#

# Total length of sequences in human Swiss-Prot as of 6/20/2021
human_proteome_aa <- 11361310


#' Detect vector overlap
#'
#' @param x Vector of start and end indices for first sequence
#' @param y Vector of start and end indices for second sequence
#'
#' @return Merged vector if there is an overlap or 0 if no overlap
#'
#' @examples
#' is_overlap(c(1,10), c(8, 20))
is_overlap <- function(x, y) {
  # Checks if there is an overlap in two numeric vectors. If there is, returns
  # contiguous vector
  if(x[1] >= y[1] & x[1] <= y[2]){
    # Start index is between start and end of contig
    
    if(x[2] >= y[1] & x[2] <= y[2]){
      # End index is between start and end of contig
      # In this case, the proteoform is fully encapsulated by contig. No change is needed
      return(y)
    } else {
      # End index is not between start and end of contig
      # There is an overlap but the end of the proteoform overhangs. End index of contig should be updated
      return(c(y[1], x[2]))
    }
  } else {
    # Start index is not between start and end of contig
    if(x[2] >= y[1] & x[2] <= y[2]){
      # End index is between start and end of contig
      # In this case, there is an overlap but the start of the proteoform overhangs. Start index of contig should be updated
      return(c(x[1], y[2]))
    } else {
      # There is no overlap, so return false
      return(0)
    }
  }
}

#' Calculate proteome coverage contributed by a single protein
#'
#' @param df Collection of proteoform start and end indices
#'
#' @return Sum of proteoform contig lengths
proteome_coverage <- function(df) {
  # We don't want to waste time merging things that are identical
  indicies <- select(df, Start_Index, End_Index) %>% unique()
  
  # Merge the proteoforms into contigs
  contigs <- assemble_contigs(indicies)
  
  # sometimes the first round misses some merges. need to figure that out
  # in the meantime, running assembly again seems to get everything
  contigs <- assemble_contigs(contigs)
  
  # Count the length of contigs
  return(sum_contigs(contigs))
}


#' Proteoform contig assembler
#' 
#' Initiates assembly of a collection of proteoforms into contiguous sequences.
#' Uses the first proteoform as a query and calls \code{make_contigs}.
#'
#' @param proteoforms Collection of proteoform indices to merge
#'
#' @return DataFrame of non-overlapping proteoform contigs
#'
#' @examples
assemble_contigs <- function(proteoforms) {
  new_query <- c(proteoforms$Start_Index[1], proteoforms$End_Index[1])
  return(make_contigs(new_query, proteoforms[-1,]))
}


#' Proteoform contig merging engine
#' 
#' Rather than performing sequence alignments, proteoform contigs are assembled
#' based on their position in the parent protein sequence.
#'
#' @param query Vector of query proteoform indices
#' @param proteoforms Dataframe of proteoform indices to compare the query to.
#' Must contain two columns labelled \code{Start_Index} and \code{End_Index}.
#'
#' @return Dataframe of proteoform indices. If the \code{query} overlaps with
#' any of the entries in \code{proteoforms} then a merge is performed. If no
#' overlap is detected then the \code{query} is appended to the end of
#' \code{proteoforms}
make_contigs <- function(query, proteoforms) {
  # if proteoforms is empty, return the query
  if(nrow(proteoforms) == 0){
    return(data.frame(Start_Index = query[1], End_Index = query[2]))
  }

  merged <- FALSE
  for(i in seq(1, nrow(proteoforms))){
    current <- c(proteoforms$Start_Index[i], proteoforms$End_Index[i])
    test <- is_overlap(query, current)
    if(length(test) == 2){
      # Vector is returned so perform a merge
      proteoforms$Start_Index[i] <- test[1]
      proteoforms$End_Index[i] <- test[2]
      
      # In version 1, we'll only perform one merge per cycle
      merged <- TRUE
      break
    }
  }
  
  if(merged) {
    # if a merge was performed, recur with new proteoform list
    new_query <- c(proteoforms$Start_Index[1], proteoforms$End_Index[1])
    return(make_contigs(new_query, proteoforms[-1,]))
  } else {
    # if no merge was performed, hold onto the old query and recur with rest
    new_query <- c(proteoforms$Start_Index[1], proteoforms$End_Index[1])
    new_prot <- make_contigs(new_query, proteoforms[-1,])
    new_prot[nrow(new_prot) + 1,] <- query
    return(new_prot)
  }
}

#' Sum length of contigs
#'
#' @param df Dataframe containing the Start and End indices of proteoform contigs
#'
#' @return Length of all contigs
sum_contigs <- function(df) {
  df <- df %>%
    mutate(Len = End_Index - Start_Index + 1)
  return(sum(df$Len, na.rm = TRUE))
}
