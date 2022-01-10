library(tidyverse)
library(DBI)

## Helper functions
# Counts the number of fragments that bear a modification
modIsCovered <- function(mod_site, fragment_df) {
  fragments <- 0
  
  # check B ions
  frag_sites <- fragment_df %>% filter(IonTypeId == "B")
  if(nrow(frag_sites) > 0) {
    frag_sites <- frag_sites$CleavageSiteIndex
    # iterate through the modification sites and count the number of fragsites
    # that cover it
    fragments <- sapply(mod_site, function(x) sum(x <= frag_sites))
  }
  
  # check Y ions
  frag_sites <- fragment_df %>% filter(IonTypeId == "Y")
  if(nrow(frag_sites) > 0) {
    frag_sites <- frag_sites$CleavageSiteIndex
    fragments <- sapply(mod_site, function(x) sum(x >= frag_sites))
  }
  
  # list of fragment counts that support each modification
  #return(sum(fragments))
  return(fragments)
}

# Perform SQL query and count fragments
countSupportingFragments <- function(row){
  
  # get sites of sequence variance
  variants <- dbGetQuery(con, 'SELECT AppliedSequenceFeature.StartIndex, B.ChemicalProteoformId, B.Id
    FROM (SELECT * FROM Hit WHERE Hit.Id == :x)
    AS B LEFT JOIN BiologicalProteoform ON BiologicalProteoform.ChemicalProteoformId=B.ChemicalProteoformId
    LEFT JOIN AppliedSequenceFeature ON BiologicalProteoform.Id=AppliedSequenceFeature.BiologicalProteoformId
    WHERE AppliedSequenceFeature.ReplacementSequence IS NOT NULL', params = list(x = row['HitId']))
  
  if(nrow(variants)){
    varied <- TRUE
  } else {
    varied <- FALSE
  }
  
  # get sites of modifications from hit report table
  mods <- row['Modification_Codes']
  if(mods != ""){
    mods <- str_split(mods, pattern="\\|")[[1]]
    mods <- str_match(mods, "@([0-9]+)")[,2] # it seems the modification site is alway one less than it should be
    mods <- as.numeric(mods)
    modified <- TRUE
  } else {
    modified <- FALSE
  }
  
  # get sites of fragmentation
  matchingIon <- dbGetQuery(con, 'SELECT * FROM MatchingIon WHERE HitId == :x', params = list(x = row['HitId']))
  
  if(modified) {
    if(varied) {
      toReturn <- data.frame(varied = varied, modified = modified)
      toReturn$frags_per_mod <- list(modIsCovered(c(mods, variants$StartIndex),matchingIon))
    } else {
      toReturn <- data.frame(varied = varied, modified = modified)
      toReturn$frags_per_mod <- list(modIsCovered(mods, matchingIon))
    }
  } else {
    if(varied) {
      toReturn <- data.frame(varied = varied, modified = modified)
      toReturn$frags_per_mod <- list(modIsCovered(variants$StartIndex, matchingIon))
    } else{
      toReturn <- data.frame(varied = c(varied), modified = c(modified), frags_per_mod = 0)
    }
  }
  toReturn$num_fragments <- sapply(toReturn$frags_per_mod, min)
  return(toReturn)

}

# Import data
tissues <- c("Heart", "Kidney", "Lung", "SmInt", "Spleen")

heart_hits <- read.csv("tdreport_dumps/Heart_Hits_List_Props.csv")
kidney_hits <- read.csv("tdreport_dumps/Kidney_Hits_List_Props.csv")
lung_hits <- read.csv("tdreport_dumps/Lung_Hits_List_Props.csv")
smint_hits <- read.csv("tdreport_dumps/SmInt_Hits_List_Props.csv")
spleen_hits <- read.csv("tdreport_dumps/Spleen_Hits_list_props.csv")

# The hit report includes hits for proteins that do not pass FDR on the protein level
# We'll filter based on the PFRs listed on the Proteoform tab in TDviewer
heart_pfr_view <- read.csv("tdreport_dumps/Heart_PFR_View.csv")
kidney_pfr_view <- read.csv("tdreport_dumps/Kidney_PFR_View.csv")
lung_pfr_view <- read.csv("tdreport_dumps/Lung_PFR_View.csv")
smint_pfr_view <- read.csv("tdreport_dumps/SmInt_PFR_View.csv")
spleen_pfr_view <- read.csv("tdreport_dumps/Spleen_PFR_View.csv")

# Apply filter
heart_hits <- filter(heart_hits, PFR %in% heart_pfr_view$PFR)
kidney_hits <- filter(kidney_hits, PFR %in% kidney_pfr_view$PFR)
lung_hits <- filter(lung_hits, PFR %in% lung_pfr_view$PFR)
smint_hits <- filter(smint_hits, PFR %in% smint_pfr_view$PFR)
spleen_hits <- filter(spleen_hits, PFR %in% spleen_pfr_view$PFR)

# Heart search
con <- dbConnect(RSQLite::SQLite(), "tdreport_dumps/Heart.tdReport")
frag_support <- apply(heart_hits, 1, countSupportingFragments)
x <- do.call(rbind, frag_support)
heart_hits <- cbind(heart_hits, x)
dbDisconnect(con)

# Kidney search
con <- dbConnect(RSQLite::SQLite(), "tdreport_dumps/Kidney.tdReport")
frag_support <- apply(kidney_hits, 1, countSupportingFragments)
x <- do.call(rbind, frag_support)
kidney_hits <- cbind(kidney_hits, x)
dbDisconnect(con)

# Lung search
con <- dbConnect(RSQLite::SQLite(), "tdreport_dumps/Lung.tdReport")
frag_support <- apply(lung_hits, 1, countSupportingFragments)
x <- do.call(rbind, frag_support)
lung_hits <- cbind(lung_hits, x)
dbDisconnect(con)

# Small Intestine search
con <- dbConnect(RSQLite::SQLite(), "tdreport_dumps/SmInt.tdReport")
frag_support <- apply(smint_hits, 1, countSupportingFragments)
x <- do.call(rbind, frag_support)
smint_hits <- cbind(smint_hits, x)
dbDisconnect(con)

# Spleen search
con <- dbConnect(RSQLite::SQLite(), "tdreport_dumps/Spleen.tdReport")
frag_support <- apply(spleen_hits, 1, countSupportingFragments)
x <- do.call(rbind, frag_support)
spleen_hits <- cbind(spleen_hits, x)
dbDisconnect(con)


# Combine all results
hits <- rbind(heart_hits, kidney_hits, lung_hits, smint_hits, spleen_hits)

# Get metadata from the filename
hits <- hits %>%
  separate(col = File_Name, into = c('Date', 'NetId', 'Project', 'Separation', 'Tissue', 'Biorep', 'Fraction', 'Injection'), extra = "drop", convert = TRUE) %>%
  mutate(Separation = ifelse(Separation == "PLRPS", "LC", "CZE")) %>%
  mutate(Biorep = as.numeric(str_sub(Biorep, 7))) %>%
  mutate(Fraction = as.numeric(str_sub(Fraction, 5)))

save.image(file = "preprocessed.RData")
