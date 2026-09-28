

#Written by Sam Majoros and Tyler A. Elliott


#working environment
setwd("D:/Documents_prune/ALAB_R/Sam/CDIPTERA_2026/")
#Setting up Data----

#Load packages
library(seqinr)
library(dplyr)
library(ggplot2)
library(viridis)
library(grDevices)
library(reshape2)
library(janitor)
library(stringr)
library(readr)
library(tidyverse)
library(gridExtra)
library(ggtext)

#Read in the fasta file 
SeqData <- read.fasta("Diptera_fasta_Mar172026.fas")

#Pull out the additional information 
SeqInfo <- getAnnot(SeqData)
#Pull of the sequence
Seqs <- getSequence(SeqData)

#Combine into a dataframe
SeqInfo <- do.call("rbind", SeqInfo)
SeqInfo <- as.data.frame(SeqInfo)
SeqInfo <- data.frame(do.call("rbind", strsplit(as.character(SeqInfo$V1), "|", fixed = TRUE)))

#Concatenate the sequence vectors
Seqs <- lapply(1:length(Seqs), function(i){
  paste(Seqs[[i]], collapse = "")
})
#Add sequences to the dataframe 
SeqInfo$nucleotides <- Seqs
#Unlist the column
SeqInfo$nucleotides <- unlist(SeqInfo$nucleotides)

colnames(SeqInfo)

#Rename columns 

names(SeqInfo)[names(SeqInfo) == "X1"] <- "processid"
names(SeqInfo)[names(SeqInfo) == "X2"] <- "Name"

#clean processid column
SeqInfo$processid<-gsub(">", "", SeqInfo$processid)


#read in the specimen information
SpecInfo<-read.delim("DS_CDIPETRA_spec_info_Mar172026.txt", sep='\t')

#join to SeqInfo

SeqInfo<-left_join(SeqInfo, SpecInfo, by='processid')

#read in the BIN uniqueness data from BOLD
UniqueBIN<-read.delim("unique_bin_Mar172026.txt", sep='\t')
#remove duplicate rows
UniqueBIN %>% distinct()->UniqueBIN

#join to SeqInfo
SeqInfo<-left_join(SeqInfo, UniqueBIN, by ='BIN')

#Unique BINs raw
SeqInfo %>% filter(Unique =="YES") %>% dplyr::select(BIN) %>% distinct() %>% tally()
#129 BINs
#unique BIN specimens raw
SeqInfo %>% filter(Unique =="YES") %>% distinct() %>% tally()
#260 specimens

#set empty data in rows and columns to NA

SeqInfo <- SeqInfo %>%
  mutate(
    across(is.character, ~na_if(., ""))
  )

#add in the EOSD land use classification data
EOSD<-read.delim("EOSD_sites.txt")
#make lat and lon columns numeric
EOSD <- EOSD %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon))
#ensure same columns in SeqInfo are numeric
SeqInfo <- SeqInfo %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon))
#join
SeqInfo<-left_join(SeqInfo, EOSD, by=c("lat", "lon"))

#standardization of Collectors
#summarize and standardize
SeqInfo %>% dplyr::select(Collectors) %>% distinct()->Collectors

SeqInfo$Collectors <- gsub("(?i)malaise|malaise trap|a-frame malaise 2", "Malaise", SeqInfo$Collectors)

#J Wang Jinjing Wang
SeqInfo$Collectors <- gsub("(?i)J Wang", "Jinjing Wang", SeqInfo$Collectors)
#Arctic & Boreal Entomology Arctic Ecology Arctic Ecology 2009 Arctic Ecology 2010 Boreal & Arctic Entomology course
SeqInfo$Collectors <- gsub("(?i)Arctic & Boreal Entomology|Arctic Ecology|Arctic Ecology 2009|Arctic Ecology 2010|Boreal & Arctic Entomology course", "Field Course", SeqInfo$Collectors)
# E.Stur E. Stur
SeqInfo$Collectors <- gsub("(?i)E. Stur", "E.Stur", SeqInfo$Collectors)
#E.Sturt and T.Ekreum T. Ekrem, E.Stur E.Stur, T. Ekrem
SeqInfo$Collectors <- gsub("(?i)E.Sturt and T.Ekreum|T. Ekrem, E.Stur|E.Stur, T. Ekrem", "E.Stur, T.Ekrem", SeqInfo$Collectors)
#R.E.Roughley R. E. Roughley R.E. Roughley
SeqInfo$Collectors <- gsub("(?i)R.E. Roughley|R. E. Roughley", "R.E.Roughley", SeqInfo$Collectors)
#P.D.N. Hebert P.Hebert
SeqInfo$Collectors <- gsub("(?i)P.Hebert", "P.D.N. Hebert", SeqInfo$Collectors)
#B Laforest Brandon Laforest
SeqInfo$Collectors <- gsub("(?i)Brandon Laforest", "B Laforest", SeqInfo$Collectors)
#J.Skevington J. Skevington
SeqInfo$Collectors <- gsub("(?i)J. Skevington", "J.Skevington", SeqInfo$Collectors)
#Emily Noel E Noele
SeqInfo$Collectors <- gsub("(?i)E Noele", "Emily Noel", SeqInfo$Collectors)
#P. Kevan/T. Woodcock P. Kevan & T. Woodcock
SeqInfo$Collectors <- gsub("(?i)P. Kevan & T. Woodcock|P. Kevan/T. Woodcock", "P. Kevan, T. Woodcock", SeqInfo$Collectors)
#A. Wiskin Andrew Wiskin
SeqInfo$Collectors <- gsub("(?i)A. Wiskin", "Andrew Wiskin", SeqInfo$Collectors)
#J.Dwaard J.deWaard
SeqInfo$Collectors <- gsub("(?i)J.Dwaard", "J.deWaard", SeqInfo$Collectors)
SeqInfo$Collectors <- gsub("(?i)J.Skevington,A.Renaud", "J.Skevington, A. Renaud", SeqInfo$Collectors)
SeqInfo$Collectors <- gsub("(?i)J.Lankshear, G.Lundie, J.McGowan", "J.Lankshear, J.McGowan, G.Lundie", SeqInfo$Collectors)
SeqInfo$Collectors <- gsub("(?i)A. Thielman", "A.Thielman", SeqInfo$Collectors)
SeqInfo$Collectors <- gsub("(?i)A. Renaud", "A.Renaud", SeqInfo$Collectors)
SeqInfo$Collectors <- gsub("(?i)J.M.Cumming", "JM Cumming", SeqInfo$Collectors)

SeqInfo  %>% dplyr::select(Collectors) %>% distinct()->Collectors #46 rows


#Numbers of specimens/collector
# CollectorSummary <- SeqInfo %>%
#   count(Collectors)

#Collectors summary

ggplot(data = SeqInfo, aes(x = Collectors)) +
  geom_bar() +
  theme_minimal() +
  xlab("Collector") +
  ylab("Count") +
  ggtitle("Number of Rows by Collector") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#SamplingProtocol standardization and finding this data from Notes in some cases

# Make SamplingProtocol2 based on matches in the Notes column
#Pan trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)Pantrap|Pan trap|pan trap|pan traps", SeqInfo$Notes), 
                                    "Pan Trap", 
                                    SeqInfo$SamplingProtocol)
#Sweeping variants
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)Handsweep|Hand sweep|sweep net|Sweep net|Sweeping|sweeping", SeqInfo$Notes), 
                                    "Sweep Net", 
                                    SeqInfo$SamplingProtocol)

SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)Sweeping", SeqInfo$SamplingProtocol), 
                                    "Sweep Net", 
                                    SeqInfo$SamplingProtocol)


#malaise trap, Malaise trap, Malaise Trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)malaise trap|Malaise trap|Malaise Trap|a-frame malaise 2|malaise|A-frame malaise 2", SeqInfo$Notes), 
                                    "Malaise", 
                                    SeqInfo$SamplingProtocol)

#UV trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)UV light", SeqInfo$Notes), 
                                    "UV Light Trap", 
                                    SeqInfo$SamplingProtocol)

#light trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)Light Trap", SeqInfo$Notes), 
                                    "Light Trap", 
                                    SeqInfo$SamplingProtocol)

#fan trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)fan trap", SeqInfo$Notes), 
                                    "Fan Trap", 
                                    SeqInfo$SamplingProtocol)

#baited trap
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)liver trap|Horse Fly trap", SeqInfo$Notes), 
                                    "Baited Trap", 
                                    SeqInfo$SamplingProtocol)

#pitfall
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)pitfall trap pitfall|pitfall|Pitfall", SeqInfo$Notes), 
                                    "Pitfall Trap", 
                                    SeqInfo$SamplingProtocol)

#netting
SeqInfo$SamplingProtocol2 <- ifelse(grepl("(?i)netting", SeqInfo$Notes), 
                                    "Dip Net", 
                                    SeqInfo$SamplingProtocol)



SeqInfo <- SeqInfo %>%
  mutate(
    SamplingProtocol2 =
      recode(
        SamplingProtocol2,
        "Sweeping" = "Sweep Net"
      )
  )

SeqInfo %>%
  filter(Collectors == "Jinjing Wang",
         str_detect(CollectionDate, "2010"))->dfJinjing0 #7679
dfJinjing0 %>% filter(SamplingProtocol2 =="Malaise Trap") %>% tally() #3854/7679 = 50% Malaise Trap
dfJinjing0 %>% filter(SamplingProtocol2 =="Hand Collection") %>% tally() #736/7679 = 9.6%
dfJinjing0 %>% filter(SamplingProtocol2 =="Sweep Net") %>% tally() #1774/7679 23.1%
dfJinjing0 %>% filter(is.na(SamplingProtocol2)) %>% tally() #410/7679 5.3%
dfJinjing0 %>% filter(SamplingProtocol2 =="Pan Trap") %>% tally() #900/7679 11.7%

dfJinjing0 %>% filter(SamplingProtocol2 =="Light Trap") %>% tally() #14/7679 0.18%

dfJinjing0 %>% dplyr::select(lat, lon) %>% distinct() %>% tally()

SeqInfo %>%
  
  filter(
    Collectors == "Jinjing Wang",
    str_detect(CollectionDate, "2010")
  ) %>%
  
  group_by(
    lat,
    lon,
    SamplingProtocol2
  ) %>%
  
  summarise(
    n_records = n(),
    .groups = "drop_last"
  ) %>%
  
  mutate(
    percent =
      100 *
      n_records /
      sum(n_records)
  ) %>%
  
  ungroup() ->
  
  dfJinjing4


dfJinjing4 %>% filter(is.na(SamplingProtocol2)) %>% tally() #39

dfJinjing4 %>% filter(is.na(SamplingProtocol2)) %>% summary(percent)
dfJinjing4 %>% filter(SamplingProtocol2=="Malaise Trap") %>% tally() #20
dfJinjing4 %>% filter(SamplingProtocol2=="Pan Trap") %>% tally() #17
dfJinjing4 %>% filter(SamplingProtocol2=="Hand Collection") %>% tally() #40
dfJinjing4 %>% filter(SamplingProtocol2=="Sweep Net") %>% tally() #42
dfJinjing4 %>% filter(SamplingProtocol2=="Light Trap") %>% tally()



dfJinjing4_summary <-
  dfJinjing4 %>%
  
  group_by(SamplingProtocol2) %>%
  
  summarise(
    
    min_records = min(n_records, na.rm = TRUE),
    max_records = max(n_records, na.rm = TRUE),
    mean_records = mean(n_records, na.rm = TRUE),
    median_records = median(n_records, na.rm = TRUE),
    
    min_percent = min(percent, na.rm = TRUE),
    max_percent = max(percent, na.rm = TRUE),
    mean_percent = mean(percent, na.rm = TRUE),
    median_percent = median(percent, na.rm = TRUE),
    
    .groups = "drop"
    
  )

write_tsv(dfJinjing4_summary, "Jinjing_2010_collection_summary3.tsv")

write_tsv(dfJinjing4, "Jinjing_2010_collection_summary2.tsv")

SeqInfo %>%
  filter(Collectors == "Jinjing Wang",
         str_detect(CollectionDate, "2010")) %>%
  group_by(lat, lon, SamplingProtocol2, CollectionDate) %>%
  summarise(n_dates = n_distinct(CollectionDate),
            n_records =n(),
            .groups = "drop")->dfJinjing


dfJinjing3 <-
  SeqInfo %>%
  
  filter(
    Collectors == "Jinjing Wang",
    str_detect(CollectionDate, "2010")
  ) %>%
  
  group_by(
    lat,
    lon,
    CollectionDate,
    SamplingProtocol2
  ) %>%
  
  summarise(
    
    n_records = n(),
    
    .groups = "drop_last"
    
  ) %>%
  
  mutate(
    
    n_protocols =
      n_distinct(
        SamplingProtocol2
      )
    
  ) %>%
  
  ungroup()

dfJinjing3 <-
  dfJinjing3 %>%
  
  mutate(
    CollectionDate =
      as.Date(
        CollectionDate,
        format = "%d-%b-%Y"
      )
  ) %>%
  
  arrange(
    lat,
    CollectionDate
  )

sampling_overview <- dfJinjing3 %>%
  mutate(
    SamplingProtocol2 = replace_na(SamplingProtocol2, "NA")
  ) %>%
  pivot_wider(
    id_cols = c(lat, lon, CollectionDate),
    names_from = SamplingProtocol2,
    values_from = n_records,
    values_fill = 0,
    values_fn = sum
  ) %>%
  mutate(
    site_id = paste0(dense_rank(paste(lat, lon)), "_siteid")
  ) %>%
  relocate(site_id, .before = lat)


site_summary <- sampling_overview %>%
  group_by(site_id, lat, lon) %>%
  summarise(
    n_collection_dates = n_distinct(CollectionDate),
    .groups = "drop"
  )

site_summary

#write_tsv(dfJinjing, "Jinjing_2010_collection_summary.tsv")

SeqInfo %>%
  filter(Collectors == "Jinjing Wang",
         str_detect(CollectionDate, "2010")) %>%
  filter(SamplingProtocol2=="Malaise Trap") %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(Species),
            n_bins = n_distinct(BIN),
            n_records =n(),
            .groups = "drop")->dfJinjing2

dfJinjing2 %>% dplyr::mutate(seq_intensity = n_records / n_species)->dfJinjing2


#raw specimen count
summary(SeqInfo)
#16786
#raw species count
SeqInfo %>% filter(!is.na(Species)) %>% dplyr::select(Species) %>% distinct() %>% tally()
#1216 species
#without species
SeqInfo %>% filter(is.na(Species)) %>% tally()
#4067 specimens without species
#raw genus count
SeqInfo %>% filter(!is.na(Genus)) %>% dplyr::select(Genus) %>% distinct() %>% tally()
#446 genera
#without genus
SeqInfo %>% filter(is.na(Genus)) %>%  tally()
#1199 specimens without genus 

#raw family count
SeqInfo %>% filter(!is.na(Family)) %>% dplyr::select(Family) %>% distinct() %>% tally()
#68 families
SeqInfo %>% filter(is.na(Family)) %>% tally()
#0 specimens without a family 

#BIN count
SeqInfo %>% filter(!is.na(BIN)) %>% dplyr::select(BIN) %>% distinct() %>% tally()
#2225 BINs
#specimens without BINs
SeqInfo %>% filter(is.na(BIN)) %>% tally()
#169 specimens without BINs

#institution data
SeqInfo %>% group_by(Storage.Institution) %>%
  summarise(rec_count = n())
  







#Filter out records with incorrect coordinates 
SeqInfo2 <- SeqInfo %>%
  filter(lat > 50 | lat < 40)
#16779

#Merge, Split, Matches Figure-----

#Filter out records that were not assigned to a BIN
SeqInfo2 <- SeqInfo2 %>%
  filter(!is.na(BIN)) #16610
#Remove records with no species name 
SeqInfo2 <- SeqInfo2 %>%
  filter(!is.na(Species)) #12608

#Filter out species with only one sequence
SeqInfo2 <- SeqInfo2 %>%
  group_by(Species) %>%
  filter(n() > 1) #12293

#BIN number

SeqInfo2 %>% ungroup() %>% dplyr::select(BIN) %>% distinct() %>% tally() #983 BINs

#Find number of records per species
RecordCount <- SeqInfo %>% filter(!is.na(Species)) %>%
  group_by(Species) %>%
  dplyr::count()

#Find number of BINs each species is in 
BINCount <- SeqInfo %>% filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarise(count = n_distinct(BIN))

#Find count of species found in more than one BIN (Splits)
length(which(BINCount$count >1)) #143

BINCount %>% filter(count >1) %>% dplyr::select(Species) %>% distinct() %>% tally()


#BINs with no species
SeqInfo %>% filter(is.na(Species)) %>% filter(!is.na(BIN)) %>% dplyr::select(BIN) %>%
  distinct() %>% tally()

#Count of BINs in multi-BIN species
BINCount %>% filter(!is.na(Species)) %>% filter(count >1) %>%
  summarise(total = sum(count)) #335

#Find how many species is in each BIN
SpecCount <- SeqInfo %>% filter(!is.na(BIN)) %>% filter(!is.na(Species)) %>%
  group_by(BIN) %>%
  summarise(count = n_distinct(Species))

#Filter for BINs with more than one species
OverOne <- filter(SpecCount, SpecCount$count >1) #40 rows
#Find the specimens in these BINs
OverOneData <- SeqInfo[SeqInfo$BIN %in% OverOne$BIN, ] #1506 rows
#Find the number of unique species (Merges)
OverOneData %>% filter(!is.na(Species)) %>% dplyr::select(Species) %>% distinct() %>% tally() #85


#multi-species BINs
length(unique(OverOneData$BIN)) #40

#Find how many matches there are
BIN_match <- filter(BINCount, BINCount$count <2)
Spec_match <- filter(SpecCount, SpecCount$count <2)
#Select only the BINs and species found in these two dataframes
match <- SeqInfo[SeqInfo$Species %in% BIN_match$Species, ]
match <- match[match$BIN %in% Spec_match$BIN, ]

#Unique rows
match %>% dplyr::select(BIN, Species) %>% distinct() %>% tally() #990

#matches where species are represented by more than on sequence

match2 <- match %>% 
  group_by(Species) %>%
  summarise(
    n_rows = n(),
    .groups = "drop"
  ) 

match2 %>% filter(n_rows >1) %>% tally()



#Plot a bar graph 
#counts <- c(823, 46, 70) #hard coded numbers
counts <- c(990, 85, 143)
Categories <- c("Matches", "Merges", "Splits")
#Create table
plotdata <- data.frame(counts, Categories)
#Create barplot
plotdata %>%
  ggplot(aes(x=counts, y= reorder(Categories, counts))) +
  geom_bar(stat="identity", fill= viridis(3)) +
  xlab("Counts") +
  ylab("Category") +
  theme_bw()


#Other Summary Stats-----

#Find number of unique species
length(unique(SeqInfo$Species))


#Finding BINs that do not have any species id
#Remove sequences without a BIN
SeqInfoBIN <- SeqInfo %>%
  filter(!grepl(">", BIN))
#Group by BIN and filter for BINs with only one sequence
NoSpecies <- SeqInfoBIN %>%
  group_by(BIN) %>%
  filter(length(unique(Species)) < 2)
NoSpecies <- ungroup(NoSpecies)
#Find BINs with no species name
NoSpecies <- filter(NoSpecies, NoSpecies$Species == "")
#Find number of BINs with no species name
length(unique(NoSpecies$BIN))


#Filter down main dataset to just these species 
NoSpecies <- filter(SeqInfoBIN, SeqInfo$BIN %in% NoSpecies$BIN)
#Find number of genus per BIN 
GenusNumber <- NoSpecies %>%
  group_by(BIN) %>%
  summarise(count = n_distinct(Genus))
#They all only have one. 
#Reduce dataset down to just genus and BIN
OneGenus <- (NoSpecies[, c("Genus","BIN")])
#Reduce down to unique rows
OneGenus <- distinct(OneGenus)
#Find number of BINs without a genus id
length(which(!(OneGenus$Genus == "")))

#Find number of families per BIN 
FamilyNumber <- NoSpecies %>%
  group_by(BIN) %>%
  summarise(count = n_distinct(Family))
#They all only have one. 
#Reduce dataset down to just family and BIN
OneFamily <- (NoSpecies[, c("Family","BIN")])
#Reduce down to unique rows
OneFamily <- distinct(OneFamily)
#Find number of BINs without a genus id
length(which(!(OneFamily$Family == "")))




#Sampling over years----

#Create a figure showing the collection occurring in each region over time
#Filter out records without collection info 
dfCollection <- SeqInfo %>%
  filter(!CollectionDate == "")

#Split collection data so the year is in its own column
dates <- data.frame(do.call("rbind", strsplit(as.character(dfCollection$CollectionDate), "-", fixed = TRUE)))
#Merge dates with the rest of the data
dfCollection <- cbind(dfCollection, dates)
#Rename columns 
names(dfCollection)[names(dfCollection) == "X1"] <- "Day"
names(dfCollection)[names(dfCollection) == "X2"] <- "Month"
names(dfCollection)[names(dfCollection) == "X3"] <- "Year"

#Find count for each Year
YearCount <- dfCollection %>%
  group_by(Year) %>%
  dplyr::count()

#Plot counts
ggplot(data=YearCount, aes(x=Year, y=n)) +
  geom_bar(stat="identity", fill= viridis(7)) +
  xlab("Year") +
  ylab("Number of Records") +
  theme_bw()




#Canada Map----

#Modified from: https://teng.pub/technical/2020/1/7/drawing-canada-maps-in-r 

library(sf)
#library(rgdal)
library(geojsonio)
#library(spdplyr)
library(rmapshaper)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(magrittr) 
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggplot2)



# Lambert Conformal Conic projection for Canada
canada_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lon_0=-91.52 +datum=NAD83 +units=m +no_defs"

# Canada provinces
canada <- ne_states(
  country = "canada",
  returnclass = "sf"
)

# Transform projection
canada <- st_transform(canada, crs = canada_crs)

# Create province label coordinates
province_labels <- canada %>%
  st_point_on_surface() %>%
  cbind(st_coordinates(.))

# Churchill coordinates
churchill <- data.frame(
  name = "Churchill",
  lon = -94.164963,
  lat = 58.68410
)

churchill_sf <- st_as_sf(
  churchill,
  coords = c("lon", "lat"),
  crs = 4326
)

# Transform Churchill point
churchill_sf <- st_transform(
  churchill_sf,
  crs = canada_crs
)

# Extract projected coordinates for Churchill label
churchill_coords <- churchill_sf %>%
  cbind(st_coordinates(.))

# ---- Rotated prairie/west labels ----

rotated_labels <- province_labels %>%
  filter(name %in% c(
    "British Columbia",
    "Alberta",
    "Saskatchewan",
    "Manitoba"
  ))

# ---- Atlantic provinces ----

atlantic_points <- data.frame(
  
  # leader line anchor points
  lon = c(-57, -64, -64, -66),
  lat = c(53, 45, 46.3, 46.5),
  
  # adjusted label positions
  label_lon = c(-54, -63, -62, -69),
  label_lat = c(55, 42.8, 48, 43.5),
  
  name = c(
    "Newfoundland\nand Labrador",
    "Nova Scotia",
    "Prince Edward\nIsland",
    "New Brunswick"
  )
)

# Convert line anchor points
atlantic_sf1 <- st_as_sf(
  atlantic_points,
  coords = c("lon", "lat"),
  crs = 4326
) %>%
  st_transform(crs = canada_crs)

atlantic_sf2 <- st_as_sf(
  atlantic_points,
  coords = c("label_lon", "label_lat"),
  crs = 4326
) %>%
  st_transform(crs = canada_crs)

atlantic_labels <- cbind(
  atlantic_points["name"],
  st_coordinates(atlantic_sf1),
  st_coordinates(atlantic_sf2)
)

colnames(atlantic_labels) <- c(
  "name",
  "point_x",
  "point_y",
  "label_x",
  "label_y"
)

# ---- Remaining labels ----

regular_labels <- province_labels %>%
  filter(!name %in% c(
    "British Columbia",
    "Alberta",
    "Saskatchewan",
    "Manitoba",
    "Newfoundland and Labrador",
    "Nova Scotia",
    "Prince Edward Island",
    "New Brunswick"
  ))

# Northwest Territories line break
regular_labels <- regular_labels %>%
  mutate(
    name = ifelse(
      name == "Northwest Territories",
      "Northwest\nTerritories",
      name
    )
  )

# Manual label adjustments
regular_labels <- regular_labels %>%
  mutate(
    
    # move Ontario upward
    Y = ifelse(
      name == "Ontario",
      Y + 200000,
      Y
    ),
    
    # move Quebec left
    X = ifelse(
      name == "Québec",
      X - 200000,
      X
    ),
    
    # shift Northwest Territories left
    X = ifelse(
      name == "Northwest\nTerritories",
      X + 100000,
      X
    )
    
  )

# ---- Plot ----

#Figure 1b

ggplot() +
  
  # provinces
  geom_sf(
    data = canada,
    aes(fill = name),
    color = "grey40",
    linewidth = 0.2,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  
  # regular labels
  geom_text(
    data = regular_labels,
    aes(X, Y, label = name),
    size = 3,
    fontface = "bold"
  ) +
  
  # rotated prairie/west labels
  geom_text(
    data = rotated_labels,
    aes(X, Y, label = name),
    angle = 90,
    size = 3,
    fontface = "bold"
  ) +
  
  # Atlantic leader lines
  geom_segment(
    data = atlantic_labels,
    aes(
      x = point_x,
      y = point_y,
      xend = label_x,
      yend = label_y
    ),
    linewidth = 0.4
  ) +
  
  # Atlantic labels
  geom_text(
    data = atlantic_labels,
    aes(label_x, label_y, label = name),
    hjust = 0,
    size = 2.5,
    fontface = "bold"
  ) +
  
  # Churchill point
  geom_sf(
    data = churchill_sf,
    color = "red",
    size = 3
  ) +
  
  # Churchill label
  geom_text_repel(
    data = churchill_coords,
    aes(
      X,
      Y,
      label = name
    ),
    size = 4,
    fontface = "bold",
    nudge_y = 300000
  ) +
  
  coord_sf(crs = canada_crs) +
  
  theme_minimal() +
  
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 40,
      b = 10,
      l = 10
    )
  )


#now lets make Figure 1a
# Arctic projection
arctic_crs <- "+proj=laea +lat_0=90 +lon_0=-100"

# World map
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Canada, Greenland, Norway, Finland
highlight <- world %>%
  filter(admin %in% c(
    "Canada",
    "Greenland",
    "Norway",
    "Finland"
  )) %>%
  mutate(
    Region = admin
  )

# Alaska only
alaska <- world %>%
  filter(admin == "United States of America") %>%
  st_crop(
    xmin = -170,
    xmax = -130,
    ymin = 50,
    ymax = 75
  ) %>%
  mutate(
    Region = "Alaska"
  )

# Combine
highlight <- bind_rows(
  highlight,
  alaska
)

# Transform projection
world_proj <- st_transform(world, crs = arctic_crs)
highlight_proj <- st_transform(highlight, crs = arctic_crs)

# Plot
ggplot() +
  
  # Background countries
  geom_sf(
    data = world_proj,
    fill = "grey92",
    color = "black",
    linewidth = 0.2
  ) +
  
  # Highlighted regions
  geom_sf(
    data = highlight_proj,
    aes(fill = Region),
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    values = c(
      "Alaska" = "#2b6cb0",
      "Canada" = "#ef6351",
      "Greenland" = "#ecec72",
      "Norway" = "#5db7ae",
      "Finland" = "#d17ac6"
    ),
    
    breaks = c(
      "Alaska",
      "Canada",
      "Greenland",
      "Norway",
      "Finland"
    ),
    
    name = "Regions"
  ) +
  
  coord_sf(
    crs = arctic_crs,
    xlim = c(-7000000, 7000000),
    ylim = c(-7000000, 5000000)
  ) +
  
  theme_minimal() +
  
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    panel.grid.major = element_line(
      color = "white",
      linewidth = 0.3
    ),
    
    panel.background = element_rect(
      fill = "grey85",
      color = NA
    ),
    
    legend.title = element_text(face = "bold"),
    
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )


#Figure 1C, of the sampling sites around Churchill



library(mapview)

#Filter out records without coordinates 
MapData <- SeqInfo %>%
  filter(!lat == "") %>%
  filter(!Species == "")

#convert coords to numeric
MapData$lat <- as.numeric(MapData$lat)
MapData$lon <- as.numeric(MapData$lon)

#Plot Figure 1c
mapview(MapData, xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE)

#need PID, Lat, Lon, Exact_Site

SpecInfo %>% dplyr::select(processid, lat, lon, ExactSite)->MapInfo

#set Lat/Lon to numeric
MapInfo$lat <- as.numeric(MapInfo$lat)
MapInfo$lon <- as.numeric(MapInfo$lon)
#filter records without coords
MapInfo %>% filter(!lat=="") %>% filter(!lon=="")->MapInfo
#filter out incorrect coords
MapInfo %>% filter(lat >50) ->MapInfo

#make summary info

#making a new one that summarizes specimen number by unique gps site and exact site info
MapInfo %>% 
  dplyr::group_by(lat, lon, ExactSite) %>% dplyr::summarise(specimen_count=n()) ->MapInfoSiteGroups

#give the sitegroups an index that makes searching on the map easier 
MapInfoSiteGroups$index <- 1:nrow(MapInfoSiteGroups)

#write_tsv(MapInfoSiteGroups, "MapInfoSiteGroups.tsv")




#lets read in the land use classification data for our coordinates


#association between land use codes and Earth Observation for Sustainable Development categories
LandClassLeg<-read.delim("EOSD_land_use_legend.txt")


#land use codes associated to our GPS coordinates via ARCGIS
LandClassCoord<-read.delim("Churchill_coords_luc.txt")

#join together
LandClass<-inner_join(LandClassCoord,LandClassLeg, by="landuse" )
#subset the data

LandClass %>% dplyr::select(Lat, Lon, Description....RG)->LandClass

#rename some columns
LandClass %>% dplyr::rename(Luc_Classification = Description....RG)->LandClass
LandClass %>% dplyr::rename(lat = Lat)->LandClass
LandClass %>% dplyr::rename(lon = Lon)->LandClass

LandClass %>% distinct()->LandClass #267 entries


#Lets massage our coordinates a bit before we join them to our data 


#Lets join to our  MapInfoSiteGroups

MapInfoSiteGroups2 <- MapInfoSiteGroups %>%
  mutate(
    lat = round(lat, 4),
    lon = round(lon, 4)
  ) %>%
  full_join(
    LandClass %>%
      mutate(
        lat = round(lat, 4),
        lon = round(lon, 4)
      ),
    by = c("lat", "lon")
  )

#two additional rows

#Make a new column with standardized names

MapInfoSiteGroups2 %>% mutate (Site2 = case_when((Luc_Classification == "Bryoids" ~"Bryoids"),
                                                 (Luc_Classification == "Coniferous-open" ~"Coniferous"),
                                                 (Luc_Classification == "Coniferous-sparse" ~"Coniferous"),
                                                 (Luc_Classification == "Developed?" ~"Exposed"),
                                                 (Luc_Classification == "Exposed/Barren Land" ~"Exposed"),
                                                 (Luc_Classification == "Rock/Rubble" ~"Exposed"),
                                                 (Luc_Classification == "Shrub Tall" ~"Shrubland"),
                                                 (Luc_Classification == "Water" ~"Water"),
                                                 (Luc_Classification == "Wetland-herb" ~"Wetland"),
                                                 (Luc_Classification == "Wetland-shrub" ~"Wetland"),
                                                 (Luc_Classification == "Wetland-treed" ~"Wetland")))->MapInfoSiteGroups2

#summarize again

MapInfoSiteGroups3 <- MapInfoSiteGroups2 %>%
  group_by(lat, lon, Site2) %>%
  mutate(
    specimen_count = sum(specimen_count, na.rm = TRUE)
  ) %>%
  ungroup()


#investigate and consolidate any duplicates

MapInfoSiteGroups3 %>% dplyr::select(lat, lon, specimen_count, Site2) %>% get_dupes()


#write_csv(MapInfoSiteGroups3, "EOSD_site_classification2.csv")

#For visualization purposes we will exclude certain sites that are far away or with too few or no records
MapInfoSiteGroups3 %>% filter(!is.na(index))->MapInfoSiteGroups3


#Figure 1C upper part
mapview(MapInfoSiteGroups3, xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE)

#some sites we need to exclude for the map because they have too few records
#Lat      Lon
#58.7536 -94.9128
#58.1248 -92.9530
#57.6130 -93.9635
#57.5758 	-93.8702
#58.6570 -93.2240
#58.6630 -93.1910
#58.6620 -93.1910
#58.6180 -94.8110 1 record
#58.7340 -93.1090 1 record
#58.3720 -93.4950 1 record 
#58.1300 -92.9000 1 record
#58.1290 -92.8650 2 records
#58.1270 -92.8650 20 records
#58.1280 -92.8610 14 records
#58.1280 -92.8560 3 records  


MapInfoSiteGroups3 %>% filter(lat != 58.7536  & lon != -94.9128) %>%
  filter(lat != 58.1248  & lon != -92.9530) %>%
  filter(lat != 57.6130  & lon != -93.9635) %>%
  filter(lat != 57.5758  & lon != -93.8702) %>%
  filter(lat != 58.6570  & lon != -93.2240) %>%
  filter(lat != 58.6630  & lon != -93.1910) %>%
  
  filter(lat != 58.6180  & lon != -94.8110) %>%
  filter(lat != 58.7340  & lon != -93.1090) %>%
  filter(lat != 58.3720  & lon != -93.4950) %>%
  filter(lat != 58.1300  & lon != -92.9000) %>%
  filter(lat != 58.1290  & lon != -92.8650) %>%
  filter(lat != 58.1280  & lon != -92.8610) %>%
  filter(lat != 58.1280  & lon != -92.8560) %>%
  
  filter(lat != 58.6620  & lon != -93.1910) -> MapInfoSiteGroups3_subset

library(mapview)
library(sp)
library(sf)
#library(rgdal) 
library(raster)
library(ggplot2)
library(ggpattern)
library(gridExtra)
library(leaflet.extras)
library(ggrepel)


#Figure 1C lower part
#mapview(MapInfoSiteGroups3, cex = 2, crs = 4269, grid = TRUE, xcol = "Lon", ycol = "Lat", zcol ="inferred", col.regions=color_palette, layer.name ="Classification Source") +
mapview(MapInfoSiteGroups3_subset, cex = 10, crs = 4269, grid = FALSE, xcol = "lon", ycol = "lat", zcol ="Site2", alpha =3, layer.name ="Classification" )



#Supplementary Figure 1, number specimens from each year


SeqInfo %>%
  filter(!is.na(CollectionDate)) %>%
  mutate(
    CollectionDate = as.Date(CollectionDate, format = "%d-%b-%Y"),
    Year = as.integer(format(CollectionDate, "%Y"))
  ) %>%
  dplyr::count(Year) %>%
  ggplot(aes(x = factor(Year), y = n, fill = factor(Year))) +
  geom_col() +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(
    x = "Year",
    y = "Number of Specimens"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Hybotidae
# 502 samples
# 31 BINs
# 16 species, 12 interim names
# 61 no species
# 
# Dolichopodidae
# 898 samples 
# 71 BINs
# 32 species, 35 interim names 
# 3 no species
# 
# Empididae
# 438 samples
# 41 BINs
# 17 species, 21 interim names 
# 1 with no species




#Collection Methods----



#Find count for each collection method

SeqInfo <- SeqInfo %>%
  mutate(
    SamplingProtocol2 = case_when(
      SamplingProtocol2 %in% c("malaise", "A-frame malaise 2") ~ "Malaise Trap",
      SamplingProtocol2 == "pan traps" ~ "Pan Trap",
      SamplingProtocol2 == "Sweeping" ~ "Sweep Net",
      TRUE ~ SamplingProtocol2
    )
  )

MethodData <- SeqInfo %>%
  group_by(SamplingProtocol2) %>%
  dplyr::count()

MethodData  <- MethodData  %>%
  mutate(SamplingProtocol2 = ifelse(is.na(SamplingProtocol2), "Unknown", SamplingProtocol2))

#Order the counts
MethodData <- MethodData[order(MethodData$n), ]


#Create barplot for Figure 2
ggplot(data = MethodData, aes(x = n, y = reorder(SamplingProtocol2, n), fill = n)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "viridis") +
  xlab("Log Specimen Count") +
  ylab("Collection Method") +
  scale_x_log10() +
  theme_bw()

#identification

IDMethod<- SeqInfo %>%
  group_by(IDMethod, Identifier) %>%
  dplyr::count()

SeqInfo <- SeqInfo %>%
  mutate(
    Identifier = case_when(
      Identifier %in% c("BOLD Identification System []", "BOLD Identification System [2017]",
                        "BOLD Identification System [2019]", "BOLD Identification System [2016]",
                        "BOLD Identification System [2020]", "BOLD Identification System [2023]") ~ "BOLD ID Engine",
      TRUE ~ Identifier
    )
  )


IDPerson <-SeqInfo %>%
  group_by(Identifier) %>%
  dplyr::count()



  


#Number of records and species per family----

#Create graph showing number of records per family
#How many families are in the dataset
length(unique(SeqInfo$Family)) #68
#Find count for each family
FamilyCount <- SeqInfo %>%
  group_by(Family) %>%
  dplyr::count()

SeqInfo %>% filter(Family=="Chironomidae") %>% filter(!is.na(Species)) %>% distinct() %>% tally()
#1681

#Order the counts
FamilyCount <- FamilyCount[order(FamilyCount$n), ]


#plot for Figure 3a
ggplot(FamilyCount, aes(x = n, y = reorder(Family, n), fill = n)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "viridis", guide ="none") +
  xlab("Number of Records") +
  ylab("Family") +
  theme_bw()


family_protocol_counts <- SeqInfo %>%
  ungroup() %>%
  mutate(
    SamplingProtocol2 = replace_na(SamplingProtocol2, "NA")
  ) %>%
  dplyr::count(SamplingProtocol2, Family, name = "n") %>%
  arrange(SamplingProtocol2, desc(n))

family_protocol_counts


#Filter out rows with no species name
# SeqInfo <- SeqInfo %>%
#   filter(!Species == "")
#Find number of species per family
SpecCount <- SeqInfo %>% filter(!is.na(Species)) %>%
  group_by(Family) %>%
  summarise(count = n_distinct(Species))

#Order the counts
SpecCount <- SpecCount[order(SpecCount$count), ]


#plot for Figure 3b
ggplot(SpecCount, aes(x = count, y = reorder(Family, count), fill = count)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "viridis", guide ="none") +
  xlab("Number of Species") +
  ylab("Family") +
  theme_bw()

#Unique BINs per Family----

#How many families are in the dataset
length(unique(SeqInfo$Family)) #68


#Find BIN count for each family
FamBINCount <- SeqInfo %>% filter(!is.na(BIN)) %>%
  group_by(Family) %>%
  summarise(count = n_distinct(BIN))

#Unique BIN info from BOLD
#Mar 17th 2026
#raw data says 129 Unique BINs with 260 members, includes ones without species names and singletons
#2096 non-unique BINs including singletons and ones with species names

UniqBINFam<-SeqInfo %>% filter(!is.na(BIN)) %>%
  group_by(Family, Unique) %>%
  summarise(count = n_distinct(BIN))

UniqBINFam <- UniqBINFam %>%
  mutate(
    Unique = toupper(Unique),
    Unique = recode(Unique,
                    "YES" = "Unique",
                    "NO"  = "Non-Unique")
  )

UniqBINFam <- UniqBINFam[order(UniqBINFam$count), ]

UniqBINFam2 <- UniqBINFam %>%
  group_by(Family) %>%
  mutate(total_count = sum(count)) %>%
  ungroup()

#Plot Figure 4
ggplot(UniqBINFam2, aes(x = count,
                        y = reorder(Family, total_count),
                        fill = Unique)) +
  geom_bar(stat = "identity", color = "black") +
  xlab("Number of BINs") +
  ylab("Family") +
  theme_bw() +
  scale_fill_manual(
    name = "BIN type",
    values = c(
      "Non-Unique" = "#40E0D0", #40E0D0 turquoise  #0072B2 blue
      "Unique" = "#FF00FF"  #FF00FF magenta #D55E00 orange
    )
  )



##Supplementary figure of BIN number versus specimen number for each species 
#SpecInfo$BIN, SpecInfo$Species

SpecInfo$Species[
  trimws(SpecInfo$Species) == ""
] <- NA

SpecInfo_summary <- SpecInfo %>%
  group_by(Species) %>%
  summarise(
    n_BINs = n_distinct(BIN, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  )
SpecInfo_summary %>% filter(!is.na(Species))->SpecInfo_summary

SpecInfo_summary %>% filter(n_rows > 1)->SpecInfo_summary2


ggplot(
  SpecInfo_summary2,
  aes(
    x = n_rows,
    y = n_BINs
  )
) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  geom_smooth(
    method = "lm",
    se = TRUE, 
    color = 'orange'
  ) +
  labs(
    x = "Number of Specimens",
    y = "Number of BINs"
  )


model <- lm(
  n_BINs ~ n_rows,
  data = SpecInfo_summary
)

summary(model)

#lets filter out rows where n_BINS is 5 or greater and where n_rows is 200 or greater
#to see if this changes the relationship

SpecInfo_summary %>% filter(n_BINs <= 4) %>% filter(n_rows <= 200)-> SpecInfo_summary2

ggplot(
  SpecInfo_summary2,
  aes(
    x = n_rows,
    y = n_BINs
  )
) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Number of Specimens",
    y = "Number of BINs"
  )

model2 <- lm(
  n_BINs ~ n_rows,
  data = SpecInfo_summary2
)

summary(model2)




#Species Accumulation Curve----

library(iNEXT)
library(ggpubr)

#Plot accumulation curve for BINs
#Filter out records that were not assigned to a BIN
BINCurve <- SeqInfo %>%
  filter(!is.na(BIN))
#Create a new dataframe with BINs
dfCurve <- (BINCurve[, c("BIN", "Species")])
#Add column with the number of each records of each BIN at each site
dfCurve <- dfCurve %>%
  group_by(BIN) %>%
  mutate(count=n())
#Remove uneeded column
dfCurve <- (dfCurve[, c("BIN", "count")])
#Filter down to only unique rows 
dfCurve <- filter(unique(dfCurve))
#Plot the Curve
R_Curve_BIN <- iNEXT(dfCurve$count, q=0, datatype="abundance")
# ggiNEXT(R_Curve_BIN) + theme_bw()

BINCurvePlot<-ggiNEXT(R_Curve_BIN) +
  theme_bw() +
  xlab("Number of Individuals") +
  ylab("BIN Count")
#Calculate the Chao Richness
ChaoRichness(dfCurve$count)

#     Observed Estimator Est_s.e. 95% Lower 95% Upper
# 1     2225   3343.872  102.766  3159.897  3564.052

#Plot species accumulation curve for Species 
SpeciesCurve <- SeqInfo %>%
  filter(!is.na(Species))
dfCurve2 <- (SpeciesCurve[, c("BIN", "Species")])
#Add column with the number of each records of each BIN at each site
dfCurve2 <- dfCurve2 %>%
  group_by(Species) %>%
  mutate(count=n())
#Remove uneeded column
dfCurve2 <- (dfCurve2[, c("Species", "count")])
#Filter down to only unique rows 
dfCurve2 <- filter(unique(dfCurve2))
#Plot the Curve
R_Curve_Species <- iNEXT(dfCurve2$count, q=0, datatype="abundance")
# ggiNEXT(R_Curve_Species,) + theme_bw()

SpecCurvePlot<-ggiNEXT(R_Curve_Species) +
  theme_bw() +
  xlab("Number of Individuals") +
  ylab("Species Count")
#Calculate the Chao Richness
ChaoRichness(dfCurve2$count)

#    Observed Estimator Est_s.e. 95% Lower 95% Upper
# 1     1216  1524.559   45.365  1447.665  1626.975

#Plot the two curves side by side
# figure <- ggarrange(ggiNEXT(SpecCurvePlot), ggiNEXT(BINCurvePlot),
#                     labels = c("A", "B"))

#now plot Figure 5
figure <- ggarrange(
  SpecCurvePlot, BINCurvePlot,
  labels = c("A", "B"),
  ncol = 2, nrow = 1
)
plot(figure)


#Range of Distances----

#Read in data


BarcodeGapData <- read.delim("barcode_gap_analysis_Mar182026.txt", sep='\t')


#Remove records with NAs
BarcodeGapData <- BarcodeGapData %>%
  filter(!Mean.Intra.Sp == "N/A")
#Plot a histogram
hist(as.numeric(BarcodeGapData$Mean.Intra.Sp), col=rgb(0,0,1,0.4), xlab = "Genetic Distance (%)", ylab = "Frequency")
#Add the 2nd variable
hist(as.numeric(BarcodeGapData$Distance.to.NN), col = rgb(1,0,0,0.2), breaks = seq(from=0, to=25, by=1), add=TRUE)
#Add a legend
legend('topright', c('Mean Intra-specific Distance', 'Distance to Nearest Neighbour'),
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.2)))

#Distances by Family----

#Read in data
BarcodeGapData <- read.delim("barcode_gap_analysis_Mar182026.txt", sep='\t')



#Remove records with NAs
BarcodeGapData <- BarcodeGapData %>%
  filter(!Mean.Intra.Sp == "N/A")

#Separate dataframe by family 
Family_List <- split(BarcodeGapData, f = BarcodeGapData$Family)
#Find average NN distance
Distances <- lapply(1:length(Family_List), function(i){
  mean(Family_List[[i]]$Distance.to.NN)
})
#Combine values into one vector
Distances <- unlist(Distances)
#Get intra-specific divergence
Divergences <- lapply(1:length(Family_List), function(i){
  mean(as.numeric(Family_List[[i]]$Mean.Intra.Sp))
})
Divergences <- unlist(Divergences)
#Extract family names
FamNames <- names(Family_List)
#Create new dataframe
MeanFamValues <- data.frame(FamNames, Distances, Divergences)

#Format the values
MeanFamValues <- reshape2::melt(MeanFamValues, id.vars="FamNames")

#information to highlight certain families

highlight_fams <- c("Coelopidae", "Calliphoridae", "Canacidae",
                    "Ulidiidae", "Lonchopteridae", "Chamaemyiidae",
                    "Bolitophilidae")

#bolding them and colouring them black 
MeanFamValues <- MeanFamValues %>%
  mutate(
    FamLabel = ifelse(
      FamNames %in% highlight_fams,
      paste0("<b style='color:#000000;'>", FamNames, "</b>"),
      FamNames
    )
  )


#Plot the figure for Figure 6
ggplot(MeanFamValues,
       aes(x = value,
           y = reorder(FamLabel, value),
           fill = variable)) +
  geom_bar(stat = "identity", color = "black") +
  xlab("Genetic Distance (%)") +
  ylab("Family") +
  theme_bw() +
  scale_fill_manual(
    name = "Metric",
    values = c(
      "Distances" = "#40E0D0",
      "Divergences" = "#FF00FF"
    ),
    labels = c(
      "Distances" = "Mean Distance to\nNearest Neighbour",
      "Divergences" = "Mean Intraspecific\nDivergence"
    )
  ) +
  theme(
    axis.text.y = element_markdown()
  )

#Distance to NN vs Max Intra-Specific----



#Read in data
BarcodeGapData <- read.delim("barcode_gap_analysis_Mar182026.txt", sep='\t')



#Replicate the figure from BOLD for Figure 7
ggplot(BarcodeGapData, aes(x=Max.Intra.Sp, y=Distance.to.NN)) +
  geom_point(color = "#56B4E9") +
  geom_abline(colour="#0072B2", size=0.75) +
  xlab("Max Intraspecific Distance (%)") +
  ylab("Distance to Nearest Neighbour (%)") +
  theme_bw()

#Find number of points on each side of the line 
y_hat <- BarcodeGapData$Max.Intra.Sp * 1 + 0
Resid <- BarcodeGapData$Distance.to.NN - y_hat
PosResid <- Resid[Resid >= 0]
NegResid <- Resid[Resid < 0]
length(PosResid)
length(NegResid)



#Phylogenetic Tree----

library(ape)
library(dplyr)
library(foreach)
library(phangorn)
library(Biostrings)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DECIPHER")

library(DECIPHER)
library(tibble)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
library(ggtree)

#Filter for species with more than 10 records 
#Remove rows with no species name 
TreeInfo <- SeqInfo %>%
  filter(!is.na(Species))
#Filter so only top 50 species are included including ties 
TreeData <- TreeInfo %>%
  add_count(Species) %>%
  filter(dense_rank(-n) < 41)

TreeData %>% dplyr::select(Species) %>% distinct() %>% tally()
#51 species due to ties

#remove sequences with greater than 1% gaps or Ns
TreeData <- TreeData %>%
  mutate(
    seq_length = nchar(nucleotides),
    n_bad = stringr::str_count(nucleotides, "[Nn-]"),
    prop_bad = n_bad / seq_length
  ) %>%
  filter(prop_bad <= 0.01) #filters to 3694

TreeData %>% dplyr::select(Species) %>% distinct() %>% tally()
#51 species

TreeData %>% dplyr::select(Family, Species) %>% distinct()->TreeDataSumm

#set up study region
region <- "Churchill"
#Set up study Taxa
taxa <- "Diptera"

#Set up the model that should be used to select the best model for maximum likelihood tree generation.
TreeModel <- "BIC"

#Set up number of intervals of discrete gamma distribution (k-value) for the maximum likelihood tree generation to be used in the function pml.
K <- 4 

#Set distance matrix model to be used in the function dist.dna.
model <- "raw"

#Set up a function to trim the sequences using a reference sequence. 
RefSeqTrim <- function(x) {
  #Create data frame for reference sequence 
  #This reference sequence was taken from BOLD for Diptera (Record ID: ACGAZ1590-12, BIN:AAA1222).
  dfRefSeq <- data.frame(taxa= taxa, nucleotides= reference_sequence$Diptera)
  colnames(dfRefSeq)[2] <- "nucleotides"
  #Set nucleotides to a character 
  dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides)
  #Trim sequences to 620bp
  dfRefSeq$nucleotides <- substr(dfRefSeq$nucleotides, 20, nchar(dfRefSeq$nucleotides)-19)
  #Check sequence length
  dfRefSeq$seqLength <- nchar(dfRefSeq$nucleotides)
  #Ensure sequences are of character type
  alignmentSeqs <- as.character(x$nucleotides)
  #Name according to processid 
  names(alignmentSeqs) <- x$processid
  #Pull out reference sequence
  alignmentref <- as.character(dfRefSeq$nucleotides[1])
  #Name reference sequence 
  names(alignmentref) <- "Reference"
  #Put sequences together for alignment
  alignmentSeqsPlusRef <- append(alignmentref, alignmentSeqs)
  #Convert to DNAStringSet 
  DNAStringSet2 <- DNAStringSet(alignmentSeqsPlusRef)
  #Remove gaps
  DNAStringSet2 <- RemoveGaps(DNAStringSet2)
  #Run alignment 
  alignment2 <- AlignSeqs(DNAStringSet2)
  #Find stop and start positions in reference 
  refSeqPos <- which(alignment2@ranges@NAMES=="Reference")
  refSeqPos <- alignment2[refSeqPos]
  refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
  refSeqPosStart <- as.numeric(refSeqPosStart)
  refSeqPosEnd <- nchar(dfRefSeq$nucleotides[1]) + refSeqPosStart
  refSeqPosEnd <- as.numeric(refSeqPosEnd)
  #Trim sequence
  alignment2Trimmed <- substr(alignment2, refSeqPosStart, refSeqPosEnd -1)
  #Convert to DNAStringSet
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)
  #Remove reference sequence 
  refSeqRm <- which(DNAStringSet3@ranges@NAMES=="Reference")
  dnaStringSet3 <- subset(DNAStringSet3[-refSeqRm])
  alignmentOrder <- DNAStringSet3@ranges@NAMES
  #Reorder based on alignment
  x <- x[match(alignmentOrder, x$processid), ]
  #Replace old sequences with new ones 
  trimmedSeqs <- as.character(DNAStringSet3)
  x$nucleotides <- trimmedSeqs
  #Return dataframe with new sequences 
  return(x)
}

#Set up reference sequence
reference_sequence <- list(Diptera = "AACTTTATATTTTATTTTTGGAGCTTGATCTAGAATAATTGGAACTTCTTTAAGAATATTAATTCGAATTGAATTAGGTCATCCAGGTTCCTTAATTGGAAATGACCAAATTTATAATGTAATTGTAACAGCTCATGCATTTATTATAATTTTTTTTATAGTAATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTTCCTTTAATATTAGGAGCACCAGATATAGCTTTTCCTCGAATGAATAATATAAGTTTTTGACTTCTTCCTCCTGCTTTAATACTTTTATTAACAAGTAGAATAGTAGAAAGTGGAGCTGGAACAGGATGAACAGTTTATCCTCCTTTATCATCTATTATTGCTCATGGAGGAGCATCTGTTGACTTAGCTATTTTTTCTCTTCATTTAGCAGGAATTTCTTCTATTTTAGGAGCTGTAAATTTTATTACAACTGTAATTAATATACGATCTATTGGTATTACCTTTGATCGAATACCTTTATTTGTTTGATCAGTTGCTATTACAGCCTTATTACTTTTATTATCTTTACCAGTTTTAGCTGGAGCAATTACAATATTATTAACAGATCGAAATTTAAATACATCATTTTTTGATCCTGCTGGAGGAGGAGATCCTATTTTATACCAACATTTATTT")

#Remove blank space at beginning and end of sequences
trimseq <- trimws(TreeData$nucleotides, whitespace = "-", which = "left")
trimseq <- trimws(trimseq, whitespace = "-", which = "right")
#remove Ns at beginning and end of sequences
trimseq <- gsub("^[Nn-]+", "", TreeData$nucleotides)
trimseq <- gsub("[Nn-]+$", "", trimseq)
#Add column to dataframe
TreeData$nucleotides <- trimseq

#Choose a consistent length for sequence 
# TreeData <- filter(TreeData, nchar(TreeData$nucleotides) == 669)

#Split data by species
#Split dataframe based on Bin_uri so each BIN can be analyzed individually
Species_List <- split(TreeData, f = TreeData$Species)

#Find a centroid sequence for each species.
#Find the number of processids in each species 
Species_Size <- sapply(Species_List, function(x)length(x$processid))
#Create new data frame with species and species size
dfData_Species <- data.frame(Species_Size)
dfData_Species <- rownames_to_column(dfData_Species, "Species")
#Merge TreeData and dfData_Species
dfSpecies_List <- merge(TreeData, dfData_Species, by.x="Species", by.y="Species")
#Reorder dfData_Species by Species
dfData_Species <- dfData_Species[order(dfData_Species$Species), ]

#Find number of unique Species 
SpeciesNumber <- unique(dfSpecies_List$Species)   
SpeciesNumber <- length(SpeciesNumber) 

#Create list from dfCentroid 
largeSpeciesList <- lapply(unique(dfSpecies_List$Species), function(x) dfSpecies_List[dfSpecies_List$Species == x, ])

#Extract process Id from each Species
largeSpeciesProcessid <- sapply(largeSpeciesList, function(x) (x$processid))

#Convert sequences to dnaStringSet
dnaStringSet1 <- sapply(largeSpeciesList, function(x) DNAStringSet(x$nucleotides))  
#Name dnaStringSet with processids   
for(i in seq(from=1, to=SpeciesNumber, by=1)) {
  names(dnaStringSet1[[i]]) <- largeSpeciesProcessid[[i]]
}
#Check class of dnaStringSet1
class(dnaStringSet1[[1]])
#Check names of stringsets 
names(dnaStringSet1[[1]])

#Remove gaps
dnaStringSet1 <- lapply(1:length(dnaStringSet1), function(i){
  RemoveGaps(dnaStringSet1[[i]])
})

#Run multiple sequence alignment for sequences in each Species in dnaStringSet1
# alignment1 <- foreach(i=1:SpeciesNumber) %do%
#   AlignSeqs(dnaStringSet1[[i]])

alignment1 <- lapply(dnaStringSet1, function(x) {
  if (length(x) >= 2) {
    AlignSeqs(x)
  } else {
    x
  }
})




#Convert to DNAbin format 
dnaBINCentroid <- foreach(i=1:SpeciesNumber) %do% as.DNAbin(alignment1[[i]])
#Check class
class(dnaBINCentroid[[1]])

#Calculate a pairwise distance matrix for each species
geneticDistanceCentroid <- foreach(i=1:SpeciesNumber) %do%
  dist.dna(dnaBINCentroid[[i]], model= model, as.matrix = TRUE,
           pairwise.deletion = TRUE)

#Determine centroid sequence; The sequence with the minimum average distance to all other sequences in the species.
centroidSeq <- foreach(i=1:SpeciesNumber) %do% which.min(rowSums(geneticDistanceCentroid[[i]]))
centroidSeq <- centroidSeq %>%
  unlist() %>%
  names()

#Subset dfSpecies_List by the processid on the list 
dfSpecies_List <- subset(dfSpecies_List, processid %in% centroidSeq)

#Merge with the original data set 
dfAllSeq <- merge(dfSpecies_List, TreeData, by.x="processid", by.y="processid")
#Reorganize and clean up 
dfAllSeq <- (dfAllSeq[, c("processid", "BIN.x", "Species_Size", "Species.x", "nucleotides.x", "lat.x", "lon.x")])
colnames(dfAllSeq)[2] <- "BIN"
colnames(dfAllSeq)[4] <- "Species"
colnames(dfAllSeq)[5] <- "nucleotides"
colnames(dfAllSeq)[6] <- "lat"
colnames(dfAllSeq)[7] <- "lon"

#Rerun alignment with centroid sequences, outgroup and a reference
#Create RefSeq data frame
dfRefSeq <- data.frame(taxa= taxa, nucleotides= reference_sequence$Diptera)

#name nucleotide column and set as character 
colnames(dfRefSeq)[2] <- "nucleotides"
dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides) 
#Trim references to standard 620
dfRefSeq$nucleotides <- substr(dfRefSeq$nucleotides, 20, nchar(dfRefSeq$nucleotides)-19) 
#Check sequence length
dfRefSeq$seqLength <- nchar(dfRefSeq$nucleotides)

#Set up outgroups
Outgroups <- data.frame(processid = c("Outgroup1", "Outgroup2", "Outgroup3"),
                        BIN = c("Outgroup1", "Outgroup2", "Outgroup3"),
                        Species = c("Outgroup1", "Outgroup2", "Outgroup3"),
                        nucleotides = c("AACTTTATATTTTATTTTTGGAAGATGAGCTGGATTACTCGGAACATCATTAAGAGTTCTAATTCGAACAGAACTAGGAACTACTAAGTCTCTAATTAATAATGATCAAATTTATAACGTATTAGTAACTGCCCATGCTTTTATCATGATTTTTTTCATAGTGATACCTATTATAATTGGAGGATTTGGAAATTGACTTGTTCCACTAATACTCGGGGCCCCAGATATAGCATTCCCGCGAATAAATAATATAAGATTTTGACTCCTTCCCCCCTCTTTAAATTTTTTACTTCTTAGAAGGTTAGTAGAAAATGGTACAGGAACCGGTTGAACTGTTTATCCGCCATTAGCGAGAAATATCGGACATATGGGAAGCTCTGTTGATCTTTCAATTTTTTCCTTACATATGGCTGGGATTTCTTCTATTCTAGGAGCTATTAATTTTATTACAACTTGTATAAGAATAAAACCAAAAAGAATACTTTTAGACCAAATACCCTTATTTGTTTGATCAGTATTAATTACAGCTTTACTACTACTTCTATCTTTACCAGTTTTAGCCGGAGCTATCACTATATTACTTACAGACCGAAACCTAAATACATCTTTTTTTGACCCAGCAGGAGGCGGAGACCCCATTTTATACCAACATTTATTT","AACTTTATATTTTATTTTTGGAAGATGAGCTGGATTACTCGGAACATCATTAAGAGTTCTAATTCGAACAGAACTAGGAACTACTAAGTCTCTAATTAATAATGATCAAATTTATAACGTATTAGTAACTGCCCATGCTTTTATCATGATTTTTTTCATAGTGATACCTATTATAATTGGAGGATTTGGAAATTGACTTGTTCCACTAATACTCGGGGCCCCAGATATAGCATTCCCGCGAATAAATAATATAAGATTTTGACTCCTTCCCCCCTCTTTAAATTTTTTACTTCTTAGAAGGTTAGTAGAAAATGGTACAGGAACCGGTTGAACTGTTTATCCGCCATTAGCGAGAAATATCGGACATATGGGAAGCTCTGTTGATCTTTCAATTTTTTCCTTACATATGGCTGGGATTTCTTCTATTCTAGGAGCTATTAATTTTATTACAACTTGTATAAGAATAAAACCAAAAAGAATACTTTTAGACCAAATACCCTTATTTGTTTGATCAGTATTAATTACCGCTTTACTACTACTTCTATCTTTACCAGTTTTAGCCGGAGCTATCACTATATTACTTACAGACCGAAACCTAAATACATCTTTTTTTGACCCAGCAGGAGGCGGAGACCCCATTTTATACCAACATTTATTT","AACTTTATATTTTATTTTTGGAAGATGGGCTGGATTACTCGGAACATCATTAAGAGTTCTAATTCGAACAGAACTAGGAACTACTAAGTCTCTAATTAATGATGATCAAATTTATAACGTATTAGTAACTGCCCATGCTTTTATCATGATTTTTTTCATAGTGATACCTATTATAATTGGAGGATTTGGAAATTGACTTGTTCCACTAATACTCGGGGCCCCAGATATAGCATTCCCGCGAATAAATAATATAAGATTTTGACTCCTTCCCCCCTCTTTAAATTTTTTACTTCTTAGAAGGTTAGTAGAAAATGGTACAGGAACCGGTTGAACTGTTTATCCGCCATTAGCGAGAAATATCGGACATATGGGAAGCTCTGTTGATCTTTCAATTTTTTCCTTACATATGGCTGGGATTTCTTCTATTCTAGGAGCTATTAATTTTATTACAACTTGTATAAGAATAAAACCAAAAAGAATACTTTTAGACCAAATACCCTTATTTGTTTGATCAGTATTAATTACAGCTTTACTACTACTTCTATCTTTACCAGTTTTAGCCGGAGCTATCACTATATTACTTACAGACCGAAACCTAAATACATCTTTTTTTGACCCAGCAGGAGGCGGAGACCCCATTTTATACCAACATTTATTT"))

#Ensure bin_uri is a character
dfAllSeq$BIN <- as.character(dfAllSeq$BIN)
#Add outgroup sequences to the dataframe 
dfAllSeq <- full_join(dfAllSeq, Outgroups)

#Align and trim sequences 
dfAligned <- RefSeqTrim(dfAllSeq)

#Remove NAs 
dfAligned <- dfAligned[-c(1), ]

#Convert sequences to DNAStringSet format 
dnaStringSet2 <- DNAStringSet(dfAligned$nucleotides)
#Name the stringset 
names(dnaStringSet2) <- dfAligned$Species
#Check names are assigned properly 
names(dnaStringSet2)

#Export alignment to check it in another program such as MEGA. 

writeXStringSet(dnaStringSet2, file = "CentroidAlignment_ChurchillDiptera_moreTrim.fas", format = "fasta")

#Set up the data in the correct format and generate the maximum likelihood tree
#Read in alignment in phyDat format

phylo_dat <- read.phyDat("CentroidAlignment_ChurchillDiptera_moreTrim.fas", format="fasta", type="DNA")
#Check format
class(phylo_dat)

#create a distance matrix
dm <- dist.ml(phylo_dat)

#creating NJ tree
tree <- NJ(dm)

#Run a model test to determine which model to use when creating the maximum likelihood tree 
model_test <- modelTest(phylo_dat, tree)

#create the environment 
env <- attr(model_test, "env")

#create function to find best model 
get_best_model = function(model_df){
  best_model = model_df['Model'][model_df[TreeModel] == min(model_df[TreeModel]) ]
  return(best_model)
}
#Determine best model for the data
Best_Model <- get_best_model(model_test)

#Get parameters for the model 
model_fit <- eval(get(Best_Model, env), env)

#Pull out inv value 
inv_values <- model_fit$inv

#Compute likelihood. 
ml_out = pml(tree, phylo_dat, k=K, inv = model_fit$inv)

#Drop the suffix from the model name
Best_Model2 = unlist(strsplit(Best_Model, "\\+"))[[1]]

#Compute likelihood and optimize parameters.
ml_tree = optim.pml(ml_out, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = Best_Model2)

#Root the tree using the outgroups
Tree_Rooted <- root(ml_tree$tree, outgroup = Outgroups$processid, resolve.root = TRUE)
#Plot the tree
plot(Tree_Rooted)

#Remove the outgroups
Finaltree <- drop.tip(Tree_Rooted, Outgroups$processid)
#Plot the tree. 
plot(Finaltree)

#Colour the tree by Max Intra-specific distance 
#Filter down gap data to only include species in dfAligned
#Read in data
BarcodeGapData <- read.delim("barcode_gap_analysis_Mar182026.txt", sep='\t')
#Subset dfData_BINs down to number of BINs in dfCentroid
dfTree_Dist <- subset(BarcodeGapData, BarcodeGapData$Species %in% dfAligned$Species)
#Find sequences to be removed from tree
tobeRemoved <- subset(dfAligned, !dfAligned$Species %in% BarcodeGapData$Species)
#Filter down dataset
dfAligned <- subset(dfAligned, dfAligned$Species %in% BarcodeGapData$Species)
#Merge the dataframes
dfTree_Data <- merge(dfTree_Dist, dfAligned, by = "Species")

#Remove sequence from tree
Finaltree <- drop.tip(Tree_Rooted, tobeRemoved$Species)

#Ensure data is in same order as the tree 
dfTree_Data <- dfTree_Data[match(Finaltree$tip.label, dfTree_Data$Species),]
#Add label column 
dfTree_Data$label <- Finaltree$tip.label

dfTree_Data$Mean.Intra.Sp <- as.numeric(as.character(dfTree_Data$Mean.Intra.Sp))

#Add data to tree
Tree_plus_Data <- full_join(Finaltree, dfTree_Data, by = 'label')
#Get Family information

meta <- data.frame(
  label = Tree_plus_Data@extraInfo[["Species"]],  # or whatever column holds names
  Family = Tree_plus_Data@extraInfo[["Family"]]
)

tip_data <- data.frame(
  label = Tree_plus_Data@phylo$tip.label
) %>%
  left_join(meta, by = "label")


# Attach metadata to tree
# Base tree with metadata
p <- ggtree(Tree_plus_Data) %<+% tip_data

# Get tip positions
tree_df <- p$data %>%
  filter(isTip)

# Compute family ranges
family_ranges <- tree_df %>%
  left_join(tip_data, by = "label") %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  summarise(
    ymin = min(y),
    ymax = max(y),
    .groups = "drop"
  )

# Define tree width (before scaling)
xmax_tree <- max(p$data$x)

# Define annotation panel position
family_ranges <- family_ranges %>%
  mutate(
    xmin = xmax_tree * 1.05,
    xmax = xmax_tree * 1.12
  )

# Final plots for Figure 8
p +
  xlim_tree(xmax_tree * 0.7) +
  coord_cartesian(xlim = c(0, xmax_tree * 1.25), clip = "off") +
  
  geom_tippoint(aes(colour = Max.Intra.Sp), size = 2) +
  geom_tiplab(
    aes(color = Max.Intra.Sp),
    offset = 0.01,
    fontface = "bold", size = 3.5
  ) +
  geom_treescale(x = 0, y = 50) +
  
  scale_color_viridis_c(
    name = "Maximum\nintraspecific\ndistance (%)",
    option = "viridis"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "left",
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid = element_blank(),
    
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )



#Plot with mean values instead 


ggtree(Tree_plus_Data) +
  geom_tippoint(aes(colour=as.numeric(Mean.Intra.Sp)), size=2) +
  geom_tiplab(aes(color=as.numeric(Mean.Intra.Sp))) +
  geom_treescale(x = 0, y = 50) +
  scale_color_viridis() +
  theme(legend.position = "left")


p +
  xlim_tree(xmax_tree * 0.7) +
  coord_cartesian(xlim = c(0, xmax_tree * 1.25), clip = "off") +
  
  geom_tippoint(aes(colour = Mean.Intra.Sp), size = 2) +
  geom_tiplab(
    aes(color = Mean.Intra.Sp),
    offset = 0.01,
    fontface = "bold", size = 3.5
  ) +
  geom_treescale(x = 0, y = 50) +
  
  scale_color_viridis_c(
    name = "Mean\nintraspecific\ndistance (%)",
    option = "viridis"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "left",
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid = element_blank(),
    
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )



#Comparing Churchill data to other regions----

#Load packages
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)
library(readr)


#Filter out records that were not assigned to a BIN
SeqRegInfo <- SeqInfo %>%
  filter(!is.na(BIN))


#set temp directory for temporary files

#set a temporary R temp directory for large files
#adjust this as needed
# dir.create("D:/R_temp", showWarnings = FALSE)
# 
# Sys.setenv(VROOM_TEMP_PATH = "D:/R_temp")
# 
# #we will use the new BOLD API to get our data
# #https://boldsystems.org/data/api/ 
# 
# 
# #Alaska
# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Alaska")
# )
# 
# query_id <- content(res)$query_id
# query_id
# #eAFLT823KijKL8vMS07VLy5JLEm1csxJLM5OtC5JrLDKL0pJLbJyySwoSS1KtM7JzM0sSU0BAABLE7c= 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# download.file(json_url, "Alaska.json", mode = "wb")
# 
# #parse JSON to dataframe
# dfAlaska <- jsonlite::stream_in(file("Alaska.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfAlaska <- dfAlaska %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# 
# 
# dfAlaska <- dfAlaska %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# write_tsv(dfAlaska, "AlaskaData_Mar2026.tsv")

dfAlaska<-read.delim("AlaskaData_Mar2026.tsv")

#What BINs are shared between Alaska and our project? 
sharedBINs <- intersect(SeqRegInfo $BIN, dfAlaska$bin_uri)
#261 BINs shared between Churchill and Alaska 
#4416 records in Alaska

dfAlaska %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#918 total BINs
#261/918 = 0.2843137 in Alaska

SeqRegInfo %>% filter(!is.na(BIN)) %>% dplyr::select(BIN) %>% distinct() %>% tally()
#2225
#261/2225 0.1173034 in Churchill 

#filter to records with BINs
#and subset to specific columns for saving later
dfAlaska %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfAlaska
  




#Yukon
# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Yukon")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Yukon.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfYukon <- jsonlite::stream_in(file("Yukon.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfYukon <- dfYukon %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfYukon <- dfYukon %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# write_tsv(dfYukon, "YukonData_Mar2026.tsv")

dfYukon<-read.delim("YukonData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfYukon$bin_uri)
#1404 BINs shared between  
#110427 records in Alaska

dfYukon %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#6957 total BINs
#261/918 = 0.2843137 in Yukon

SeqRegInfo %>% filter(!is.na(BIN)) %>% dplyr::select(BIN) %>% distinct() %>% tally()
#2225
#261/2225 0.1173034 in Churchill

dfYukon %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfYukon





#Northwest Territories
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Northwest Territories")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("NWT.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfNWT <- jsonlite::stream_in(file("NWT.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfNWT <- dfNWT %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfNWT <- dfNWT %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# write_tsv(dfNWT, "NWTData_Mar2026.tsv")

dfNWT<-read.delim("NWTData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfNWT$bin_uri)
#1067 BINs shared between  
#42143 records in NWT

dfNWT %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#3761 total BINs

dfNWT %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfNWT



#Nunavut
# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Nunavut")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Nunavut.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfNunavut <- jsonlite::stream_in(file("Nunavut.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfNunavut <- dfNunavut %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfNunavut <- dfNunavut %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# write_tsv(dfNunavut, "NunavutData_Mar2026.tsv")
# 

dfNunavut<-read.delim("NunavutData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfNunavut$bin_uri)
#1067 BINs shared between  
#66626 records 

dfNunavut %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#1918 total BINs

dfNunavut %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfNunavut


#Greenland

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:country/ocean:Greenland")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Greenland.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfGreenland <- jsonlite::stream_in(file("Greenland.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfGreenland <- dfGreenland %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfGreenland <- dfGreenland %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# 
# write_tsv(dfGreenland, "GreenlandData_Mar2026.tsv")


dfGreenland<-read.delim("GreenlandData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfGreenland$bin_uri)
#95 BINs shared between  
#28303 records 

dfGreenland %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#338 total BINs

dfGreenland %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfGreenland

#Finland
# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:country/ocean:Finland")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Finland.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfFinland <- jsonlite::stream_in(file("Finland.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfFinland <- dfFinland %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfFinland  <- dfFinland  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# write_tsv(dfFinland, "FinlandData_Mar2026.tsv")

dfFinland<-read.delim("FinlandData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfFinland$bin_uri)
#440 BINs shared between  
#37272 records 

dfFinland %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#6725 total BINs

dfFinland %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfFinland


#Norway

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:country/ocean:Norway")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Norway.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfNorway <- jsonlite::stream_in(file("Norway.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfNorway <- dfNorway %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfNorway  <- dfNorway  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# write_tsv(dfNorway, "NorwayData_Mar2026.tsv")


dfNorway<-read.delim("NorwayData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfNorway$bin_uri)
#417 BINs shared between  
#106241 records 

dfNorway %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#6343 total BINs

dfNorway %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfNorway

#Manitoba for Southern Manitoba
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Manitoba")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Manitoba.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfManitoba <- jsonlite::stream_in(file("Manitoba.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfManitoba <- dfManitoba %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfManitoba <- dfManitoba %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# #filter for southern records
# dfManitoba %>% filter(!is.na(lat)) %>% filter(lat<54)->dfManitoba
# 
# 
# 
# 
# 
# write_tsv(dfManitoba, "SouthernManitobaData_Mar2026.tsv")

dfManitoba<-read.delim("SouthernManitobaData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfManitoba$bin_uri)
#3801 BINs shared between  
#9154 records 

dfManitoba %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#6343 total BINs

dfManitoba %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfManitoba

#Ontario for Southern Ontario
# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Ontario")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Ontario.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfOntario <- jsonlite::stream_in(file("Ontario.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfOntario  <- dfOntario  %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfOntario  <- dfOntario  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# #filter for southern records
# dfOntario  %>% filter(!is.na(lat)) %>% filter(lat<52)->dfOntario 
# 
# 
# 
# 
# 
# write_tsv(dfOntario, "SouthernOntarioData_Mar2026.tsv")

dfOntario<-read.delim("SouthernOntarioData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfOntario$bin_uri)
#936 BINs shared between  
#422415 records 

dfOntario %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#14934 total BINs


dfOntario %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfOntario

#British Columbia for Southern British Columbia

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:British Columbia")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("BC.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfBC <- jsonlite::stream_in(file("BC.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfBC  <- dfBC  %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfBC  <- dfBC  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# #filter for southern records
# dfBC  %>% filter(!is.na(lat)) %>% filter(lat<54)->dfBC
# 
# 
# 
# 
# 
# write_tsv(dfBC, "SouthernBCData_Mar2026.tsv")


dfBC<-read.delim("SouthernBCData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfBC$bin_uri)
#936 BINs shared between  
#141148 records 

dfBC %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#8392 total BINs

dfBC %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfBC

#Alberta for Southern Alberta

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Alberta")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Alberta.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfAlberta <- jsonlite::stream_in(file("Alberta.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfAlberta  <- dfAlberta  %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfAlberta  <- dfAlberta  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# #filter for southern records
# dfAlberta  %>% filter(!is.na(lat)) %>% filter(lat<54)->dfAlberta
# 
# 
# 
# 
# 
# write_tsv(dfAlberta, "SouthernAlbertaData_Mar2026.tsv")

dfAlberta<-read.delim("SouthernAlbertaData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfAlberta$bin_uri)
#1255BINs shared between  
#122481 records 

dfAlberta %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#8002 total BINs

dfAlberta %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfAlberta

#Saskatchewan for Southern Saskatchewan

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Saskatchewan")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("Saskatchewan.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfSaskatchewan <- jsonlite::stream_in(file("Saskatchewan.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfSaskatchewan  <- dfSaskatchewan  %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfSaskatchewan <- dfSaskatchewan  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# #filter for southern records
# dfSaskatchewan  %>% filter(!is.na(lat)) %>% filter(lat<54)->dfSaskatchewan
# 
# 
# 
# 
# 
# write_tsv(dfSaskatchewan, "SouthernSaskatchewanData_Mar2026.tsv")

dfSaskatchewan<-read.delim("SouthernSaskatchewanData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfSaskatchewan$bin_uri)
#742 BINs shared between  
#54741 records 

dfSaskatchewan %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#4382 total BINs

dfSaskatchewan %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfSaskatchewan

#New Brunswick

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:New Brunswick")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("NB.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfNB <- jsonlite::stream_in(file("NB.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfNB  <- dfNB  %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfNB <- dfNB  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# 
# 
# write_tsv(dfNB, "NBData_Mar2026.tsv")

dfNB<-read.delim("NBData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfNB$bin_uri)
#421 BINs shared between  
#30748 records 

dfNB %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#3208 total BINs

dfNB %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfNB

#Nova Scotia

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Nova Scotia")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("NS.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfNS <- jsonlite::stream_in(file("NS.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfNS  <- dfNS %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfNS <- dfNS  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# 
# 
# write_tsv(dfNS, "NSData_Mar2026.tsv")

dfNS<-read.delim("NSData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfNS$bin_uri)
#423 BINs shared between  
#47968 records 

dfNS %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#3208 total BINs

dfNS %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfNS


#Prince Edward Island

# 
# res <- GET(
#   "https://portal.boldsystems.org/api/query",
#   query = list(query = "tax:order:Diptera;geo:province/state:Prince Edward Island")
# )
# 
# query_id <- content(res)$query_id
# query_id
# 
# 
# #make the URL and download the JSON
# 
# json_url <- paste0(
#   "https://portal.boldsystems.org/api/documents/",
#   query_id,
#   "/download?format=json"
# )
# 
# options(timeout = 600)
# 
# res <- GET(
#   json_url,
#   write_disk("PEI.json", overwrite = TRUE),
#   progress()
# )
# 
# #parse JSON to dataframe
# dfPEI <- jsonlite::stream_in(file("PEI.json"), flatten = TRUE)
# 
# #clean out internal tabs and newline characters
# 
# dfPEI <- dfPEI %>%
#   mutate(across(where(is.character),
#                 ~ str_replace_all(., "[\t\r\n]", " ")))
# 
# #convert the coord column to lat and lon and convert to numeric
# 
# dfPEI <- dfPEI  %>%
#   mutate(
#     lat = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[1]),
#     lon = map_dbl(coord, ~ if (is.null(.x)) NA_real_ else .x[2])
#   )
# 
# 
# 
# 
# write_tsv(dfPEI, "PEIData_Mar2026.tsv")

dfPEI<-read.delim("PEIData_Mar2026.tsv")

#What BINs are shared with Churchill? 
sharedBINs <- intersect(SeqRegInfo$BIN, dfPEI$bin_uri)
#163 BINs shared between  
#17918 records 

dfPEI %>% filter(!is.na(bin_uri)) %>% dplyr::select(bin_uri) %>% distinct() %>% tally()
#3208 total BINs

dfPEI %>% filter(!is.na(bin_uri)) %>% 
  dplyr::select(processid, bin_uri, country.ocean, province.state)->dfPEI


#Combine all dataframes together to create supplementary file
AllRegions <- rbind(dfAlaska, dfNWT, dfNunavut, dfGreenland, dfFinland, dfNorway, dfManitoba, dfOntario, dfBC, dfAlberta, dfSaskatchewan, dfNB, dfNS, dfPEI)

#Save to a tsv file
#write_tsv(AllRegions,"BINsforRegionComparison.tsv")

#Create dataset of BINs shared with Churchill

AllRegions_filtered <- AllRegions %>%
  dplyr::filter(bin_uri %in% SeqInfo$BIN)

#write_tsv(AllRegions_filtered,"Supplementary_File_4_Shared_BINs_for_Region_Comparison.tsv")


#rm(dfAlaska, dfNWT, dfNunavut, dfGreenland, dfFinland, dfNorway, dfManitoba, dfOntario, dfBC, dfAlberta, dfSaskatchewan, dfNB, dfNS, dfPEI)

#Plot the number of BINs shared with our project 

#read in data for this plot, summarized from the dat above
plotdata<-as.data.frame(read.delim("SharedBIN_Mar2026.txt", sep='\t'))


#Add column that corrects for sampling in the region
plotdata$correctedvalue <- plotdata$sharedBINs/plotdata$NumberofRecords

#add columns to show the percentage of shared BINs based on Churchill and each Region

# get Churchill denominator
church_total <- plotdata %>%
  filter(Regions == "Churchill") %>%
  pull(totalBINs)

# add new columns
plotdata <- plotdata %>%
  mutate(
    sharedRegionpercent = sharedBINs / totalBINs,
    sharedChurchpercent = sharedBINs / church_total
  )

plotdata2 <- plotdata %>% filter(Regions!= "Churchill")

#Create barplot for Figure 10a
plotdata2 %>%
  ggplot(aes(x=sharedBINs, y= reorder(Regions, sharedBINs), fill=Category)) +
  geom_bar(stat="identity") +
  xlab("Number of Shared BINs") +
  ylab("Region") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")


#Figure 9a
plotdata2 %>%
  mutate(
    Region_label = paste0(Regions, " (", totalBINs, ")")
  ) %>%
  ggplot(aes(x = sharedRegionpercent,
             y = reorder(Region_label, sharedRegionpercent),
             fill = Category)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  xlab("Proportion of Shared BINs") +
  ylab("Region (Total BINs)") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")


#Plot with the corrected values
 plotdata2 %>%
   ggplot(aes(x=correctedvalue, y= reorder(Regions, correctedvalue), fill=Category)) +
   geom_bar(stat="identity") +
   xlab("Number of Shared BINs/Number of Specimens") +
   ylab("Region") +
   theme_bw() +
   scale_fill_brewer(palette = "Blues")



#Alpha and beta diversity----
library(vegan)



AllRegions<-as.data.frame(read.delim("BINsforRegionComparison.tsv", sep='\t'))

#groups for metrics

AllRegions <- AllRegions %>%
  mutate(Site = case_when(
    country.ocean %in% c("Greenland", "Finland", "Norway") ~ country.ocean,
    TRUE ~ province.state
  ))

#add in Churchill data
#extract relevant columns and rename them to fit other data
  SeqInfo %>% filter(!is.na(BIN)) %>% 
  dplyr::select(processid, BIN, Country, Province)->SeqInfoDiv
  
  SeqInfoDiv <-  SeqInfoDiv %>%
    dplyr::rename(
      bin_uri = BIN,
      country.ocean = Country,
      province.state = Province
    )
  
#add in Site column 
  SeqInfoDiv <-SeqInfoDiv %>%
    mutate(Site ="Churchill") %>%
    ungroup() %>%
    dplyr::select(processid, bin_uri, country.ocean, province.state, Site)
  
  #bind to AllRegions
  
  AllRegions<-rbind(AllRegions, SeqInfoDiv)
  
  AllRegions <-  AllRegions %>%
    dplyr::rename(
      Region = Site)
#rename Site to Region
  

#create a community matrix

comm <- AllRegions %>%
  dplyr::count(Region, bin_uri) %>%
  pivot_wider(
    names_from = bin_uri,
    values_from = n,
    values_fill = 0
  )

library(tibble)

comm_mat <- comm %>%
  column_to_rownames("Region") %>%
  as.matrix()

#extract country metadata

metadata <- AllRegions %>%
  distinct(Region, country.ocean)

metadata <- metadata[match(rownames(comm_mat), metadata$Region), ]

#regional diversity
alpha <- data.frame(
  Region = rownames(comm_mat),
  Richness = specnumber(comm_mat),
  Shannon = diversity(comm_mat, index = "shannon"),
  Simpson = diversity(comm_mat, index = "simpson")
)

#rarefaction for unequal sampling

min_depth <- min(rowSums(comm_mat))




set.seed(123)

n_iter <- 100   # number of rarefaction runs

rare_results <- replicate(n_iter, {
  comm_rare <- rrarefy(comm_mat, sample = min_depth)
  
  data.frame(
    Richness = specnumber(comm_rare),
    Shannon = diversity(comm_rare),
    Simpson = diversity(comm_rare, index = "simpson")
  )
}, simplify = FALSE)


#comm_rare <- rrarefy(comm_mat, sample = min_depth)

rare_df <- do.call(rbind, rare_results)
rare_df$Region <- rep(rownames(comm_mat), times = n_iter)


alpha_rare_summary <- rare_df %>%
  group_by(Region) %>%
  summarise(
    RarefiedRichness_mean = mean(Richness),
    RarefiedRichness_sd   = sd(Richness),
    RarefiedShannon_mean  = mean(Shannon),
    RarefiedShannon_sd    = sd(Shannon),
    RarefiedSimpson_mean  = mean(Simpson),
    RarefiedSimpson_sd    = sd(Simpson)
  )

alpha <- alpha %>%
  dplyr::left_join(alpha_rare_summary, by = "Region")
# 
# alpha$RarefiedRichness <- specnumber(comm_rare)
# alpha$RarefiedShannon <- diversity(comm_rare)


#relative abundance
comm_rel <- decostand(comm_mat, method = "total")

#Bray-Curtis, abundance based
bray <- vegdist(comm_rel, method = "bray")

#extract distances to Churchill

bray_mat <- as.matrix(bray)

churchill_dist <- bray_mat["Churchill", ]

churchill_df <- data.frame(
  Region = names(churchill_dist),
  BrayCurtis_to_Churchill = as.numeric(churchill_dist)
)

churchill_df <- churchill_df %>%
  filter(Region != "Churchill")




#Barplot for Figure10b

# plot_alpha_beta <- alpha %>%
#   select(Site, RarefiedShannon) %>%
#   left_join(churchill_df, by = "Site")
# 
# ggplot(plot_alpha_beta,
#        aes(x = RarefiedShannon, y = BrayCurtis_to_Churchill)) +
#   geom_point(size = 3) +
#   geom_text(aes(label = Site), vjust = -0.5) +
#   theme_bw() +
#   xlab("Rarefied Shannon Diversity") +
#   ylab("Beta Diversity relative to Churchill")
# 
# 
# #long format
# plot_long <- plot_alpha_beta %>%
#   pivot_longer(
#     cols = c(RarefiedShannon, BrayCurtis_to_Churchill),
#     names_to = "Metric",
#     values_to = "Value"
#   )
# 
# #faceted barplot
# ggplot(plot_long,
#        aes(x = Value,
#            y = reorder(Site, Value))) +
#   geom_col(fill = "steelblue") +
#   facet_wrap(~ Metric, scales = "free_x",
#              labeller = labeller(
#                Metric = c(
#                  BrayCurtis_to_Churchill = "Beta Diversity relative to Churchill",
#                  RarefiedShannon = "Rarefied Shannon Diversity"
#                )
#              )) +
#   theme_bw() +
#   xlab("Value") +
#   ylab("Site")





alpha_ci <- rare_df %>%
  group_by(Region) %>%
  summarise(
    mean = mean(Shannon),
    lower = quantile(Shannon, 0.025),
    upper = quantile(Shannon, 0.975)
  )

plot_df <- alpha_ci %>%
  left_join(churchill_df, by = "Region")





plot_long <- plot_df %>%
  pivot_longer(
    cols = c(mean, BrayCurtis_to_Churchill),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    lower_ci = ifelse(Metric == "mean", lower, NA),
    upper_ci = ifelse(Metric == "mean", upper, NA)
  )


category_df <- plotdata2 %>%
  dplyr::select(Region = Regions, Category) %>%
  bind_rows(
    data.frame(Region = "Churchill", Category = "Northern")
  )

plot_long <- plot_long %>%
  left_join(category_df, by = "Region")

 #Figure 9b
ggplot(plot_long,
       aes(x = Value,
           y = reorder(Region, Value),
           fill = Category)) +
  
  geom_col() +
  
  geom_errorbar(
    aes(xmin = lower_ci, xmax = upper_ci),
    width = 0.2,
    na.rm = TRUE
  ) +
  
  facet_wrap(~ Metric, scales = "free_x",
             labeller = labeller(
               Metric = c(
                 BrayCurtis_to_Churchill = "Beta Diversity relative to Churchill",
                 mean = "Rarefied Shannon Diversity (95% CI)"
               )
             )) +
  
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  xlab("Value") +
  ylab("Region")








#EOSD site classification summaries
#Figure 10

SeqRegInfo %>% dplyr::select(EOSDSiteClassification) %>% distinct()

# Wetland
# Exposed
# Water
# Bryoids
# Coniferous
# <NA>
# Shrubland

SeqRegInfo %>% dplyr::select(processid, Family, EOSDSiteClassification) %>%
  group_by(Family, EOSDSiteClassification) %>%
  summarise(specimen_count = n_distinct(processid))->FamilySiteRecordCount

FamilySiteRecordCount$specimen_count<-as.numeric(FamilySiteRecordCount$specimen_count)


SeqRegInfo %>% dplyr::select(processid, Family, EOSDSiteClassification) %>%
  group_by(EOSDSiteClassification) %>%
  summarise(specimen_count = n_distinct(processid))->SiteRecordCount

SiteRecordCount$specimen_count<-as.numeric(SiteRecordCount$specimen_count)

#Bryoids
FamilySiteRecordCount %>% ungroup() %>% filter(EOSDSiteClassification=="Bryoids") ->Bryoids
#top 10 families
Bryoids10 <- Bryoids %>%
  slice_max(order_by = specimen_count, n = 10)



#Families
Bryoids10 %>% arrange(desc(specimen_count)) %>% dplyr::select(Family)

#families and hexcode colours for them 
# 1 Chironomidae   "#C409A2"
# 2 Muscidae       "#E59B11"
# 3 Syrphidae      "#dc4242"
# 4 Culicidae      "#7E4E75"
# 5 Stratiomyidae  "#B931D8"
# 6 Tipulidae     "#8E0A43"
# 7 Limoniidae    "#EA0F6E"
# 8 Oestridae     "#e7ba24"
# 9 Sciomyzidae   "#7BF025"
# 10 Sphaeroceridae "#6E6036"

#make families as factors to make ordering them easier
Bryoids10$Family <- factor(Bryoids10$Family, levels = c("Chironomidae", "Muscidae",
                                                        "Syrphidae" , "Culicidae" ,
                                                        "Stratiomyidae" , "Tipulidae" ,
                                                        "Limoniidae" , "Oestridae" ,
                                                        "Sciomyzidae" , "Sphaeroceridae"  ))


bp<- ggplot(Bryoids10, aes(x="", y=specimen_count, fill=Family))+
  geom_bar(width = 1, stat = "identity")
#bp

pie <- bp + coord_polar("y", start=0) +  theme(axis.text = element_blank(),
                                               axis.ticks = element_blank(),
                                               panel.grid  = element_blank())
pie_data<-ggplot_build(pie)
pie_data$plot$data$Family
#make a named vector for colour to family matching
cols1 <- c("Chironomidae" = "#C409A2", "Muscidae" = "#E59B11",
           "Syrphidae" = "#dc4242", "Culicidae" = "#7E4E75",
           "Stratiomyidae" = "#B931D8", "Tipulidae" = "#8E0A43",
           "Limoniidae" = "#EA0F6E", "Oestridae" = "#e7ba24",
           "Sciomyzidae" = "#7BF025", "Sphaeroceridae" = "#6E6036")

#build the pie chart for Bryoids

brypie<-pie + scale_fill_manual(values =cols1, labels=c("Chironomidae", "Muscidae",
                                                        "Syrphidae" , "Culicidae" ,
                                                        "Stratiomyidae" , "Tipulidae" ,
                                                        "Limoniidae" , "Oestridae" ,
                                                        "Sciomyzidae" , "Sphaeroceridae" ),
                                breaks=c("Chironomidae", "Muscidae",
                                         "Syrphidae" , "Culicidae" ,
                                         "Stratiomyidae" , "Tipulidae" ,
                                         "Limoniidae" , "Oestridae" ,
                                         "Sciomyzidae" , "Sphaeroceridae" )) +
  ggtitle("Bryoid Sites")  +
  xlab("") +
  ylab("") +
  theme(legend.text = element_text(size=7), legend.title = element_text(size=10), legend.key.size = unit(0.4, 'cm'))

brypie

#Coniferous

FamilySiteRecordCount %>% ungroup() %>% filter(EOSDSiteClassification=="Coniferous") ->Coniferous
#top 10 families
Coniferous10 <- Coniferous %>%
  slice_max(order_by = specimen_count, n = 10)


#Families
Coniferous10 %>% arrange(desc(specimen_count)) %>% dplyr::select(Family)

#families and hexcode colours for them 

# 1 Cecidomyiidae  "#5d3d9c"
# 2 Muscidae       "#E59B11"
# 3 Mycetophilidae "#6942EE"
# 4 Dolichopodidae "#12B9D3"
# 5 Chironomidae   "#C409A2"
# 6 Chloropidae    "#98D50C"
# 7 Phoridae       "#12D3C7"
# 8 Pipunculidae   "#808080"
# 9 Anthomyiidae   "#D9740F"
# 10 Hybotidae  "#1771CF"

#make families as factors to make ordering them easier
Coniferous10$Family <- factor(Coniferous10$Family, levels = c("Cecidomyiidae", "Muscidae",
                                                  "Mycetophilidae", "Dolichopodidae", 
                                                  "Chironomidae", "Chloropidae",
                                                  "Phoridae", "Pipunculidae",
                                                  "Anthomyiidae", "Hybotidae" ))

bp2<- ggplot(Coniferous10, aes(x="", y=specimen_count, fill=Family))+
  geom_bar(width = 1, stat = "identity")
#bp

pie2 <- bp2 + coord_polar("y", start=0) +  theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.grid  = element_blank())
pie2_data<-ggplot_build(pie2)
pie2_data$plot$data$Family

#make a named vector for colour to family matching
cols2<- c("Cecidomyiidae" ="#5d3d9c", "Muscidae" ="#E59B11",
          "Mycetophilidae" ="#6942EE", "Dolichopodidae" ="#12B9D3", 
          "Chironomidae" ="#C409A2", "Chloropidae" ="#98D50C",
          "Phoridae" ="#12D3C7", "Pipunculidae" ="#808080",
          "Anthomyiidae" ="#D9740F", "Hybotidae" ="#1771CF" )


Conpie<-pie2 + scale_fill_manual(values =cols2, labels=c("Cecidomyiidae", "Muscidae",
                                                         "Mycetophilidae", "Dolichopodidae", 
                                                         "Chironomidae", "Chloropidae",
                                                         "Phoridae", "Pipunculidae",
                                                         "Anthomyiidae", "Hybotidae"), breaks=c("Cecidomyiidae", "Muscidae",
                                                                                                "Mycetophilidae", "Dolichopodidae", 
                                                                                                "Chironomidae", "Chloropidae",
                                                                                                "Phoridae", "Pipunculidae",
                                                                                                "Anthomyiidae", "Hybotidae")) + ggtitle("Coniferous Sites")  +
  xlab("") + ylab("")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=10), legend.key.size = unit(0.4, 'cm'))



Conpie

#Exposed


FamilySiteRecordCount %>% ungroup() %>% filter(EOSDSiteClassification=="Exposed") ->Exposed
#top 10 families
Exposed10 <- Exposed %>%
  dplyr::slice_max(order_by = specimen_count, n = 10)


#Families
Exposed10 %>% arrange(desc(specimen_count)) %>% dplyr::select(Family)

#families and hexcode colours for them 

# 1 Chironomidae   "#C409A2"
# 2 Muscidae       "#E59B11"
# 3 Chloropidae     "#98D50C"
# 4 Phoridae       "#12D3C7"
# 5 Dolichopodidae  "#12B9D3"
# 6 Syrphidae      "#dc4242"
# 7 Anthomyiidae   "#D9740F"
# 8 Cecidomyiidae  "#5d3d9c"
# 9 Ceratopogonidae "#F362D9"
# 10 Scathophagidae  "#F6C10E"

#make families as factors to make ordering them easier
Exposed10$Family <- factor(Exposed10$Family, levels = c("Chironomidae", "Muscidae",
                                                              "Chloropidae", "Phoridae", 
                                                              "Dolichopodidae", "Syrphidae",
                                                              "Anthomyiidae", "Cecidomyiidae",
                                                              "Ceratopogonidae", "Scathophagidae" ))

bp3<- ggplot(Exposed10, aes(x="", y=specimen_count, fill=Family))+
  geom_bar(width = 1, stat = "identity")
#bp

pie3 <- bp3 + coord_polar("y", start=0) +  theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.grid  = element_blank())
pie3_data<-ggplot_build(pie3)
pie3_data$plot$data$Family

#make a named vector for colour to family matching
cols3<- c("Chironomidae" ="#C409A2", "Muscidae" ="#E59B11",
          "Chloropidae" ="#98D50C", "Phoridae" ="#12D3C7", 
          "Dolichopodidae" ="#12B9D3", "Syrphidae" ="#dc4242",
          "Anthomyiidae" ="#D9740F", "Cecidomyiidae" ="#5d3d9c",
          "Ceratopogonidae" ="#F362D9", "Scathophagidae" ="#F6C10E")

Exppie<-pie3 + scale_fill_manual(values =cols3, labels=c("Chironomidae", "Muscidae",
                                                         "Chloropidae", "Phoridae", 
                                                         "Dolichopodidae", "Syrphidae",
                                                         "Anthomyiidae", "Cecidomyiidae",
                                                         "Ceratopogonidae", "Scathophagidae"), breaks=c("Chironomidae", "Muscidae",
                                                                                                        "Chloropidae", "Phoridae", 
                                                                                                        "Dolichopodidae", "Syrphidae",
                                                                                                        "Anthomyiidae", "Cecidomyiidae",
                                                                                                        "Ceratopogonidae", "Scathophagidae")) + ggtitle("Exposed Sites")  +
  xlab("") + ylab("")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=10), legend.key.size = unit(0.4, 'cm'))



Exppie

#water

FamilySiteRecordCount %>% ungroup() %>% filter(EOSDSiteClassification=="Water") ->Water
#top 10 families
Water10 <- Water %>%
  dplyr::slice_max(order_by = specimen_count, n = 10)


#Families
Water10 %>% arrange(desc(specimen_count)) %>% dplyr::select(Family)

#families and hexcode colours for them 

# 1 Muscidae       "#E59B11" 
# 2 Chironomidae  "#C409A2"
# 3 Dolichopodidae "#12B9D3"
# 4 Syrphidae     "#dc4242" 
# 5 Anthomyiidae   "#D9740F"
# 6 Chloropidae   "#98D50C"
# 7 Culicidae     "#7E4E75"
# 8 Sciomyzidae   "#7BF025"
# 9 Phoridae      "#12D3C7"
# 10 Scathophagidae "#F6C10E"

#make families as factors to make ordering them easier
Water10$Family <- factor(Water10$Family, levels = c("Muscidae", "Chironomidae",
                                                        "Dolichopodidae", "Syrphidae", 
                                                        "Anthomyiidae", "Chloropidae",
                                                        "Culicidae", "Sciomyzidae",
                                                        "Phoridae", "Scathophagidae" ))


bp4<- ggplot(Water10, aes(x="", y=specimen_count, fill=Family))+
  geom_bar(width = 1, stat = "identity")
#bp

pie4 <- bp4 + coord_polar("y", start=0) +  theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.grid  = element_blank())
pie4_data<-ggplot_build(pie4)
pie4_data$plot$data$Family

#make a named vector for colour to family matching
cols4<- c("Muscidae" ="#E59B11", "Chironomidae" ="#C409A2",
          "Dolichopodidae" ="#12B9D3", "Syrphidae" ="#dc4242", 
          "Anthomyiidae" ="#D9740F", "Chloropidae" ="#98D50C",
          "Culicidae" ="#7E4E75", "Sciomyzidae" ="#7BF025",
          "Phoridae" ="#12D3C7", "Scathophagidae" ="#F6C10E")

Watpie<-pie4 + scale_fill_manual(values =cols4, labels=c("Muscidae", "Chironomidae",
                                                         "Dolichopodidae", "Syrphidae", 
                                                         "Anthomyiidae", "Chloropidae",
                                                         "Culicidae", "Sciomyzidae",
                                                         "Phoridae", "Scathophagidae"), breaks=c("Muscidae", "Chironomidae",
                                                                                                 "Dolichopodidae", "Syrphidae", 
                                                                                                 "Anthomyiidae", "Chloropidae",
                                                                                                 "Culicidae", "Sciomyzidae",
                                                                                                 "Phoridae", "Scathophagidae")) + ggtitle("Water Sites")  +
  xlab("") + ylab("")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=10), legend.key.size = unit(0.4, 'cm'))



Watpie


#Wetland

FamilySiteRecordCount %>% ungroup() %>% filter(EOSDSiteClassification=="Wetland") ->Wetland
#top 10 families
Wetland10 <- Wetland %>%
  dplyr::slice_max(order_by = specimen_count, n = 10)


#Families
Wetland10 %>% arrange(desc(specimen_count)) %>% dplyr::select(Family)

#families and hexcode colours for them 
# 1 Muscidae      "#E59B11"
# 2 Chironomidae  "#C409A2"
# 3 Dolichopodidae"#12B9D3"
# 4 Anthomyiidae  "#D9740F"
# 5 Culicidae     "#7E4E75"
# 6 Simuliidae    "#914583"
# 7 Ephydridae    "#DED40A"
# 8 Syrphidae     "#dc4242"
# 9 Sciaridae     "#7649D1"
# 10 Chloropidae  "#98D50C"

#make families as factors to make ordering them easier
Wetland10$Family <- factor(Wetland10$Family, levels = c("Muscidae", "Chironomidae",
                                                    "Dolichopodidae", "Anthomyiidae", 
                                                    "Culicidae", "Simuliidae",
                                                    "Ephydridae", "Syrphidae",
                                                    "Sciaridae", "Chloropidae" ))

bp5<- ggplot(Wetland10, aes(x="", y=specimen_count, fill=Family))+
  geom_bar(width = 1, stat = "identity")
#bp

pie5 <- bp5 + coord_polar("y", start=0) +  theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.grid  = element_blank())
pie5_data<-ggplot_build(pie5)
pie5_data$plot$data$Family


#make a named vector for colour to family matching
cols5<- c("Muscidae" ="#E59B11", "Chironomidae" ="#C409A2",
          "Dolichopodidae" ="#12B9D3", "Anthomyiidae" ="#D9740F", 
          "Culicidae" ="#7E4E75", "Simuliidae" ="#914583",
          "Ephydridae" ="#DED40A", "Syrphidae" ="#dc4242",
          "Sciaridae" ="#7649D1", "Chloropidae" ="#98D50C")

Wetland_pie<-pie5 + scale_fill_manual(values =cols5, labels=c("Muscidae", "Chironomidae",
                                                              "Dolichopodidae", "Anthomyiidae", 
                                                              "Culicidae", "Simuliidae",
                                                              "Ephydridae", "Syrphidae",
                                                              "Sciaridae", "Chloropidae"), breaks=c("Muscidae", "Chironomidae",
                                                                                                    "Dolichopodidae", "Anthomyiidae", 
                                                                                                    "Culicidae", "Simuliidae",
                                                                                                    "Ephydridae", "Syrphidae",
                                                                                                    "Sciaridae", "Chloropidae")) + ggtitle("Wetland Sites")  +
  xlab("") + ylab("")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=10), legend.key.size = unit(0.4, 'cm'))



Wetland_pie


grid.arrange(brypie, Conpie, Exppie, Watpie, Wetland_pie,  nrow = 3, ncol=2)