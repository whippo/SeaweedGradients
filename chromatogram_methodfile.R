  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #                                                                                ##
  # Chromatogram Method File                                                       ##
  # Data are current as of 2020-10-19                                              ##
  # Data source: OIMB                                                              ##
  # R code prepared by Ross Whippo                                                 ##
  # Last updated 2020-10-19                                                        ##
  #                                                                                ##
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # SUMMARY:
  # Script to merge chromatagram .csv files together to generate manual entries for 
  # identifying peaks. 
  
  # Required Files (check that script is loading latest version):
  # *.csv files exported from OpenChrom
  
  # Associated Scripts:
  # FILE.R
  
  # TO DO 
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # TABLE OF CONTENTS                                                            ####
  #                                                                                 +
  # RECENT CHANGES TO SCRIPT                                                        +
  # LOAD PACKAGES                                                                   +
  # READ IN AND PREPARE DATA                                                        +
  # MANIPULATE DATA                                                                 +
  #                                                                                 +
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # RECENT CHANGES TO SCRIPT                                                     ####
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # 2020-10-19 Script Created
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # LOAD PACKAGES                                                                ####
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  library(tidyverse)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # READ IN AND PREPARE DATA                                                     ####
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  peakexport <- read.csv("/home/ross/Desktop/GradientsOpenChrom/std4_250ng.mL_342020_1_peakexport.csv")
  peaklist <- read.csv("/home/ross/Desktop/GradientsOpenChrom/std4_250ng.mL_342020_1_peaklist.csv")
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # MANIPULATE DATA                                                              ####
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # PEAK DETECTOR SCRIPT
  
  # separate out m.z. and intensity for each FA
  traces <- peakexport %>%
    select(Name, m.z, intensities)
  # generate temporary list of column names
  tracecols <- paste("Name", 1:999, sep="")
  # separate out m.z. values
  trace_mz <- separate(traces, m.z, into = tracecols, sep = " ")
  # remove excess from ':' onward in each m.z. column (retains only those values)
  mz_short <- as.data.frame(lapply(trace_mz[tracecols], gsub, pattern = ":.*", replacement = ""))
  # attach FA name column to mz
  mz_short$Name <- trace_mz$Name
  # make mz long
  mz_long <- mz_short %>%
    pivot_longer(-Name, names_to = "x", values_to = "m.z")
  # separate out intensities
  trace_int <- separate(traces, m.z, into = tracecols, sep = " ")
  # remove excess from ':' before in each m.z. column (retains only those values)
  int_short <- as.data.frame(lapply(trace_int[tracecols], gsub, pattern = ".*:", replacement = ""))
  # attach FA name column to int
  int_short$Name <- trace_int$Name
  # make int long
  int_long <- int_short %>%
    pivot_longer(-Name, names_to = "y", values_to = "int")
  # join long mz and int
  mz_int <- bind_cols(mz_long, int_long)
  # remove place holder columns
  mz_int <- mz_int %>%
    select(Name, m.z, int)
  # remove NAs
  mz_int <- mz_int %>%
    drop_na()
  # make int numeric
  mz_int$int <- as.numeric(as.character(mz_int$int))
  
  # create list of top 20 m.z for each FA
  FA_toptrace <- mz_int %>%
    arrange(desc(int)) %>%
    group_by(Name) %>%
    slice(1:20)
  
  # list of names for columns 
  FA_toptrace$trace <- "trace"
  num <- rep(1:20, n_distinct(FA_toptrace$Name))
  FA_toptrace$num <- num
  FA_toptrace <- FA_toptrace %>%
    unite(rank, trace, num)
  # create new column of traces for each FA
  FA_trace <- FA_toptrace %>%
    pivot_wider(-int, names_from = rank, values_from = m.z)
  
  # remove decimal
  trace_nodec <- as.data.frame(lapply(FA_trace, gsub, pattern = "\\..*", replacement = ""))
  # unite all traces into single column
  trace_nodec <- trace_nodec %>%
    unite(traces, trace_1:trace_20, sep = ",")
  # unite these traces with original .peak export
  newpeakexport <- full_join(peakexport, trace_nodec, by = "Name")
  
  # Join data frames by Name
  joinedpeaks <- full_join(newpeakexport, peaklist, by = "Name") 
  peakdetector <- joinedpeaks %>%
    select(Start.RT, Stop.RT, traces)
  # add optimization statement
  peakdetector$optimize <- "true"
  # add valley valley option
  peakdetector <- peakdetector %>%
    add_column(modus = "VV", .after = "Stop.RT")
  # add refernce column 
  peakdetector$ref <- "Reference" 
  
  write_delim(peakdetector, path = "~/Desktop/testdelim.txt", delim = "|", col_names = FALSE)
  
  
  # PEAK IDENTIFIER SCRIPT
  
  peakIDer <- joinedpeaks %>%
    select(Start.RT, Stop.RT, Name, traces)
  peakIDer <- peakIDer %>%
    add_column(CAS = "na", .after = "Name")
  peakIDer <- peakIDer %>%
    add_column(comment = "na", .after = "CAS")
  peakIDer <- peakIDer %>%
    add_column(contributor = "Ross Whippo", .after = "comment")
  peakIDer <- peakIDer %>%
    add_column(referenceId = "na", .after = "contributor")
  
  write_delim(peakIDer, path = "~/Desktop/peakIDer.txt", delim = "|", col_names = FALSE)
  

############### SUBSECTION HERE

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####

# RUN WITH TOP 5 NOT 20

# create list of top 10 m.z for each FA
FA_toptrace <- mz_int %>%
  arrange(desc(int)) %>%
  group_by(Name) %>%
  slice(1:5)

# list of names for columns 
FA_toptrace$trace <- "trace"
num <- rep(1:5, n_distinct(FA_toptrace$Name))
FA_toptrace$num <- num
FA_toptrace <- FA_toptrace %>%
  unite(rank, trace, num)
# create new column of traces for each FA
FA_trace <- FA_toptrace %>%
  pivot_wider(-int, names_from = rank, values_from = m.z)

# remove decimal
trace_nodec <- as.data.frame(lapply(FA_trace, gsub, pattern = "\\..*", replacement = ""))
# unite all traces into single column
trace_nodec <- trace_nodec %>%
  unite(traces, trace_1:trace_5, sep = ",")
# unite these traces with original .peak export
newpeakexport <- full_join(peakexport, trace_nodec, by = "Name")

# Join data frames by Name
joinedpeaks <- full_join(newpeakexport, peaklist, by = "Name") 
peakdetector <- joinedpeaks %>%
  select(Start.RT, Stop.RT, traces)
# add optimization statement
peakdetector$optimize <- "true"
# add valley valley option
peakdetector <- peakdetector %>%
  add_column(modus = "VV", .after = "Stop.RT")
# add refernce column 
peakdetector$ref <- "Reference" 
# peakdetector$traces <- ""

write_delim(peakdetector, path = "~/Desktop/testdelim.txt", delim = "|", col_names = FALSE)


# PEAK IDENTIFIER SCRIPT

peakIDer <- joinedpeaks %>%
  select(Start.RT, Stop.RT, Name, traces)
peakIDer <- peakIDer %>%
  add_column(CAS = "na", .after = "Name")
peakIDer <- peakIDer %>%
  add_column(comment = "na", .after = "CAS")
peakIDer <- peakIDer %>%
  add_column(contributor = "Ross Whippo", .after = "comment")
peakIDer <- peakIDer %>%
  add_column(referenceId = "na", .after = "contributor")
# peakIDer$traces <- ""

write_delim(peakIDer, path = "~/Desktop/peakIDer.txt", delim = "|", col_names = FALSE)
