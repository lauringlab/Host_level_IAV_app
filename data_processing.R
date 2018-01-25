require(magrittr)
require(tidyverse)
if(!("package:tidyverse" %in% search())){
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(readr)
  require(rlang)
}
require(HIVEr)

set.seed(42) # Set seed so randomization is reproducible

# ------------------------- Read in data ------------------------------------------------
meta<-read_csv("../Host_level_IAV_evolution/data/reference/all_meta.sequence_success.csv")

intra<-read_csv("../Host_level_IAV_evolution/data/processed/secondary/Intrahost_all.csv")
intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra
intra<-subset(intra,freq1<0.5) # Just to make sure


trans<-read_csv("../Host_level_IAV_evolution/data/processed/secondary/trans_freq.csv")
trans<-trans %>% filter(freq1<0.5) %>%
  mutate(Endpoint="Persistent") %>%
  mutate(Endpoint = if_else(freq1==0,"Arisen",Endpoint)) %>%
  mutate(Endpoint = if_else(freq2==0,"Lost",Endpoint))
trans$Endpoint<-factor(trans$Endpoint,levels = c("Persistent","Arisen","Lost"),ordered = T)

trans <- trans %>% mutate(
  within_host_time = abs(transmission-collect1)+abs(collect2-transmission)
  )

# Now I will work to get longitudinal pairs for people in the transmission subset. 
# The first step is to subset the transmission data into just the meta data for each pair
# and rename the columns according to whom the data comes from. The data starts as 1 row
# for each mutation so the distinct calapses this.

trans_long<- trans %>% select(HOUSE_ID,Donor_ENROLLID=ENROLLID1,
                              Recipient_ENROLLID=ENROLLID2,
                              Donor_onset=onset1,Recipient_onset=onset2,
                              transmission) %>% distinct()

# Now we will get the SPECID for the donor by joining with the meta data
# we are only interested in samples that qualified for snv identification.
# NB : left join includes all combinations so if there are two SPECID in meta, 
# which is present as two rows, then both are kept, again as two rows. Key columns 
# added and the SPECID are the spread to 1 row with _home and _clinic columns.
#  Each time we rename the added SPECID from the meta data according to
#  whose ENROLLID we are using to join.
trans_long<-trans_long %>%left_join(filter(meta,snv_qualified==T),
              by=c("Donor_ENROLLID"="ENROLLID"))%>% 
            select(HOUSE_ID=HOUSE_ID.x,Donor_ENROLLID,
                                      Recipient_ENROLLID,
                                      Donor_onset,Recipient_onset,
                                      transmission,Donor_SPECID=SPECID) %>%
        mutate(Donor_sample=if_else(condition = grepl("HS",Donor_SPECID),
                              true = "Donor_home",
                              false = "Donor_clinic")) %>%
        spread(key = Donor_sample,value = Donor_SPECID)


# Now we do the same for the recipient.
trans_long<- left_join(trans_long,filter(meta,snv_qualified==T),
                            by=c("Recipient_ENROLLID"="ENROLLID")) %>%
            select(HOUSE_ID=HOUSE_ID.x,Donor_ENROLLID,
                                      Recipient_ENROLLID,
                                      Donor_onset,Recipient_onset,
                                      transmission,Donor_home,Donor_clinic,
                                      Recipient_SPECID=SPECID) %>% 
          mutate(Recipient_sample=if_else(condition = grepl("HS",Recipient_SPECID),
                           true = "Recipient_home",
                           false = "Recipient_clinic")) %>%
          spread(key = Recipient_sample,value = Recipient_SPECID)

# Now we will add the collection dates, by joining with meta and renaming the columns
# 

trans_long <- left_join(trans_long,select(filter(meta,snv_qualified==T),SPECID,collect), # Donor home
                        by=c("Donor_home" = "SPECID")) %>%
  rename("Donor_home_collect" = "collect") %>%
  left_join(.,select(filter(meta,snv_qualified==T),SPECID,collect), # Donor clinic
            by=c("Donor_clinic" = "SPECID")) %>%
  rename("Donor_clinic_collect" = "collect") %>%
  left_join(.,select(filter(meta,snv_qualified==T),SPECID,collect), # Recipient home
            by=c("Recipient_home" = "SPECID")) %>%
  rename("Recipient_home_collect" = "collect")%>%
  left_join(.,select(filter(meta,snv_qualified==T),SPECID,collect), # Recipient home
            by=c("Recipient_clinic" = "SPECID")) %>%
  rename("Recipient_clinic_collect" = "collect")


#---------------- Adding iSNV comparisons -------------------------

qual<-read_csv("../Host_level_IAV_evolution/data/processed/secondary/qual.snv.csv",
                col_types = list(
                ENROLLID= col_character(),
                SPECID = col_character(),
                LAURING_ID = col_character(),
                Id = col_character()
                )) # read in quality variant calls from all 

# Donor comparisions
DD<- trans_long  %>% filter(.,!(is.na(Donor_home)) & !(is.na(Donor_clinic))) %>%
  rowwise() %>%
    do(get_freqs(c(.$Donor_home,.$Donor_clinic),qual)) %>% 
    rename(Donor_home=SPECID1,Donor_clinic=SPECID2,
    Donor_home_freq = freq1, Donor_clinic_freq=freq2) %>%
    left_join(.,trans_long)

RR<-trans_long  %>% filter(.,!(is.na(Recipient_home)) & !(is.na(Recipient_clinic))) %>% rowwise() %>%
      do(get_freqs(c(.$Recipient_home,.$Recipient_clinic),qual)) %>% 
      rename(Recipient_home=SPECID1,Recipient_clinic=SPECID2,
          Recipient_home_freq = freq1, Recipient_clinic_freq=freq2) %>%
      left_join(.,trans_long)

DhRc <-trans_long  %>% filter(.,!(is.na(Donor_home)) & !(is.na(Recipient_clinic))) %>% rowwise() %>%
          do(get_freqs(c(.$Donor_home,.$Recipient_clinic),qual)) %>% 
          rename(Donor_home=SPECID1,Recipient_clinic=SPECID2,
          Donor_home_freq = freq1, Recipient_clinic_freq=freq2) %>%
          left_join(.,trans_long)

DhRh <- trans_long  %>% filter(.,!(is.na(Donor_home)) & !(is.na(Recipient_home))) %>% rowwise() %>%
          do(get_freqs(c(.$Donor_home,.$Recipient_home),qual)) %>% 
          rename(Donor_home=SPECID1,Recipient_home=SPECID2,
          Donor_home_freq = freq1, Recipient_home_freq=freq2) %>%
          left_join(.,trans_long)

RhDc <- trans_long  %>% filter(.,!(is.na(Recipient_home)) & !(is.na(Donor_clinic))) %>% rowwise() %>%
        do(get_freqs(c(.$Recipient_home,.$Donor_clinic),qual)) %>% 
        rename(Recipient_home=SPECID1,Donor_clinic=SPECID2,
        Recipient_home_freq = freq1, Donor_clinic_freq=freq2) %>%
        left_join(.,trans_long)

RcDc <- trans_long  %>% filter(.,!(is.na(Recipient_clinic)) & !(is.na(Donor_clinic))) %>% rowwise() %>%
      do(get_freqs(c(.$Recipient_clinic,.$Donor_clinic),qual)) %>% 
      rename(Recipient_clinic=SPECID1,Donor_clinic=SPECID2,
      Recipient_clinic_freq = freq1, Donor_clinic_freq=freq2) %>%
      left_join(.,trans_long)


# These joins are not working - there can't be NAs in both donor columns for freq There are mutations with 
# NA in some row and then the mutations are repeated with the real frequency filled in. 
# I think the solution is in asigning the by variables and dropping some.
full_join(DD,DhRc)->isnv_tl

full_join(isnv_tl,RcDc)->isnv_tl

full_join(isnv_tl,DhRh)->isnv_tl

full_join(isnv_tl,RhDc)->isnv_tl

full_join(isnv_tl,RR)->isnv_tl


# This seems to fix it
isnv_tl<- isnv_tl %>%
group_by(mutation,Donor_ENROLLID,Recipient_ENROLLID) %>%
  mutate(
    Donor_home_freq = 
      ifelse(length(which(!is.na(Donor_home_freq)))>0,
        unique(Donor_home_freq[which(!is.na(Donor_home_freq))]),
        unique(Donor_home_freq)),
  Donor_clinic_freq = 
      ifelse(length(which(!is.na(Donor_clinic_freq)))>0,
        unique(Donor_clinic_freq[which(!is.na(Donor_clinic_freq))]),
        unique(Donor_clinic_freq)),

  Recipient_home_freq = 
    ifelse(length(which(!is.na(Recipient_home_freq)))>0,
      unique(Recipient_home_freq[which(!is.na(Recipient_home_freq))]),
      unique(Recipient_home_freq)),
  Recipient_clinic_freq = 
    ifelse(length(which(!is.na(Recipient_clinic_freq)))>0,
      unique(Recipient_clinic_freq[which(!is.na(Recipient_clinic_freq))]),
      unique(Recipient_clinic_freq))) %>%
  distinct()


# Collection??
isnv_tl<- isnv_tl %>% rowwise() %>% mutate(pair = paste0(Donor_ENROLLID,"-",Recipient_ENROLLID)) 
# Think about frequency cut off.
write.csv(x = filter(isnv_tl,Donor_home_freq<0.5 | Donor_clinic_freq<0.5),"~/Documents/Analysis/Influenza_host_level_app/freq.wide.csv")

# Make long form
filter(isnv_tl,Donor_home_freq<0.5|Donor_clinic_freq<0.5) %>% gather(sample,freq,ends_with("freq")) %>% 
  mutate(sample = 
           gsub(pattern = "_freq",x = sample,replacement = "")) ->long_freq
filter(isnv_tl,Donor_home_freq<0.5|Donor_clinic_freq<0.5) %>% gather(sample,collection,ends_with("collect"))%>%
  mutate(sample = 
           gsub(pattern = "_collect",x = sample,replacement = ""))->long_collect
left_join(select(long_freq,-ends_with("collect")),select(long_collect,-ends_with("freq")))->isnv_long
isnv_long <- isnv_long %>%  mutate(day = collection-Donor_onset,
                                   sample_class=gsub(pattern = "(.*)_.*",x=sample,replacement = "\\1",perl = T),
                                   pair = paste0(Donor_ENROLLID,"-",Recipient_ENROLLID)) 
write.csv(x = isnv_long,"~/Documents/Analysis/Influenza_host_level_app/freq.long.csv")


