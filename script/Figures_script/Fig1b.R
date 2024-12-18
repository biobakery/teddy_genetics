setwd("~/ddong/TEDDY_project/")
library(tidyverse)
library(data.table)

##antibiotic 
pencillin_amox_ind <- c("MED00002", "MED00003", "MED00108", "MED01081", "MED00679", 
                        "MED00062", "MED00106", "MED00107", "MED00113", "MED00176", 
                        "MED00225", "MED00286", "MED00605", "MED00712")
pencillin_ind <- c("MED00679", "MED00062", "MED00106","MED00107", "MED00113", 
                   "MED00176", "MED00225", "MED00286", "MED00605", "MED00712")
macrolide_ind <- c("MED00005", "MED00041", "MED00114", "MED00154", "MED00481", 
                   "MED00640", "MED01162")
cephalosporins_ind <- c("MED00060", "MED00072", "MED00095", "MED00098", "MED00161", 
                        "MED00194", "MED00220", "MED00226", "MED00241", "MED00311", 
                        "MED00317", "MED00330", "MED00404", "MED00532", "MED00619", 
                        "MED00717", "MED01091")
amoxicillin_ind <- c("MED00002", "MED00003", "MED00108", "MED01081")
carbacephem_ind <- c("MED00538", "MED00544")
other_ind <- c("MED00089", "MED000093", "MED00129", "MED00192", "MED00290", 
               "MED00296", "MED00318", "MED00350", "MED00435", "MED00453", 
               "MED00561", "MED00598", "MED00621", "MED00638", "MED00667", 
               "MED00677", "MED00895", "MED01052", "MED01067", "MED01157", 
               "MED01179", "MED00628", "MED00674", "MED00538", "MED00544", 
               "MED00120",	"MED00248", "MED00145", "MED00368", "MED00437",
               "MED00591", "MED01047", "MED01057", "MED00238", "MED00199",
               "MED00190", "MED00014", "MED00135", "MED00242", "MED00362", 
               "MED00479", "MED00790", "MED01059", "MED00165", "MED00324", 
               "MED00405", "MED00566", "MED01151", "MED00081", "MED01193", 
               "MED01115", "MED00138", "MED00271",	"MED00845", "MED01089", 
               "MED00080", "MED00508", "MED00743")
antibiotics <- c(pencillin_amox_ind, pencillin_ind, macrolide_ind, cephalosporins_ind, 
                 amoxicillin_ind, carbacephem_ind)


# subject level metadata
subjects <- fread("./data/TEDDY_metadata/metadata_subject.txt",header = T)
metadata <- fread("./data/TEDDY_metadata/2020-06-23/mp201_subject_level_masked.csv") %>%
  rename(subject_id=1) %>%
  select(1:19, 128:142) %>%
  mutate(country = recode(country, `1`="USA", `2`="Finland", `3`="Germany", `4`="Sweden")) %>%
  as.data.frame()

samples <- fread("./data/data_derived/metadata_for_samples.tsv", header=T) %>%
  column_to_rownames("sample_mask_id") %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  rename(subject_id=mask_id, mgx_age=subject_age_days)
samples <- samples %>% left_join(subjects)
metadata <- merge(samples, metadata, by="subject_id") %>%
  mutate(current_t1d_diag = case_when(
    mgx_age >= t1d_age ~ "Y" ,
    mgx_age < t1d_age ~ "N")) %>%
  select(-t1d_age) %>%
  mutate(current_t1d_diag = tidyr::replace_na(current_t1d_diag, "N")) %>%
  mutate(current_persist_conf_ab = case_when(
    mgx_age >= persist_conf_age ~ "Y",
    mgx_age < persist_conf_age ~ "N")) %>%
  select(-persist_conf_age) %>%
  mutate(current_persist_conf_ab = tidyr::replace_na(current_persist_conf_ab, "N")) %>%
  mutate(current_brst_fed = case_when(
    mgx_age < time_to_brstfed_stop ~ "Y",
    mgx_age >= time_to_brstfed_stop ~ "N"
  )) %>%
  select(-time_to_brstfed_stop) %>%
  mutate(current_brst_fed = tidyr::replace_na(current_brst_fed, "N")) %>%
  mutate(ever_brstfed = recode(ever_brstfed, `0`="N", `1`="Y")) %>%
  mutate(t1d = recode(t1d, `0`="N", `1`="Y")) %>%
  mutate(persist_conf_ab = recode(persist_conf_ab, `0`="N", `1`="Y")) %>%
  rename(abx_consumed=AntibioticsConsumed) %>%
  mutate(abx_consumed = tidyr::replace_na(abx_consumed, 0)) %>%
  mutate(abx_consumed = recode(abx_consumed, `0`="N", `1`="Y")) %>%
  mutate(delivery = recode(delivery, `1`="vaginal", `2`="c_section", `3`="vaginal", `4`="c_section"))

# childhood medications

meds <- read.csv("./data/data_derived/metadata_cleaned/mp201_child_medications.tsv", sep="\t", header=T) %>%
  rename(subject_id=1, abx_age=2, medcode=4, duration=6, ongoing=7) %>%
  select(1:2, 4, 6:7) %>%
  mutate(subject_id = gsub("b\\'|'", "", subject_id)) %>%
  mutate(medcode = gsub("b\\'|'", "", medcode)) %>%
  mutate(class = case_when(
    medcode %in% pencillin_amox_ind ~ "pencillin_amox_ind",
    medcode %in% pencillin_ind ~ "pencillin_ind",
    medcode %in% macrolide_ind ~ "macrolide_ind",
    medcode %in% cephalosporins_ind ~ "cephalosporins_ind",
    medcode %in% amoxicillin_ind ~ "amoxicillin_ind",
    medcode %in% carbacephem_ind ~ "carbacephem_ind"
  )) %>%
  mutate(antibiotics = case_when(medcode %in% antibiotics ~ "antibiotics")) %>%
  drop_na() #%>%
#filter(abx_age <= 1095)
head(meds)

#

md <- metadata %>%
  rename(mgx_age=3) %>%
  merge(., meds, by="subject_id") %>%
  mutate(days = mgx_age - abx_age) %>%
  filter(days >= 0 & days <= 30) %>%
  mutate(abx_last_30_days = "Y") %>%
  select(sample_id, abx_last_30_days) %>% unique()

metadata_updated <- merge(metadata, md, by="sample_id", all.x=TRUE) %>%
  mutate(abx_last_30_days = tidyr::replace_na("N"))# %>%
#filter(mgx_age <= 1095)



######solid food
metadata_subject_all<- fread("./data/TEDDY_metadata/2020-06-23/mp201_subject_level_masked.csv")
solidfood <- c("oat", "rye", "buckwheat_and_millet", "fruits_and_berries", "Fruit_or_berries", 
               "potatoes", "sweet_potatoes_or_yams", "carrots", "spinach", "beets", "peas_green_beans",
               "Turnip_parsnip_artichoke_rutabag", "cabbages", "squash_pumpkin",
               "tomato_or_tomato_sauce", "corn", "other_vegetable", "rice", "wheat", "barley",
               "pork_beef", "poultry", "other_kinds_of_meat", "apple_sauce_or_apple_juice", 
               "sausage_hot_dogs", "fish_or_other_seafood", "egg", "Commercial_baby_food_containing",
               "Soy_milk_and_other_soy_products", "rice_milk", "goat_milk","gluten", "cereals", 
               "nongluten_cereals",    "noncorn_cereals",  "nongluten_noncorn_cereals", 
               "rice_ricemilk",  "root_vegetables", "total_roots")

metadata_solid <- metadata_subject_all  %>% select( one_of(solidfood))
metadata_solid <- apply(metadata_solid, 1, function(x){ min(x)})
metadata_subject_all$age_startsolid <- metadata_solid
metadata_solid <- metadata_subject_all %>% select(MP201_MaskID,age_startsolid,time_to_brstfed_stop)


##############################
metadata_subject <- metadata %>% select(subject_id,delivery, abx_consumed, ever_brstfed,Sex,delivery,brst_fed,Country,Clinical_Center,Probiotic) %>% unique()
metadata_subject <- metadata_subject %>% left_join(metadata_solid,by = c("subject_id" = 'MP201_MaskID'))
library(table1)
table(metadata_subject$brst_fed)
metadata_subject_NA <- metadata_subject[apply(metadata_subject, 1, function(x) any(is.na(x))), ]
metadata_subject$time_to_brstfed_stop[is.na(metadata_subject$time_to_brstfed_stop)] <- mean(metadata_subject$time_to_brstfed_stop,na.rm= T)
table1(~ Sex + delivery + age_startsolid + ever_brstfed +time_to_brstfed_stop+abx_consumed + Probiotic+Clinical_Center| Country, data=metadata_subject)
save(metadata_subject,metadata,file = "./output/preprpreprocess/metadata_all.Rdata")


load("~/ddong/TEDDY_project//output/preprpreprocess/metadata_all.Rdata")
pdf('~/ddong/TEDDY_project/figures/Fig1b.pdf',width = 6,height = 4.5)
ggplot(metadata,aes(x= mgx_age, fill = Country ))+geom_histogram()+theme_cowplot()+
  ylab("Number of samples") +xlab("Time (Days)")
dev.off()

