get_triocomb_data <- function(id, output_type= "DMSO_collapsed_data"){
  dir <-"C:/Users/wenyu/Documents/work group Jing/Drugcomb/synergyfinder/synergyfinderplus_showcase" 
  url.base <- "https://tripod.nih.gov/matrix-client/rest/matrix/export/%s?default"
  
  library(rvest)
  #load data
  url <- sprintf(url.base, id)
  file.name <- paste0(as.character(id), ".zip")
  if (! (file.name %in% list.files(paste0(dir, "/data/raw_data/"))  ) ){  
    download.file(url = url, destfile = paste0(dir, "/data/raw_data/", file.name), mode = "wb")}
  
  unzip(zipfile = paste0(dir, "/data/raw_data/", file.name), 
        exdir = paste0(dir, "/data/unzip/", id))
  calc <- read.csv(paste0("./data/unzip/", id, "/calc.csv"), stringsAsFactors = FALSE)
  meta <- read.csv(paste0("./data/unzip/", id, "/metadata.csv"), stringsAsFactors = FALSE)
  response <- read_csv(file = paste0("./data/unzip/", id, "/responses.csv"))
  response <- response %>% 
    select(-Layer) %>% 
    mutate(Replicate= str_split(Replicate,pattern = "default",n = 2,simplify = T)[,1]) %>% 
    mutate(Replicate= as.integer(Replicate))
  
  library(rvest)
  url_thirddrug= paste0("https://matrix.ncats.nih.gov/matrix-client/rest/matrix/blocks/", id, "/table", collapse = "")
  thirddrug_info <- read_html(url_thirddrug) 
  
  tmp <- thirddrug_info %>% 
    html_nodes(".idcell") %>%
    html_text()
  idx <- seq.int(from = 3,to = length(tmp),by = 3)
  tmp1 <- tmp[idx]
  
  tmp2 <- str_split(string = as.matrix(tmp1),pattern = "uM", n=2,simplify =T)
  
  tmp3 <- as_tibble(data.frame(tmp2)) %>% 
    mutate(serial= 1:length(idx)) %>% 
    mutate(X3= str_detect(X1, "\\("))

  # DMSO data points  
  tmp4 <- tmp3 %>% 
    filter(X1== "DMSOUnknown") %>% 
    mutate(drug3= "DMSO") %>% 
    mutate(conc3= 0) %>% 
    select(serial, drug3, conc3)
  # deall with drugs whose name does not have a space in it
  tmp5 <- tmp3 %>% filter(!X3) %>% 
    filter(X1 != "DMSOUnknown") %>%
    mutate(X4= str_split(X1, pattern = " ")) %>%
    mutate(drug3= map_chr(.x = X4,
                          .f = function(x) {
                            x1= unlist(x)
                            tmp <-str_c(x1[-length(x1):-(length(x1)-1)],
                                        collapse = " ")
                            return(tmp)
                          } )) %>%
    mutate(conc3= map_chr(.x = X4,
                          .f = function(x) {
                            x1= unlist(x)
                            tmp <- x1[ length(x1)-1]
                            return(tmp)
                          } )) %>%
    select(serial, drug3, conc3) %>% 
    mutate_at(.vars = 3,.funs = as.numeric)

  #deal with those drug who has a spaec in their names  
  tmp6 <- tmp3 %>% 
    filter(X3) %>% 
    mutate(X4= str_split(X1, pattern = "\\) ")) %>% 
    mutate(drug3= map_chr(.x = X4,.f = function(x) {
      tmp <- x[[1]]
      tmp1 <- paste0(tmp,  ")")
      return(tmp1)
    } )) %>% 
    mutate(conc3= map_dbl(.x = X4,.f = function(x) {
      tmp <- as.numeric(x[[2]])
      return(tmp)
    } ))  %>% 
    select(serial, drug3, conc3)
  
  # combine all three parts together
  webtable_drug3 <- bind_rows(tmp4, tmp5, tmp6) %>% 
    mutate(conc_unit3= "uM") %>% 
    arrange(serial) 
  
  # webtable can be directly combined with calc because I mannually checked that the MedianExcess column is the same  (also in same order)
  # between two table. calc table proivde the block an replicate id and can then be used to join with meta table
  # now I have this enhanced metat which also have additional columns for drug3
  meta_enhanced <- webtable_drug3 %>% 
    bind_cols(calc %>% select(BlockId,	Replicate)) %>% 
    inner_join(meta,by = "BlockId") %>% 
    select(-serial) %>% 
    distinct() 

  #expand the enhanced meta table to release the folded conc1 and conc2 
  meta_expanded <- meta_enhanced %>% 
    rename( drug1= RowName) %>% 
    rename(drug2=ColName) %>% 
    rename(conc1= RowConcs) %>%  
    rename(conc2= ColConcs) %>% 
    rename(conc_unit1= RowConcUnit) %>% 
    rename(conc_unit2= ColConcUnit) %>% 
    separate_rows(conc1,sep = ",") %>% 
    separate_rows(conc2,sep= ",") %>% 
    mutate(conc1= as.numeric(conc1)) %>% 
    mutate(conc2= as.numeric(conc2)) %>%  
    select(BlockId,Replicate,drug1,drug2,drug3,conc1,conc2,conc3,conc_unit1,conc_unit2,conc_unit3,  RowSid, ColSid)
  
  # now to join this expanded meta table with the response table 
  # I have to generate so called row and column id 
  # as the response table does not use the drug names and concentration to indicate data point.
  # Instead response table use block to indicate drug pair and row/column id to indicate concentration
  id <- rep(1:10, each = nrow(response)/10)
  
  complete <- meta_expanded %>% 
    arrange(desc(conc1)) %>% 
    mutate(Row= id) %>% 
    arrange(desc(conc2)) %>% 
    mutate(Col= id) %>% 
    inner_join(response %>% rename(response = Value), by = c("BlockId", "Col", "Row", "Replicate")) %>%
    select(- Col, - Row) %>% 
    arrange(BlockId, Replicate, desc(conc3), desc(conc2), desc(conc1)) 

  # now deal with the DMSO issue, basically use response of rows where third drug is DMSO as zero concentration response for all the other drug 3
  
  DMSO <- complete %>% filter(drug3 == "DMSO")
  
  other_drug3 <- complete %>% filter(drug3 != "DMSO") %>% pull(drug3) %>% unique()
  
  complete_enhanced <- DMSO %>% 
    mutate(drug3= str_c(other_drug3,collapse = ",")) %>% 
    separate_rows(drug3,sep = ",") %>% 
    bind_rows(complete %>% filter(drug3 != "DMSO")) %>% 
    distinct() %>%
    select(- BlockId)
  
  blockid= group_indices(.data = complete_enhanced, drug1,drug2,drug3)
  complete_enhanced <- complete_enhanced %>%  mutate(block_id= blockid)
  
  if (output_type== "DMSO_collapsed_data"){return(complete_enhanced)}
  if (output_type== "raw_data"){return(complete)}
}


get_synergy_summary <- 
  function(id){
    dir <-"C:/Users/wenyu/Documents/work group Jing/Drugcomb/synergyfinder/synergyfinderplus_showcase" 
    trio_set_file_name <-  paste0(dir,"/data/trio_data/", id, ".csv")
    
    # if (paste0(id, ".csv") %in% list.files(paste0(dir,"/data/trio_data/"))){
    #   data1 = read_csv(trio_set_file_name)
    # } else {
      data1 <- get_triocomb_data(id) 
      write_csv(data1 , file = trio_set_file_name )
    # }
    
    res_response = ReshapeData(data1, data_type = "viability")
    res_hsa <- CalculateSynergy(res_response , method = c("HSA"))
    # res_Bliss <- CalculateSynergy(res_response , method = c("Bliss"))
    hsa_summary <- res_hsa$synergy_scores %>% 
      select(block_id, HSA_synergy) %>% 
      group_by(block_id) %>% 
      summarize(median= median(HSA_synergy) ,
                mean= mean(HSA_synergy), 
                max= max(HSA_synergy),
                min= min(HSA_synergy)
      ) %>% 
      ungroup()
    # bliss_summary <- res_Bliss$synergy_scores %>% 
    #   select(block_id, Bliss_synergy) %>% 
    #   group_by(block_id) %>% 
    #   summarize(median= median(Bliss_synergy) ,
    #             mean= mean(Bliss_synergy), 
    #             max= max(Bliss_synergy)
    #   )
    # res <- c(hsa=list(hsa_summary), bliss=list(bliss_summary))
    return(hsa_summary)
  } 