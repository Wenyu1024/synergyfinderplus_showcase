get_triocomb_data <- function(id){
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
    mutate(id= 1:length(idx)) %>% 
    mutate(X3= str_detect(X1, "\\("))
  
  tmp4 <- tmp3 %>% 
    filter(X1== "DMSOUnknown") %>% 
    mutate(drug3= "DMSO") %>% 
    mutate(conc3= 0) %>% 
    select(id, drug3, conc3)
  
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
    select(id, drug3, conc3) %>% 
    mutate_at(.vars = 3,.funs = as.numeric)
  
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
    select(id, drug3, conc3)
  
  id <- rep(1:10, each = nrow(response)/10)
  complete <- bind_rows(tmp4, tmp5, tmp6) %>% 
    mutate(conc_unit3= "uM") %>% 
    arrange(id) %>% 
    bind_cols(calc %>% select(BlockId,	Replicate)) %>% 
    inner_join(meta,by = "BlockId") %>% 
    select(-id) %>% 
    distinct() %>% 
    separate_rows(RowConcs,sep = ",") %>% 
    separate_rows(ColConcs,sep= ",") %>% 
    rename( drug1= RowName) %>% 
    rename(drug2=ColName) %>% 
    rename(conc1= RowConcs) %>%  
    rename(conc2= ColConcs) %>% 
    rename(conc_unit1= RowConcUnit) %>% 
    rename(conc_unit2= ColConcUnit) %>% 
    arrange(conc1,decreasing=T) %>% 
    mutate(Row= id) %>% 
    arrange(conc2,decreasing=T) %>% 
    mutate(Col= id) %>% 
    left_join(response %>% rename(response = Value)) %>%
    mutate(conc1= as.numeric(conc1)) %>% 
    mutate(conc2= as.numeric(conc2)) %>%  
    select(BlockId,Replicate,drug1,drug2,drug3,conc1,conc2,conc3,response,conc_unit1,conc_unit2,conc_unit3,  RowSid, ColSid)
  
  
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
  
  return(complete_enhanced)
}