# NCATS dataset https://matrix.ncats.nih.gov/ "Malaria TACT" project all 3-agent combinations
# Download files
url.base <- "https://tripod.nih.gov/matrix-client/rest/matrix/export/%s?default"

# change following line with the adday IDs on NCATA Matrix webpage
id <- 10023
url <- sprintf(url.base, id)
file.name <- paste0(as.character(id), ".zip")
download.file(url = url, destfile = paste0("./data/raw_data/", file.name), 
              mode = "wb")
# unzip(paste0("./data/raw_data/", file.name), exdir = paste0("./data/unzip/", id))

meta <- read.csv(paste0("./data/unzip/", id, "/metadata.csv"), stringsAsFactors = FALSE)
response <- read.csv(paste0("./data/unzip/", id, "/responses.csv"), stringsAsFactors = FALSE)

# Change following line with the concentrations for 3rd drug (available on webpage)
meta$Conc3 <- c(
  0.75, 0.375, 0.1875, 0.938, 0.469, 0.0234, 0.0117, 0.0059,
  0.0029, 0.0015, 0.0007, 0
)
# Change following line with the 3rd drug's name and lengt 
meta$Drug3 <- rep("Piperaquine", length(meta$Conc3))

data <- NULL
%>% 
for (i in 1:12) {
  Conc1 <- as.numeric(unlist(strsplit(meta$RowConcs[i], ",")))
  Conc2 <- as.numeric(unlist(strsplit(meta$ColConcs[i], ",")))
  tmp <- response[which(response$BlockId == i), ]
  # Change "1" in the following line for different assays/blocks
  df <- data.frame(block_id = rep(1, nrow(tmp)), stringsAsFactors = FALSE)
  df$drug1 <- rep(meta$RowName[i], nrow(df))
  df$drug2 <- rep(meta$ColName[i], nrow(df))
  df$drug3 <- rep(meta$Drug3[i], nrow(df))
  df$donc1 <- rep(Conc1, each = 10)
  df$donc2 <- rep(Conc2, times = 10)
  df$donc3 <- rep(meta$Conc3[i], nrow(df))
  df$response <- tmp$Value
  df$conc_unit1 <- rep(meta$RowConcUnit[i], nrow(df))
  df$conc_unit2 <- rep(meta$ColConcUnit[i], nrow(df))
  df$conc_unit3 <- rep("uM", nrow(df))
  data <- rbind.data.frame(data, df)
}

# Run codes above on all the assays and merge all tables together.
