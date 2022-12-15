library(pracma)
library(readr)
library(dplyr)
library("tidyr")
library(data.table)
library(igraph)

setwd("~/Onedrive - Emory University/Fall_SNA/Social Network Analysis - Final Project")
raw <- fread("WS_CBS_PUB_csv_col.csv")

#-----------clean the crossholdings----------
# reporting type: domestic reporting 
# remove non-country
# remain maturity: including all years (including lower than 1 year data)
# cbs reporting basis: immediate counterparty basis
# balance sheet position: total claims
# Type of instruments: All
# remove panama because of low gdp value
raw <- raw[CBS_BANK_TYPE == "4R" & L_REP_CTY != "5A" & REM_MATURITY == "A" & CBS_BASIS == "F" & L_POSITION == "C" & L_REP_CTY != "PA"]
raw <- na.omit(raw, cols=c("L_REP_CTY", "L_CP_COUNTRY", "2021-Q4"))

setkeyv(raw, c("L_REP_CTY", "L_CP_COUNTRY"))

dt <- raw[, c("L_REP_CTY", "L_CP_COUNTRY",
              "Reporting country", "Counterparty country", 
              "2021-Q4")]

setnames(dt, c("from_country_code","to_country_code", "from_country", "to_country", "weight") )
# remove `countries` that is not a country
dt <- dt[is.na(as.numeric(substr(to_country_code, 1, 1))) == 1,]

write.csv(dt, file = "crossholding dataset.csv")
# save csv to draw graphs

#--------clean the gdps--------
gdp <- fread("gdp_imf_list.csv", header=TRUE)
# get gdp minimum of every year
setkeyv(gdp, c("country_code"))

gdp_min <- apply(gdp[,3:ncol(gdp)],2,min)
# normalize gdp per year
gdp_nml <- sweep(gdp[,3:ncol(gdp)], 2, gdp_min, FUN = '/')
gdp_nml <- cbind(gdp[,1:2], gdp_nml)

gdp_22 = gdp_nml[, c(1,2,10)]
gdp_19 = gdp_nml[, c(1,2,7)]

#------build adjacency matrix--------

debt <- graph_from_data_frame(dt, directed = TRUE)
debt <- set_edge_attr(debt, "weight", value= dt$weight)

debt_mt <- as_adjacency_matrix(debt,sparse=FALSE,attr="weight")
debt_mt <- debt_mt[rowSums(debt_mt^2)>0,]
country_list <- rownames(debt_mt)
debt_mt <- debt_mt[, colnames(debt_mt) %in% country_list]

n <- length(country_list)

# sum over every column
debt_mt_sum <- t(data.matrix(colSums(debt_mt)))

# Here we estimate the ratio of total debt held outside the issuing country by 1/3.
c <- .33

# get C bar (private holdings in a country) and save as diag matrix
C_row <- as.vector(matrix(rep(1,n), ncol = n)/debt_mt_sum)
C_diag <- diag(C_row)

# Now generate the matrix of shares

# normalizes columns to sum to 1
# imposes that only c of a country is held outside
C <- c * debt_mt %*% C_diag

# the remainder is held by private shareholders, i.e. the citizens
Chat <- (1-c)*diag(n)

# get dependency matrix A
A <- Chat%*%solve(diag(n)-C)
A_df <- data.frame(A)
# save primitive holding to csv
write_csv(A_df, "primitive_holding.csv")


# p as predicted gdp
p <- as.matrix(gdp_22[,3])
pcurrent <- as.vector(p)
beta <- 0.5

# function for wave trend by theta
contagion <- function(theta){
  # set wave number
  wave <<- 1
  first_failed_num <- 0
  
  # define theta and threshold
  v_threshold <- theta*(A %*% as.matrix(gdp_19[,3]))
  
  # set failure indicator and failed country saver
  all_failure_indicator <- matrix(rep(0),25,1)
  new_failed_countries <- c(1)
  all_failure_list <- list()
  
  # loop for waves
  while(!isempty(new_failed_countries)){
    
    a <- as.matrix((as.numeric(A %*% pcurrent < v_threshold)))
    b <- as.matrix(as.numeric(all_failure_indicator == 0))
    
    new_failure_indicator <-  a * b
    
    # find pairwise maximum
    all_failure_indicator <- pmax(all_failure_indicator, new_failure_indicator)
    
    # cut the failed country's gpd by half of threshold
    pcurrent <- pcurrent - new_failure_indicator * (v_threshold * beta)
    
    # find new_failed_countries 
    new_failed_countries <- which(new_failure_indicator == 1)
    # get failed country names
    new_failed_names = matrix(country_list[new_failed_countries], nrow = 1)
    
    if(wave == 1) {
      first_failed_num <<- length(new_failed_countries)
    }
    
    if (!isempty(new_failed_countries)){
      # cat('Wave ', wave, ' failures are ')
      # cat(new_failed_names)
      all_failure_list[[length(all_failure_list)+1]] <- list(new_failed_names)
      
    }
    else{
      break
    }
    # cat("\n","end wave", wave,"\n")
    wave <<- wave + 1
  }
  all_failed_num <<- length(which(all_failure_indicator==1))
  print(all_failure_list)
  return(all_failure_list)
}

theta = c()
wave_length = c()
first_fail = c()
all_failed = c()

i <- 1
# for loop of theta from .9 to 1.2
for (t in seq(0.9, 1.2, by = 0.01)){
  cat("start loop", i, "\n")
  contagion(t)
  theta <- append(theta, t)
  wave_length <- append(wave_length, wave-1)
  first_fail <- append(first_fail, first_failed_num)
  all_failed <- append(all_failed, all_failed_num)
  i <- i+1
}

wave.df<- data.frame(
  Theta = theta,
  Total_Wave = wave_length,
  First_Fail_Num = first_fail,
  All_Fail_Num = all_failed
)

write_csv(wave.df, "wave.csv")

theta = 0.96
simulist <- contagion(0.96)

# create full-connected graph
fc_graph <- graph_from_adjacency_matrix(debt_mt, mode = "directed", weighted=TRUE)
current_infected <- c()
layout.old = layout.fruchterman.reingold(fc_graph, niter = 1000)

plot_gif = function(simulist){
  m = 1
  while(m <= length(simulist)){
    layout(matrix(c(1, 2, 1, 3), byrow = TRUE), widths=c(3,1), heights=c(1, 1))
    current_infected <- c(current_infected, unlist(simulist[[m]]))
    V(fc_graph)$color = "white"
    V(fc_graph)$color[V(fc_graph)$name %in% current_infected] = "red"
    plot(fc_graph, edge.arrow.size=0.2, layout =layout.old)
    m = m + 1}
}

library(animation)
plot_gif(simulist)

saveGIF({
  ani.options(interval = 0.5, convert = shQuote("C:/Program Files/ImageMagick-6.9.9-Q16-HDRI/convert.exe"))
  # start the plot
  plot_gif(simulist)
}, ani.width = 1600, ani.height = 1000)


