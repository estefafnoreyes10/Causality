# #run once starts 
# # load packages, if not installed, install them
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if (!require(mgcv)) {install.packages("mgcv"); library(mgcv)}
if (!require(geometry)) {install.packages("geometry"); library(geometry)}
if (!require(furrr)) {install.packages("furrr"); library(furrr)}
if (!require(pracma)) {install.packages("pracma"); library(pracma)}
if (!require(combinat)) {install.packages("combinat"); library(combinat)}
if (!require(readxl)) {install.packages("readxl"); library(readxl)}
if (!requireNamespace("motif", quietly = TRUE)) {
  install.packages("motif")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("pcalg", quietly = TRUE)) install.packages("pcalg")
if (!require(tibble)) {
  install.packages("tibble")
  library(tibble)
}

install.packages("ggraph")
install.packages("graphlayouts")

BiocManager::install("RBGL")
library(dplyr)
library(readxl)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(RBGL)
library(Rgraphviz)
library(pcalg)
library(readr)
library(igraph)
library(motif)
library(ggraph)
library(graphlayouts)



# # Check current working directory
getwd()
setwd("/Users/estefano/Downloads/vaginal_ci/vaginal_ci")

# # Now the relative path 'data/vaginal.xlsx' should work
file_path <- "data/vaginal.xlsx"
file.exists(file_path)  # This should return TRUE if the file exists
data <- read_excel(file_path, skip = 4)
# # import the metadata about taxonomy information
taxainfo <- read_delim(file = "data/dict.txt", delim = " - ") %>% transmute(D4 = Family, D5 = Genus)
file_path_taxa<-"data/dict.txt"
data <- read_excel("data/vaginal.xlsx", skip = 4)
data <- data[, c(-4, -5, -6, -7, -8, -10)]
colnames(data)[c(1, 2, 3, 4)] <- c("sample_ID", "time_in_study", "sub_ID", "total_reads")
# # convert perceptage to absolute read counts
data <- data %>% mutate(across(.cols = 5:ncol(data), .fns = ~ . * total_reads / 100))

process_data_given_id_00 <- function(data,sub_id,taxainfo){
  #process data for one id and filter columns 
  data_time <- data %>%
    filter(sub_ID == sub_id) %>%
    select(c(2, 5:ncol(data)))
  
  #data filter by taxa 
  data_vg <- data_time %>%
    pivot_longer(cols = !time_in_study, names_to = "genus", values_to = "abundances") %>%
    pivot_wider(id_cols = genus, names_from = "time_in_study", names_prefix = "t_", values_from = "abundances") %>%
    select(-genus) %>%
    bind_cols(taxainfo, .) %>%
    filter(!grepl("Unclassified", .$D4, fixed = TRUE))
  
  # Aggregate data to the genus level
  
  #data filtered by genus
  data_genus <- data_vg %>%
    group_by(D4, D5) %>%
    summarize(across(c(everything()), sum, .names = "s_{.col}"), .groups = "keep") %>%
    unite(eff_taxo, D4, D5, sep = ";")
  
  
  #data filtered by family 
  data_family <- data_vg %>%
    group_by(D4) %>%
    summarize(across(c(-1), sum, .names = "s_{.col}"), .groups = "keep") %>%
    unite(eff_taxo, D4, sep = ";")
  
  # Calculate the total reads throughout the time period for each taxon
  presence <- data_family %>%
    transmute(total = rowSums(select(., -1)))
  
  # Filter out the taxa with too few reads
  group_candidates <- which(presence > 10^(-6) * max(presence))
  #species_names<-data_family['eff_taxo']
  #species_indices<-length(rownames(species_names))
  return (list(data_time=data_time,data_vg=data_vg,data_genus=data_genus,data_family=data_family,presence=presence,group_candidates=group_candidates))
}

process_data_all_ids_01 <- function(data, taxainfo, id_only) {
  ids_data_family_container <- list()
  for (sub_id in id_only) {
    id_processed_data <- process_data_given_id_00(data, sub_id, taxainfo)
    id_processed_data_family <- id_processed_data$data_family
    ids_data_family_container[[sub_id]] <- id_processed_data_family
  }
  return(ids_data_family_container)
}

#processing one id only data family 
id_only<-c(3)
ids_container_data_family <- process_data_all_ids_01(data,taxainfo,id_only)
id_only_data_family<-ids_container_data_family[[id_only]]

#process relative abundance 
relative_abundance_function <-function(data_family){
  #computes relative abundance
  numeric_columns <- sapply(data_family, is.numeric)
  #column_sums <- colSums(data_family[-1, numeric_columns])
  column_sums <- colSums(data_family[,numeric_columns])
  relative_abundance <- sweep(data_family[, numeric_columns], 2, column_sums, FUN = "/")
  # Combine with the non-numeric columns (assuming the first column is non-numeric)
  relative_abundance <- cbind(data_family[, !numeric_columns], relative_abundance)
  return(relative_abundance)
}


#rel abun id 1 
data_family_id_one_relative_abundance<-relative_abundance_function(id_only_data_family)

#1% filter of rel abun and also binarized data 
sum_rows_relative_abun_id_one<-rowSums(data_family_id_one_relative_abundance[, -1])
filter<-0.01 
#species rows whose sum > 1%
filtered_data_family_rel_abun_id_one<-data_family_id_one_relative_abundance[sum_rows_relative_abun_id_one>filter,]


#median binarization for each specie 
medians <-apply(id_only_data_family[,-1],1,median)  


data_family_binarized <- id_only_data_family[,-1]
data_family_binarized <- as.data.frame(apply(data_family_binarized, 2, function(x, med) as.numeric(x >= med), med = medians))
data_family_binarized_with_eff_taxo <-cbind(id_only_data_family[,1],data_family_binarized) #86x30


selected_indices<-c(2,12,22,23,35,39,40,41,42,51,58,60,61,62,66,69,74,85)
frequency_data_family_binarized<-rowSums(data_family_binarized_no_eff_taxo)

frequency_data_family_binarized
data_family_binarized_no_eff_taxo<-data_family_binarized #86x29

#transpose
data_family_binarized_with_eff_taxo_t<-t(data_family_binarized_with_eff_taxo) #30x86
data_family_binarized_no_eff_taxo_t <-t(data_family_binarized_no_eff_taxo) #29x86 filtered


transposed_data_05 <- function(transposed_data) {
  #this function cleans na values and converts all nums into numeric
  transposed_data_numeric <- apply(transposed_data, 2, function(x) {
    
    # Check the type of the column and convert appropriately
    if (is.factor(x) || is.character(x)) {
      return(as.numeric(as.character(x)))
    } else if (is.logical(x)) {
      return(as.numeric(x))
    } else if (is.numeric(x)) {
      return(x)  # Leave numeric columns as they are
    } else {
      stop("Unsupported data type found")
    }
  })
  transposed_data_numeric[is.na(transposed_data_numeric)] <- 0
  return(transposed_data_numeric)
}


#30x86
data_family_binarized_with_eff_taxo_t_clean<-transposed_data_05(data_family_binarized_with_eff_taxo_t)
#29x86
data_family_binarized_no_eff_taxo_t_clean<-transposed_data_05(data_family_binarized_no_eff_taxo_t)


pc_algo_discrete <- function(species_indices, data_family_b, data_family_b_t_clean) {
  species_indices<-as.numeric(species_indices)
  
  print("species indices testing")
  print(species_indices)
  selected_data_family_binarized <- data_family_b[species_indices, , drop = FALSE]
  
  print("selected_data_family_binarized testing")
  print(selected_data_family_binarized)
  
  
  
  print("data_family_b_t_clean testing ")
  print(data_family_b_t_clean)
  
  
  selected_data_family_binarized_t_clean <- data_family_b_t_clean[, species_indices, drop = FALSE]
  selected_data_family_binarized_t_clean_temp <- as.matrix(selected_data_family_binarized_t_clean)
  
  
  
  print("selected_data_family_binarized_t_clean_temp testing")
  print(selected_data_family_binarized_t_clean_temp)

  # Run PC algorithm if variability exists
  pc_fit_selected <- pc(suffStat = list(dm = selected_data_family_binarized_t_clean_temp, nlev = rep(2, ncol(selected_data_family_binarized_t_clean_temp)), adaptDF = FALSE),
                        indepTest = binCItest,
                        #indepTest = disCItest, 
                        skel.method = "stable",#additional info
                        labels = rownames(selected_data_family_binarized),
                        alpha = 0.01)
  
  
  
  adj_matrix <- as(pc_fit_selected@graph, "matrix")
  #par(mar = c(1, 1, 1, 1))
  plot(pc_fit_selected@graph)
  
  return(list(pc_fit = pc_fit_selected, adj_matrix = adj_matrix, tranposed_data_family_binarized = selected_data_family_binarized_t_clean_temp, graph = pc_fit_selected@graph))
}



causal_data_all_up_to_pc_10 <- function(data_family_b, data_family_b_t_clean) {
  
  print("inside causal func 10")
  result<-list()
  #selected_indices<-c(2 ,12,22,23,35,39,40,41,42,51,58,60,61,62,66,69,74,85)
  
  #remove 2,35,42,82,69,66, and 22
  
  #useful but modify for each id
  #selected_indices<-c(12,23,39,40,41,51,58,60,61,62,74)
  
  #selected_indices<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
  
  #fix this back to selected_species 
  species_amount<-nrow(data_family_b)
  selected_indices<-sample(c(1:species_amount),17,replace=FALSE)
  
  print("selected_indices")
  print(selected_indices)
  #lowest presence species
  #selected_indices<-c(12,22,23,39,40,41,42,51,58,60,61,62,66,69,74,85)
  while (TRUE) {
    max_indices <- ncol(data_family_b_t_clean) #86
    
    list_size <- sample(8:min(10, max_indices - 1), 1) # Between 3 and 10
    #
    selected_species <- sample(selected_indices, list_size, replace = FALSE)
    
    
    pc_algo_output <- pc_algo_discrete(selected_species, data_family_b, data_family_b_t_clean)
    
    print("pc algo output testing")
    print(pc_algo_output)
    
    if (!is.null(pc_algo_output)) {  # Proceed only if a valid result is obtained
      #hash_map[[id]] <- list(
      result<-list(
        pc_fit = pc_algo_output$pc_fit,
        adj_matrix = pc_algo_output$adj_matrix,
        pc_algo_data_family_binarized_t_clean = pc_algo_output$tranposed_data_family_binarized,
        graph = pc_algo_output$graph
      )
      break
    }
  }
  #return(hash_map)
  return (result)
}

repeated_call_causal_data_all_up_to_pc_101 <- function(iterations, data_family_b, data_family_b_t_clean) {
  pc_container <- list()
  
  for (iter in 1:iterations) { # Computed multiple times (e.g., 1000 times)
    print(paste(iter, "iter"))
    
    causal_output_one_id_iter <- causal_data_all_up_to_pc_10(data_family_b, data_family_b_t_clean)
    
    print("causal_output_one_id_iter testing" )
    print(causal_output_one_id_iter)
    
    pc_container[[iter]] <- list(
      pc_fit = causal_output_one_id_iter$pc_fit,
      adj_matrix = causal_output_one_id_iter$adj_matrix,
      selected_data_numeric = causal_output_one_id_iter$pc_algo_data_family_binarized_t_clean,
      graph = causal_output_one_id_iter$graph
    )
  }
  
  return(pc_container)
}

pc_results <- repeated_call_causal_data_all_up_to_pc_101(10,  # Number of iterations
                                                         #data_family_binarized_with_eff_taxo, data_family_binarized_no_eff_taxo_t_clean, 1)
                                                         data_family_binarized_no_eff_taxo, data_family_binarized_no_eff_taxo_t_clean)




