source("src/functions_smake_test.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_wildcard_n <- args[1]
args_input <- args[2]
args_output <- args[3]
#### test args####
# args_wildcard_n <- c("3")
# args_input <- c("data/smake_test/input_number.csv")
# args_output <- c("data/smake_test/3_df.csv")

input_number <- read_csv(args_input, col_names = FALSE)

sum_2000 <- as.numeric(input_number) + as.numeric(args_wildcard_n)
sum_2000 <- as.data.frame(sum_2000)
write_csv(sum_2000, args_output)