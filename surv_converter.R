## script of surv_converter() function
## Converts HMI survival data to long format. Output file has information on flies alive per day
# and number of flies with different status for survival analysis.

## Author: Luis Teixeira, Host-Microorganism Interactions lab
##version 1.0

# requires tidyverse library

surv_converter <- function(filename, save_table = TRUE) {
  
  #load file
  surv <- read.csv(filename, stringsAsFactors = FALSE, header = FALSE)
  
  
  # change entries using higher case in old format to lower case. These lines are only required to make this script compatible with older data files.
  #It can be eliminated in a new script and an error will be shown if old format is being used. Directly editing data file will also solve the problem.
  surv[surv == "Factors"] <- c("factors")
  surv[surv == "Days"] <- c("day")
  surv[surv == "Info"] <- c("info")
  surv[surv == "Date"] <- c("date")
  
  #find rows which separate info, factors, and survival data blocks.
  factors_row <- match("factors",surv[,1])
  day_row <- match("day",surv[,1])
  
  # Get date, last day of survival, total number of flies, names of factors
  date <- surv[match("date",surv[,1]), 2]
  lastday <- as.numeric(surv[nrow(surv), 1])
  sum_control <- sum(as.numeric(surv[day_row+1, 2:ncol(surv)]))
  factors_names <- surv[(factors_row+1):(day_row-1), 1]
  
  #print total number of flies
  print(paste(("total number of flies is "),sum_control,sep=""))
  
  # quality check #####
  
  #check some file formatting
  if (surv[1,1] != "info") {stop("There is an error in  file format")}
  if (is.na(factors_row)) {stop("There is an error in  file format, no cell with *factors*")}
  if (is.na(day_row)) {stop("There is an error in  file format, no cell with *days*")}
  
  if (day_row < factors_row) {stop("There is an error in  file format")}
  if (day_row == (factors_row+1)) {stop("There is an error in  file format - factors are not defined")}
  
  
  #print row and column number to check in case of troubleshooting
  print(paste("input file has ", nrow(surv), " rows and ", ncol(surv), " columns", sep=""))
  
  #check if number of flies increases from one day to the other
  for (j in 2:ncol(surv)) {
    for (i in (day_row+1):(nrow(surv)-1)) {
      if (as.numeric(surv[i,j]) < as.numeric(surv[i+1,j])) {
        stop(paste("Spontaneous generation of flies at day ", surv[i+1,1],", column ", j, sep=""))
      }
    }
  }
  
  #check if days start at 0 and are consecutive numbers till last day
  if (surv[day_row+1, 1] != 0) {stop("There is an error in  file format - day 0 data is not shown")}
  
  for (i in (day_row+1):(nrow(surv)-1)) {
    if ((as.numeric(surv[i+1,1]) - as.numeric(surv[i,1])) != 1 ) {
      stop(paste("There is an error in days of survival, they should be consecutive days. Check rows ", i," and ", i+1, sep=""))
    }
  }
  
  
  # transform data #####
  
  #transpose table, get column names from first row, delete first row, delete columns of "info" block and column "day"
  trans_surv <- as.data.frame(t(surv))
  colnames(trans_surv)=trans_surv[1,]
  trans_surv <- trans_surv[-1, ]
  trans_surv <- trans_surv[, -c(1:factors_row,day_row)]
  
  #convert survival data to numeric
  data_col <- ((day_row - factors_row):ncol(trans_surv))
  trans_surv[, data_col] <- sapply(trans_surv[, data_col], as.numeric)
  
  #check that each row is a unique combination of factors
  if (nrow(trans_surv) != nrow(unique(trans_surv[, factors_names]))) {
    stop("There is an error in factors specification - at least two column have the same combination of factors")
  }
  
  #transform wide to long
  trans_surv <- trans_surv %>% 
    pivot_longer(
      cols = all_of(data_col),
      names_to = "day_alive",
      values_to = "number_alive",
      values_drop_na = FALSE
    )
  
  #add date if exists
  if (!is.na(date)) {trans_surv <- add_column(trans_surv, date = date, .before = 1)}
  
  
  #make "day" and "survival" numeric
  trans_surv$day_alive <- as.numeric(trans_surv$day_alive)
  trans_surv$number_alive <- as.numeric(trans_surv$number_alive)
  
  
  #make day_status, number_status and status columns
  trans_surv <- trans_surv %>% add_column(day_status = NA, number_status = NA, status = NA)
  
  
  #get dead per day, add 1 to day with dead flies, set status for statistical analysis (1  = dead, 0 = alive (truncated at the last day))
  for (i in 1:(nrow(trans_surv))){
    if (trans_surv$day_alive[i] != lastday) {
      trans_surv$number_status[i] = trans_surv$number_alive[i] - trans_surv$number_alive[i+1]
      trans_surv$day_status[i] <- trans_surv$day_alive[i] + 1
      if (trans_surv$number_status[i] != 0) {
        trans_surv$status[i] <- 1
      }
    } else {
      trans_surv$number_status[i] = trans_surv$number_alive[i]
      trans_surv$day_status[i] = trans_surv$day_alive[i]
      if (trans_surv$number_status[i] != 0) {
        trans_surv$status[i] <- 0
      }
    }
  }
  
  
  # check if final number of flies is equal to inicial number of flies 
  if (sum(trans_surv$number_status) != sum_control) {  stop("There is an error number of flies in final table") }
  
  
  #write file
  if(save_table) {
    filename <- str_remove(filename, ".csv")
    write.table(trans_surv, file = paste(filename, "_converted.csv", sep=""), sep = ",", row.names = FALSE, col.names = TRUE)
  }
  
  #end
  return(trans_surv)
  
}