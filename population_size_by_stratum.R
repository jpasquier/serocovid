setwd("~/Projects/SerocoViD")
pop <- read.table(header = TRUE, text = "
  survey  stratum       N
       1        1   37490
       1        2   42538
       1        3   42700
       1        4   42171
       1        5  212336
       1        6  268898
       1        7   67412
       1        8   63097
       3        4   42140
       3        5  214498
       3        6  270755
       3        7   67462
       3        8   63605
       4        1   38206
       4        2   43179
       4        3   43666
       4        4   42434
       4        5  216839
       4        6  273416
       4        7   67019
       4        8   64770
")
saveRDS(pop, "data/population_size_by_stratum.rds", compress = "xz")
write.table(pop, "data/population_size_by_stratum.csv", sep = ";",
            row.names = FALSE)
# Sources:
# 1 data-fso/POPULATION PAR STRATES.xlsx
# 3 data-fso/POUR_PONDEARTION_JPASQUIER.xlsx
# 4 data-fso/Pour ponderation vague 4.xlsx

