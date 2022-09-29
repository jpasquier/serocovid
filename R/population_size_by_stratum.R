setwd("~/Projects/SerocoViD")
pop <- read.table(header = TRUE, text = "
  survey  stratum    age       N
       1        1   6m-4   37490
       1        2    5-9   42538
       1        3  10-14   42700
       1        4  15-19   42171
       1        5  20-39  212336
       1        6  40-64  268898
       1        7  65-74   67412
       1        8   >=75   63097
       3        4  15-19   42140
       3        5  20-39  214498
       3        6  40-64  270755
       3        7  65-74   67462
       3        8   >=75   63605
       4        1   6m-4   38206
       4        2    5-9   43179
       4        3  10-14   43666
       4        4  15-19   42434
       4        5  20-39  216839
       4        6  40-64  273416
       4        7  65-74   67019
       4        8   >=75   64770
       6        1  15-29  142257
       6        2  30-44  175251
       6        3  45-64  217352
       6        4   >=65  134481
")
saveRDS(pop, "data/population_size_by_stratum.rds", compress = "xz")
write.table(pop, "data/population_size_by_stratum.csv", sep = ";",
            row.names = FALSE)
# Sources:
# 1 data-fso/POPULATION PAR STRATES.xlsx
# 3 data-fso/POUR_PONDEARTION_JPASQUIER.xlsx
# 4 data-fso/Pour ponderation vague 4.xlsx
# 6 data-fso/MARGES_VD-VAGUE5.xlsx

