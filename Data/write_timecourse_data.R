

library(data.table)
library(magrittr)

curr_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/Yeast/Input_data/"

timecourse_norm <- fread(paste0(curr_dir, "Reduced/Granovskaia_timecourse_normalised_reduced.csv"), header = T)
harbison <- fread(paste0(curr_dir, "Reduced/harbison_marina.csv"), header = T)
ppi <- fread(paste0(curr_dir, "Reduced/ppi.csv"), header = T)

nrow(timecourse_norm)

all(timecourse_norm$V1 == harbison$V1)
all(timecourse_norm$V1 == ppi$V1)
all(ppi$V1 == harbison$V1)

timecourse_full <- fread(paste0(curr_dir, "marina_alpha_orfs.tsv"), header = T)

timecourse_reduced <- timecourse_full[timecourse_full$`# Feature` %in% timecourse_norm$V1,]

timecourse_order <- match(timecourse_norm$V1, timecourse_reduced$`# Feature`)

timecourse_ordered <- timecourse_reduced[timecourse_order, ]
colnames(timecourse_ordered) <- c("V1", colnames(timecourse_ordered)[-1])

all(timecourse_norm$V1 == harbison$V1)
all(timecourse_norm$V1 == ppi$V1)
all(timecourse_norm$V1 == timecourse_ordered$V1)
all(ppi$V1 == harbison$V1)

cols_to_keep <- which(c(T, mod(2:42, 2) == 0))

fwrite(timecourse_ordered, file = paste0(curr_dir, "Granovskaia_time_course_551.csv"))
fwrite(timecourse_ordered[, ..cols_to_keep], file = paste0(curr_dir, "Granovskaia_time_course_551_less_time_points.csv"))
