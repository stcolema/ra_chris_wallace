

library(data.table)
library(magrittr)

timecourse_norm <- fread("../Yeast/Input_data/Reduced/Granovskaia_timecourse_normalised_reduced.csv", header = T)
harbison <- fread("../Yeast/Input_data/Reduced/harbison_marina.csv", header = T)
ppi <- fread("../Yeast/Input_data/Reduced/ppi.csv", header = T)

nrow(timecourse_norm)

all(timecourse_norm$V1 == harbison$V1)
all(timecourse_norm$V1 == ppi$V1)
all(ppi$V1 == harbison$V1)

timecourse_full <- fread("../Yeast/Input_data/marina_alpha_orfs.tsv", header = T)

timecourse_reduced <- timecourse_full[timecourse_full$`# Feature` %in% timecourse_norm$V1,]

timecourse_order <- match(timecourse_norm$V1, timecourse_reduced$`# Feature`)

timecourse_ordered <- timecourse_reduced[timecourse_order, ]
colnames(timecourse_ordered) <- c("V1", colnames(timecourse_ordered)[-1])

all(timecourse_norm$V1 == harbison$V1)
all(timecourse_norm$V1 == ppi$V1)
all(timecourse_norm$V1 == timecourse_ordered$V1)
all(ppi$V1 == harbison$V1)

cols_to_keep <- which(c(T, mod(2:42, 2) == 0))

fwrite(timecourse_ordered[, ..cols_to_keep], file = "../Yeast/Input_data/Granovskaia_time_course_551.csv")
