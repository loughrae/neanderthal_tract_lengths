library(tidyverse)
gen <- round(c(56000/29, 47000/29, 65000/29), digits = 0)
stime <- seq(from = 344, to = max(gen), by = 10)
s <- seq(from = 0.001, to = 0.02, by = 0.001)
combos <- expand.grid(s = s, t1 = stime, t2 = gen)
head(combos)
str(combos)
combos <- combos[which(combos$t1 <= combos$t2),]
combos$t1 <- ifelse(combos$t1 == combos$t2, combos$t1 - 1, combos$t1)
combos$row <- 1:nrow(combos)
cores <- sample.int(208, nrow(combos), replace = TRUE)
table(cores)
combos <- cbind(combos, cores)
lis <- split(combos, combos$cores)
for (i in names(lis)) {
  lis[[i]] %>% as.data.frame() %>% select(s, t1, t2) %>% write.table(file = paste0('allcombos_nov2021_', i, '.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}
