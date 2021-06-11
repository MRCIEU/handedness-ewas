

unrelated.children <- read.table(file.path(data.dir, "alspac/related/children/inclusion_list.txt"), stringsAsFactors=F)
unrelated.children <- unrelated.children[,1]

unrelated.mothers <- read.table(file.path(data.dir, "alspac/related/mothers/inclusion_list.txt"), stringsAsFactors=F)
unrelated.mothers <- unrelated.mothers[,1]
