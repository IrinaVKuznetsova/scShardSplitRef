# vignette that depends on Internet access need to be pre-compiled and takes a
# while to run
library("devtools")
library("knitr")

install() # ensure we're using the latest version of scShardSplitRef

# Precompile the vignette
knit("vignettes/scShardSplitRefQuickStart.Rmd.orig", "vignettes/scShardSplitRefQuickStart.Rmd")

# remove file path such that vignettes will build with figures
replace <- readLines("vignettes/scShardSplitRefQuickStart.Rmd")
replace <- gsub("<img src=\"vignettes/", "<img src=\"", replace)
fileConn <- file("vignettes/scShardSplitRefQuickStart.Rmd")
writeLines(replace, fileConn)
close(fileConn)

# build vignette
build_vignettes()

# move resource files to /docs
resources <-
  list.files("vignettes/", pattern = ".png$", full.names = TRUE)
file.copy(from = resources, to = "docs", overwrite = TRUE)
