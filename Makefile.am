# Makefile.am

document:
	Rscript -e "devtools::document(roclets = c('rd', 'collate', 'namespace'))"

buildbin:
	Rscript -e "devtools::build(binary = TRUE, args = c('--preclean'))"

buildsrc:
	Rscript -e "devtools::build(args=c('--no-build-vignettes'))"

install:
	Rscript -e "devtools::install_local(force=TRUE, upgrade=FALSE)"

test:
	Rscript -e "devtools::test()"

site:
	Rscript -e "pkgdown::build_site(preview=FALSE)"


all-local: document buildsrc buildbin install test site
