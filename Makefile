### Makefile template from https://github.com/lgaborini/R-package-makefile/blob/main/Makefile
SCRIPTS := $(wildcard R/*.R)
MAN := $(wildcard man/*.Rd)

PKG := $(shell head -1 DESCRIPTION | sed 's/Package: //')
VERSION := $(shell sed -n 3p DESCRIPTION | sed 's/Version: //')
PKG_SOURCE := ../$(PKG)_$(VERSION).tar.gz
PKG_BINARY := ../$(PKG)_$(VERSION).zip


all: document data source_package binary_package packages test install pkgdown #drat

.PHONY: document data source_package binary_package packages targets all install test pkgdown #drat

inst/extdata/10x/%.csv: inst/extdata/10x
	$(shell wget --no-check-certificate "https://cdn.10xgenomics.com/raw/upload/v1655151897/support/in-line documents/$(@F)" -O $@)

inst/extdata/10x: inst/extdata
	$(shell	mkdir inst/extdata/10x)

inst/extdata: inst
	$(shell	mkdir inst/extdata)
inst:
	$(shell mkdir inst)

data: inst/extdata/10x/Dual_Index_Kit_TT_Set_A.csv inst/extdata/10x/Dual_Index_Kit_NN_Set_A.csv inst/extdata/10x/Dual_Index_Kit_NT_Set_A.csv inst/extdata/10x/Single_Index_Kit_T_Set_A.csv inst/extdata/10x/Single_Index_Kit_N_Set_A.csv

targets:
	$(info Package : $(PKG))
	$(info Version : $(VERSION))
	$(info Source  : $(PKG_SOURCE))
	$(info Binary  : $(PKG_BINARY))
	$(info Rds     : $(MAN))
	$(info Scripts : $(SCRIPTS))

# Requires GNU Make 4.3
#
# https://www.gnu.org/software/make/manual/html_node/Multiple-Targets.html
# &: groups all targets as if they were one
#
# It builds the help files at most once (instead of once per .Rd)
$(MAN) &: $(SCRIPTS)
	$(info ** Running devtools::document **)
	@R --slave -e "devtools::document(roclets = c('rd', 'collate', 'namespace'))"

document: $(MAN)

$(PKG_SOURCE): $(MAN)
	$(info ** Making source package $(PKG_SOURCE) **)
	R --slave -e "devtools::build(args = c('--no-build-vignettes'))"

$(PKG_BINARY): $(MAN)
	$(info ** Making binary package $(PKG_BINARY) **)
	R --slave -e "devtools::build(binary = TRUE, args = c('--preclean'))"

install: $(PKG_SOURCE) $(PKG_BINARY)
	$(info ** Installing package **)
	R --slave -e "devtools::install_local(force=TRUE)"

test: $(PKG_SOURCE) $(PKG_BINARY)
	$(info ** Testing package **)
	R --slave -e "devtools::test()"

source_package: $(PKG_SOURCE)

binary_package: $(PKG_BINARY)

packages: document source_package binary_package

pkgdown: document
	$(info ** Making pkgdown site **)
	R --slave -e "pkgdown::build_site(preview = FALSE)"

# drat: packages
# 	$(info ** Pushing packages to drat repo **)
# 	R --slave -e "drat::insertPackage('$(PKG_SOURCE)')"
# 	R --slave -e "drat::insertPackage('$(PKG_BINARY)')"
