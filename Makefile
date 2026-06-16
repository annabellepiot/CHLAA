.PHONY: document test check build render-vignettes clean-vignettes

R := Rscript --vanilla

document:
	$(R) -e "devtools::document(roclets = c('rd', 'namespace'))"

test:
	$(R) -e "devtools::test()"

check:
	$(R) -e "devtools::check(args = c('--no-manual'), error_on = 'warning')"

build:
	$(R) -e "devtools::build()"

render-vignettes: vignettes/fitting.Rmd

vignettes/fitting.Rmd: vignettes_src/fitting.Rmd
	rm -f vignettes/fitting.md vignettes/fitting.Rmd
	rm -rf vignettes/fitting vignettes_src/fitting
	$(R) -e "if (requireNamespace('pkgload', quietly = TRUE)) pkgload::load_all('.', quiet = TRUE) else devtools::load_all('.', quiet = TRUE); rmarkdown::render('$<', output_format = rmarkdown::md_document(variant = 'markdown_github', preserve_yaml = TRUE), output_file = 'fitting.md', output_dir = 'vignettes', quiet = FALSE, envir = new.env(parent = globalenv()))"
	mv vignettes/fitting.md vignettes/fitting.Rmd
	mv vignettes_src/fitting vignettes/fitting

clean-vignettes:
	rm -f vignettes/fitting.md
	rm -f vignettes/fitting.Rmd
	rm -rf vignettes/fitting
	rm -rf vignettes_src/fitting
