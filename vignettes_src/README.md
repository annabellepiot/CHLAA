# Heavy Source Vignettes

This folder contains source documents for heavyweight vignettes that should not
be rebuilt during routine `R CMD check`.

Render the precomputed Markdown artefacts with:

```bash
make render-vignettes
```

The Makefile renders `vignettes_src/fitting.Rmd` to `vignettes/fitting.md`,
then renames that rendered Markdown artefact to `vignettes/fitting.Rmd`.
Figures are written next to the rendered output. The `vignettes_src/` directory
itself is excluded from package builds via `.Rbuildignore`.
