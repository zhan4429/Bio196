# Bio196-22 Computing for Bioinformatics

*Linux, HPC, environments, workflows, and reproducible analysis.*

This is the course reader for Bio196-22, a six-week computing-for-bioinformatics course.
It covers Linux basics and text editing, SSH and HPC, Conda and containers, NGS and
RNA-seq, Nextflow and nf-core, version control, and reproducible reports.

## Repository Layout

```
.
├── _quarto.yml         # Quarto book configuration
├── index.qmd           # Preface
├── references.qmd      # Bibliography page
├── references.bib      # BibTeX bibliography
├── chapters/           # Chapter sources
└── images/             # Figures referenced from chapters
```

## Building

You need [Quarto](https://quarto.org/) >= 1.4 installed.

```bash
# Render the full book (HTML and PDF if a LaTeX engine is present)
quarto render

# Live-reloading preview while writing
quarto preview
```

The HTML output is written to `_book/`.

## Contributing a Chapter

1. Add the chapter `.qmd` file to `chapters/`.
2. Add or update any figures in `images/`.
3. Add bibliography entries to `references.bib` (BibTeX format).
4. Add the chapter to `_quarto.yml` under the appropriate Part.
5. Build locally and check the rendered output before committing.

## License

See `LICENSE` (or contact the author).
