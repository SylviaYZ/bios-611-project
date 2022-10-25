.PHONY: clean

clean:
	rm -rf figures
	rm -rf .created-dirs
	rm -f report.pdf
	rm -f report.html
	rm -f report.tex
	rm -f report.log
	
.created-dirs:
	mkdir -p figures
	touch .created-dirs

# Build UMAP on RNA with resolution = 0.8 for clustering
~/figures/UMAP_RNA0.8.png ~/figures/UMAP_RNA0.8.rds: .created-dirs ~/Data/cov01
	Rscript ReadIn.R

# Build report
report.html: .created-dirs ~/figures/UMAP_RNA0.8.rds
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"html_document\")"

report.pdf: .created-dirs ~/figures/UMAP_RNA0.8.rds
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"pdf_document\")"