PHONY: clean

clean:
	rm derived_data/*
	rm logs/*
	rm figures/*
	
.created-dirs:
  mkdir -p figures

UMAP_RNA0.8.png: Data/
	Rscript ReadIn.R

report.pdf: /figures/UMAP_RNA0.8.png /figures/UMAP_RNA0.8.rds
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"pdf_document\")"