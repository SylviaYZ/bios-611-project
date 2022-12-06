.PHONY: clean

clean:
	rm -rf figures
	rm -rf .created-dirs
	rm -f report.pdf
	rm -f report.html
	rm -f report.tex
	rm -f report.log
	rm -f Rplots.pdf
	
.created-dirs:
	mkdir -p figures
	touch .created-dirs

# RNA QC + RNA processing
# Build UMAP on RNA with resolution = 0.5 for clustering
~/work/figures/UMAP_RNA0.5.png\
~/work/figures/UMAP_RNA0.5.rds\
~/work/Data/Cov2RNA_clustered.rds\
~/work/figures/violin_RNA_QC.rds: .created-dirs ~/work/Data/cov02
	Rscript RNA_processing.R
	
# RNA violin plots for each cluster on given DE genes
~/work/figures/RNA_violin_T.rds\
~/work/figures/RNA_violin_T4_T8.rds\
~/work/figures/RNA_violin_B.rds\
~/work/figures/RNA_violin_NK.rds\
~/work/figures/RNA_violin_CD14mono.rds\
~/work/figures/RNA_violin_DC.rds\
~/work/figures/RNA_violin_PL.rds: .created-dirs ~/work/Data/Cov2RNA_clustered.rds
	Rscript RNA_cluster_violin_plots.R
	
# RNA DE analysis for each cluster using Wilcoxon rank sum test
~/work/figures/RNA_DEcluster_heatmap.rds: .created-dirs ~/work/Data/Cov2RNA_clustered.rds
	Rscript RNA_DE_Wilcoxon.R
	
# Explore ADT
~/work/figures/ADT_hist_CD3.rds\
~/work/figures/RNA_hist_CD3.rds\
~/work/Data/ADT_10sample.rds: .created-dirs ~/work/Data/Cov2RNA_clustered.rds
	Rscript ADT_explore.R
	
# Processing ADT
~/work/Data/ADT_10sample_normalzied.rds\
~/work/figures/ADT_violin_T.rds\
~/work/figures/ADT_violin_T4_T8.rds\
~/work/figures/ADT_violin_T4naive.rds\
~/work/figures/ADT_violin_T4memeff.rds\
~/work/figures/ADT_violin_B.rds\
~/work/figures/ADT_violin_NK.rds\
~/work/figures/ADT_violin_CD14mono.rds\
~/work/figures/ADT_violin_CD14mono.rds: .created-dirs ~/work/Data/Cov2RNA_clustered.rds
	Rscript ADT_processing.R

# Build report
report.html: .created-dirs ~/work/figures/UMAP_RNA0.5.rds\
	~/work/figures/violin_RNA_QC.rds\
	~/work/figures/RNA_violin_T.rds\
	~/work/figures/RNA_violin_T4_T8.rds\
	~/work/figures/RNA_violin_B.rds\
  ~/work/figures/RNA_violin_NK.rds\
  ~/work/figures/RNA_violin_CD14mono.rds\
  ~/work/figures/RNA_violin_DC.rds\
  ~/work/figures/RNA_DEcluster_heatmap.rds\
  ~/work/figures/ADT_hist_CD3.rds\
  ~/work/figures/RNA_hist_CD3.rds\
  ~/work/Data/ADT_10sample.rds\
  ~/work/Data/ADT_10sample_normalzied.rds\
	~/work/figures/ADT_violin_T.rds\
	~/work/figures/ADT_violin_T4_T8.rds\
	~/work/figures/ADT_violin_T4naive.rds\
	~/work/figures/ADT_violin_T4memeff.rds\
	~/work/figures/ADT_violin_B.rds\
	~/work/figures/ADT_violin_NK.rds\
	~/work/figures/ADT_violin_CD14mono.rds\
	~/work/figures/ADT_violin_DC.rds\
	~/work/UMAP_celltypeRNA.png\
	~/work/UMAP_celltypeADT.png
		R -e "rmarkdown::render(\"report.Rmd\", output_format=\"html_document\")"

	
	
	
	