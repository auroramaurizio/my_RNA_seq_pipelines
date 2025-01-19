# https://cnvkit.readthedocs.io/en/stable/rna.html

cnvkit.py import-rna *_Counts.txt \
--normal TCL_Bcell3669_9_Counts.txt TCL_Bcell3670_10_Counts.txt TCL_Bcell3671_11_Counts.txt TCL_Bcell3672_12_Counts.txt \
--gene-resource ensembl_gene_info.tsv \
--output runGSE228976norm4summary.tsv --output-dir runGSE228976norm4noGC
