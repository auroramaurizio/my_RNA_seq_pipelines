# Cellranger guidelines to insert a transgene in your reference genome 

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr#marker
From this link go to "Add a Marker Gene to the FASTA and GTF" and follow the instructions:

1) get the  FASTA sequence of the transgene ypu want to insert (e.g. TdTomato, EYFP, TCL1A etc...)

2) append it to the FASTA of your reference eg. mm10 

3) append to the GTF reference, the GTF of your transgene 

4) activate a conda env with cellranger: e.g. conda activate demux

4.1) maybe you can skip this point but I'll put it here just in case just in case https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
cellranger mkgtf input.gtf output.gtf --attribute=key:allowable_value
e.g. cellranger mkgtf  TdTomato.basic.annotation.gtf TdTomato.basic.annotation.gtf.cellranger3 --attribute=key:allowable_value

5) run the command: cellranger mkref --genome=output_genome --fasta=input.fa --genes=input.gtf
input.fa è il tuo nuovo fasta con il transgene e input.gtf è il tuo nuovo gtf con il transgene
e.g.: cellranger mkref --genome=TdTomato.GRCm38.primary_assembly.cellranger3_index --fasta=TdTomato.GRCm38.primary_assembly.genome.fa --genes=TdTomato.basic.annotation.gtf.cellranger3 --nthreads 12

Done!


