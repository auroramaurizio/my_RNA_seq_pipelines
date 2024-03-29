#head file.bed 
#chr1	11868	14409	ENSG00000223972.5	0	+	11868	11868	255,0,0	1	2540	11868
#chr1	11868	12227	ENSG00000223972.5	0	+	11868	11868	255,0,0	1	358	11868
#chr1	12612	12721	ENSG00000223972.5	0	+	12612	12612	255,0,0	1	108	12612


#tree TIN
#├── PTCL_0004_RNA.bam
#├── PTCL_0004_RNA.bam.bai
#├── PTCL_0010_RNA.bam
#├── PTCL_0010_RNA.bam.bai

tin.py -i TIN -r file.bed -c 5 -n 50 -s

# calculate-tin.py with respect to tin.py takes advantage of parallelization

calculate-tin.py -i file.bam -r file.bed -c 5 -n 50 -s -p 10
