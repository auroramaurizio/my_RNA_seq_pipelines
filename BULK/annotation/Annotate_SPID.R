library("hypeR")
library("assertr")
library("tidyr")

options(warn=-1)

Annotate_SPID <- function(DGE, enrich.database  = "WikiPathway_2021_Human", output_tsv = F) {
  
  annoname = enrich.database
  annotation_table = enrichr_download(enrich.database)
  annotation_table = as.data.frame(do.call(rbind, annotation_table)) #here gives a worning ignore it
  colNames = colnames(annotation_table) # could be any number of column names here
  annotation_table['test']=col_concat(annotation_table, sep = " ")
  annotation_table['GeneID'] <- trimws(annotation_table$test, which = c("both")) #remove whitespaces
  annotation_table$term = row.names(annotation_table)
  annotation_table_sub=annotation_table[,c("term","GeneID")]
  colnames(annotation_table_sub)  
  exploded = unique(separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  #group by gene.. that is: one gene x row with associated all the descriptions (space separated)
  grouped  <- exploded %>%
    group_by(GeneID) %>%
    summarise(temp = toString(term)) %>%
    ungroup()
  
  nrow(grouped)
  colnames(grouped)=c("GeneID",annoname)
  
  merged <- merge(DGE, grouped, by="GeneID", all.x = TRUE)
  
  if (output_tsv) {
    write.table(merged, paste(annoname, "_SP_annotation.tsv", sep ="_"), quote = F, sep = "\t")
  }
  return(merged)
}

merged = Annotate_SPID(DGE, "OMIM_Disease")
head(merged)
