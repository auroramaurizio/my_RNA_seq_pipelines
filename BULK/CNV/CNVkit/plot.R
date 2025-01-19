setwd("/Users/maurizio.aurora/brendy/CNVKIT/runGSE228976norm4/ok")
getwd()

CNVKITBcell3658 <-read.table("TCL_Bcell3658_1_Counts_call.cns", header=TRUE);
CNVKITBcell3658 <- CNVKITBcell3658[CNVKITBcell3658$probes > 100,]
CNVKITBcell3658 <- CNVKITBcell3658[,c("chromosome","start","end","cn")]
CNVKITBcell3658$strand <- "*"
CNVKITBcell3658$width  <- CNVKITBcell3658$end - CNVKITBcell3658$start
CNVKITBcell3658 <- CNVKITBcell3658[CNVKITBcell3658$width >= 1000000,]
colnames(CNVKITBcell3658) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3658$seqnames <- paste0("chr", CNVKITBcell3658$seqnames)
CNVKITBcell3658 <- CNVKITBcell3658[!grepl("chrY", CNVKITBcell3658$seqnames), ]

CNVKITBcell3662 <-read.table("TCL_Bcell3662_4_Counts_call.cns", header=TRUE);
CNVKITBcell3662 <- CNVKITBcell3662[CNVKITBcell3662$probes > 100,]
CNVKITBcell3662 <- CNVKITBcell3662[,c("chromosome","start","end","cn")]
CNVKITBcell3662$strand <- "*"
CNVKITBcell3662$width  <- CNVKITBcell3662$end - CNVKITBcell3662$start
CNVKITBcell3662 <- CNVKITBcell3662[CNVKITBcell3662$width >= 1000000,]
colnames(CNVKITBcell3662) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3662$seqnames <- paste0("chr", CNVKITBcell3662$seqnames)
CNVKITBcell3662 <- CNVKITBcell3662[!grepl("chrY", CNVKITBcell3662$seqnames), ]


CNVKITBcell3668 <-read.table("TCL_Bcell3668_8_Counts_call.cns", header=TRUE);
CNVKITBcell3668 <- CNVKITBcell3668[CNVKITBcell3668$probes > 100,]
CNVKITBcell3668 <- CNVKITBcell3668[,c("chromosome","start","end","cn")]
CNVKITBcell3668$strand <- "*"
CNVKITBcell3668$width  <- CNVKITBcell3668$end - CNVKITBcell3668$start
CNVKITBcell3668 <- CNVKITBcell3668[CNVKITBcell3668$width >= 1000000,]
colnames(CNVKITBcell3668) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3668$seqnames <- paste0("chr", CNVKITBcell3668$seqnames)
CNVKITBcell3668 <- CNVKITBcell3668[!grepl("chrY", CNVKITBcell3668$seqnames), ]

CNVKITBcell3667 <-read.table("TCL_Bcell3667_7_Counts_call.cns", header=TRUE);
CNVKITBcell3667 <- CNVKITBcell3667[CNVKITBcell3667$probes > 100,]
CNVKITBcell3667 <- CNVKITBcell3667[,c("chromosome","start","end","cn")]
CNVKITBcell3667$strand <- "*"
CNVKITBcell3667$width  <- CNVKITBcell3667$end - CNVKITBcell3667$start
CNVKITBcell3667 <- CNVKITBcell3667[CNVKITBcell3667$width >= 1000000,]
colnames(CNVKITBcell3667) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3667$seqnames <- paste0("chr", CNVKITBcell3667$seqnames)
CNVKITBcell3667 <- CNVKITBcell3667[!grepl("chrY", CNVKITBcell3667$seqnames), ]


CNVKITBcell3670 <-read.table("TCL_Bcell3670_10_Counts_call.cns", header=TRUE);
CNVKITBcell3670 <- CNVKITBcell3670[CNVKITBcell3670$probes > 100,]
CNVKITBcell3670 <- CNVKITBcell3670[,c("chromosome","start","end","cn")]
CNVKITBcell3670$strand <- "*"
CNVKITBcell3670$width  <- CNVKITBcell3670$end - CNVKITBcell3670$start
CNVKITBcell3670 <- CNVKITBcell3670[CNVKITBcell3670$width >= 1000000,]
colnames(CNVKITBcell3670) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3670$seqnames <- paste0("chr", CNVKITBcell3670$seqnames)
CNVKITBcell3670 <- CNVKITBcell3670[!grepl("chrY", CNVKITBcell3670$seqnames), ]


CNVKITBcell3660 <-read.table("TCL_Bcell3660_2_Counts_call.cns", header=TRUE);
CNVKITBcell3660 <- CNVKITBcell3660[CNVKITBcell3660$probes > 100,]
CNVKITBcell3660 <- CNVKITBcell3660[,c("chromosome","start","end","cn")]
CNVKITBcell3660$strand <- "*"
CNVKITBcell3660$width  <- CNVKITBcell3660$end - CNVKITBcell3660$start
CNVKITBcell3660 <- CNVKITBcell3660[CNVKITBcell3660$width >= 1000000,]
colnames(CNVKITBcell3660) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3660$seqnames <- paste0("chr", CNVKITBcell3660$seqnames)
CNVKITBcell3660 <- CNVKITBcell3660[!grepl("chrY", CNVKITBcell3660$seqnames), ]

CNVKITBcell3665 <-read.table("TCL_Bcell3665_5_Counts_call.cns", header=TRUE);
CNVKITBcell3665 <- CNVKITBcell3665[CNVKITBcell3665$probes > 100,]
CNVKITBcell3665 <- CNVKITBcell3665[,c("chromosome","start","end","cn")]
CNVKITBcell3665$strand <- "*"
CNVKITBcell3665$width  <- CNVKITBcell3665$end - CNVKITBcell3665$start
CNVKITBcell3665 <- CNVKITBcell3665[CNVKITBcell3665$width >= 1000000,]
colnames(CNVKITBcell3665) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3665$seqnames <- paste0("chr", CNVKITBcell3665$seqnames)
CNVKITBcell3665 <- CNVKITBcell3665[!grepl("chrY", CNVKITBcell3665$seqnames), ]


CNVKITBcell3671 <-read.table("TCL_Bcell3671_11_Counts_call.cns", header=TRUE);
CNVKITBcell3671 <- CNVKITBcell3671[CNVKITBcell3671$probes > 100,]
CNVKITBcell3671 <- CNVKITBcell3671[,c("chromosome","start","end","cn")]
CNVKITBcell3671$strand <- "*"
CNVKITBcell3671$width  <- CNVKITBcell3671$end - CNVKITBcell3671$start
CNVKITBcell3671 <- CNVKITBcell3671[CNVKITBcell3671$width >= 1000000,]
colnames(CNVKITBcell3671) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3671$seqnames <- paste0("chr", CNVKITBcell3671$seqnames)
CNVKITBcell3671 <- CNVKITBcell3671[!grepl("chrY", CNVKITBcell3671$seqnames), ]

CNVKITBcell3661 <-read.table("TCL_Bcell3661_3_Counts_call.cns", header=TRUE);
CNVKITBcell3661 <- CNVKITBcell3661[CNVKITBcell3661$probes > 100,]
CNVKITBcell3661 <- CNVKITBcell3661[,c("chromosome","start","end","cn")]
CNVKITBcell3661$strand <- "*"
CNVKITBcell3661$width  <- CNVKITBcell3661$end - CNVKITBcell3661$start
CNVKITBcell3661 <- CNVKITBcell3661[CNVKITBcell3661$width >= 1000000,]
colnames(CNVKITBcell3661) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3661$seqnames <- paste0("chr", CNVKITBcell3661$seqnames)
CNVKITBcell3661 <- CNVKITBcell3661[!grepl("chrY", CNVKITBcell3661$seqnames), ]

CNVKITBcell3666 <-read.table("TCL_Bcell3666_6_Counts_call.cns", header=TRUE);
CNVKITBcell3666 <- CNVKITBcell3666[CNVKITBcell3666$probes > 100,]
CNVKITBcell3666 <- CNVKITBcell3666[,c("chromosome","start","end","cn")]
CNVKITBcell3666$strand <- "*"
CNVKITBcell3666$width  <- CNVKITBcell3666$end - CNVKITBcell3666$start
CNVKITBcell3666 <- CNVKITBcell3666[CNVKITBcell3666$width >= 1000000,]
colnames(CNVKITBcell3666) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3666$seqnames <- paste0("chr", CNVKITBcell3666$seqnames)
CNVKITBcell3666 <- CNVKITBcell3666[!grepl("chrY", CNVKITBcell3666$seqnames), ]


CNVKITBcell3669 <-read.table("TCL_Bcell3669_9_Counts_call.cns", header=TRUE);
CNVKITBcell3669 <- CNVKITBcell3669[CNVKITBcell3669$probes > 100,]
CNVKITBcell3669 <- CNVKITBcell3669[,c("chromosome","start","end","cn")]
CNVKITBcell3669$strand <- "*"
CNVKITBcell3669$width  <- CNVKITBcell3669$end - CNVKITBcell3669$start
CNVKITBcell3669 <- CNVKITBcell3669[CNVKITBcell3669$width >= 1000000,]
colnames(CNVKITBcell3669) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3669$seqnames <- paste0("chr", CNVKITBcell3669$seqnames)
CNVKITBcell3669 <- CNVKITBcell3669[!grepl("chrY", CNVKITBcell3669$seqnames), ]


CNVKITBcell3672 <-read.table("TCL_Bcell3672_12_Counts_call.cns", header=TRUE);
CNVKITBcell3672 <- CNVKITBcell3672[CNVKITBcell3672$probes > 100,]
CNVKITBcell3672 <- CNVKITBcell3672[,c("chromosome","start","end","cn")]
CNVKITBcell3672$strand <- "*"
CNVKITBcell3672$width  <- CNVKITBcell3672$end - CNVKITBcell3672$start
CNVKITBcell3672 <- CNVKITBcell3672[CNVKITBcell3672$width >= 1000000,]
colnames(CNVKITBcell3672) <- c("seqnames","start","end","cn","strand","length")
CNVKITBcell3672$seqnames <- paste0("chr", CNVKITBcell3672$seqnames)
CNVKITBcell3672 <- CNVKITBcell3672[!grepl("chrY", CNVKITBcell3672$seqnames), ]
########


cn.calls <- list(
  "Bcell3658" = CNVKITBcell3658,
  "Bcell3660" = CNVKITBcell3660,
  "Bcell3661" = CNVKITBcell3661,
  "Bcell3662" = CNVKITBcell3662,
  "Bcell3665" = CNVKITBcell3665,
  "Bcell3666" = CNVKITBcell3666,
  "Bcell3667" = CNVKITBcell3667,
  "Bcell3668" = CNVKITBcell3668,
  "Bcell3669" = CNVKITBcell3669,
  "Bcell3670" = CNVKITBcell3670,
  "Bcell3671" = CNVKITBcell3671,
  "Bcell3672" = CNVKITBcell3672
)



pdf("plot_CNVKIT_Bcell3662lrr.pdf", 15, 6)
kp <- plotKaryotype(genome="mm10", chromosomes=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"), plot.type = 4, main = "CELL3662")
#plotCopyNumberCalls(kp, CNVKITBcell3662, r0=0, r1=0.20,cn.colors = c("#0000FF", "lightblue", "white", "pink", "#FFAAAA", "red", "#bf0000","#9F1111","#911919","#751515","#431919"))
plotLRR(kp, LRR.calls.CNVKIT3662, points.cex = 2, axis.cex = 2)
#cn.cols <- getCopyNumberColors(colors = c("#0000FF", "lightblue", "white", "pink", "#FFAAAA", "red", "#bf0000","#9F1111","#911919","#751515","#431919"))
#legend("center", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols))
dev.off()


pdf("CN_summary_test_1M_cnvkit_GSE228976_normok.pdf", 20, 8)
kp <- plotKaryotype(genome="mm10", chromosomes=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19", "chrX"), plot.type = 4)
plotCopyNumberCalls(kp, cn.calls, loh.height = 0, r0=0.3,  cn.colors = c("#0000FF", "lightblue", "white", "pink", "#FFAAAA"))
cn.cols <- getCopyNumberColors(colors = c("lightblue", "white", "pink"))
plotCopyNumberSummary(kp, cn.calls, r1=0.25, direction = "out", gain.color = "red", loss.color = "blue")
legend("top", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols))
dev.off()

