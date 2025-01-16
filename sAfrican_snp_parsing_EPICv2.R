library(GenomicRanges)
library(VariantAnnotation)

options(stringsAsFactors=FALSE)

project_dir <- file.path("/scratch/chrY_varaints_in_PCa/Jenna/snp_parsing/batch2")
manifest <- file.path(project_dir, "EPIC-8v2-0_A2_no_chr0M_edit.csv.gz") #hg38
vcf <- "/scratch/chrY_varaints_in_PCa/Jenna/germline_merge/batch2/epic_germline_b2_99_nochr.vcf.gz"
report.dir <- "/scratch/chrY_varaints_in_PCa/Jenna/snp_parsing/batch2/reports"



#-------------------------------------------------------------------------------

### 1. PARSE VCF FILES ###

# Parsing VCF file
parse_vcf_results <- function(x) {
  rv <- NULL
  if (length(x$rowRanges) > 0) {
    # Variant coordinates
    temp <- as.data.frame(x$rowRanges, row.names = 1:length(x$rowRanges))
    #temp <- as.data.frame(x$rowRanges)
    rv <- temp[, c("seqnames", "start", "end")]
    rv$id <- names(x$rowRanges)
    # Variant confidence
    rv$qual <- x$QUAL
    rv$filter <- x$FILTER
    # Variant type
    #rv$type <- x$INFO$VT
    # Reference and alternative alleles
    rv$ref <- as.character(x$REF)
    rv$alt <- sapply(x$ALT, function(y) paste(y, collapse=";"))
    # Allele frequency
    rv$af <- x$INFO$AF
  }
  return(rv)     
}

# EPIC CpG probe coordinates
probe_annot <- read.csv(gzfile(manifest), skip=7, header=TRUE)
idxs_cg <- grep("^cg", probe_annot$IlmnID)
probe_annot <- probe_annot[idxs_cg, c("IlmnID", "Infinium_Design_Type", "Next_Base", "CHR", "MAPINFO", "Strand_FR")]
rownames(probe_annot) <- probe_annot$IlmnID

# Create relevant ranges: SBE, CPG, REST
type1_idxs <- probe_annot$Infinium_Design_Type=="I"
type2_idxs <- probe_annot$Infinium_Design_Type=="II"
sbe <- with(probe_annot[type1_idxs, ], GRanges(seqnames=CHR, 
                                               ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO-1, MAPINFO+2), width=1),
                                               strand=ifelse(Strand_FR=="F", "+", "-"),
                                               name=IlmnID))
cpg <- with(probe_annot, GRanges(seqnames=CHR, 
                                 ranges=IRanges(start=MAPINFO, width=2), 
                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                 name=IlmnID))
rest1 <- with(probe_annot[type1_idxs, ], GRanges(seqnames=CHR, 
                                                 ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO+2, MAPINFO-5), width=5),
                                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                                 name=IlmnID))
rest2 <- with(probe_annot[type2_idxs, ], GRanges(seqnames=CHR, 
                                                 ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO+2, MAPINFO-5), width=5),
                                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                                 name=IlmnID))
rest <- c(rest1, rest2)

#--- Retrieve variants overlapping EPIC probe regions
regions <- reduce(c(cpg, sbe, rest), ignore.strand=TRUE)
CHUNK_SIZE <- 1000
curr_idx <- 1
variants <- list()
while(curr_idx <= length(regions)) {
  print(curr_idx)
  end_idx <- min(curr_idx + CHUNK_SIZE, length(regions)+1)
  rv_vcf <- scanVcf(vcf, param=ScanVcfParam(which=regions[curr_idx:(end_idx-1)]))
  curr_variants <- lapply(rv_vcf, parse_vcf_results)
  variants <- c(variants, curr_variants)
  curr_idx <- end_idx
}
#save(variants, 
#     file=file.path(report.dir, "epic_variants_101_b2.Rdata"))
#load(file.path(report.dir, "epic_variants_101_b2.Rdata"))

CHUNK_SIZE <- 1000
curr_idx <- 2
variant_table <- variants[[1]]
while(curr_idx <= length(variants)) {
  print(curr_idx)
  end_idx <- min(curr_idx + CHUNK_SIZE, length(variants)+1)
  temp <- do.call(rbind, variants[curr_idx:(end_idx-1)])
  idxs <- sapply(temp$af, max) > 0.01
  #idxs <- temp$type != "SV" & sapply(temp$af, max) > 0.01
  #idxs <- temp$type != "SV"
  variant_table <- rbind(variant_table, temp[idxs, ])
  curr_idx <- end_idx
}
###
#colSums(is.na(variant_table))
variant_table_new <- na.omit(variant_table)
# 185 717 - 184 077 = 1640 rows omitted
###
variant_regions <- with(variant_table_new, GRanges(seqnames=seqnames, 
                                                   ranges=IRanges(start=start, end=end)))
save(variant_table_new, variant_regions, 
     file=file.path(report.dir, "epic_variants_99_b2.Rdata"))
colSums(is.na(variant_table))



#-------------------------------------------------------------------------------

### 2. CREATE SNP-AFFECTED PROBE TABLES: CpG, SBE, BODY ###

keep_cols <- c(1:6, 8:10, 13:16)
name_cols <- c("PROBE", "DESIGN_TYPE", "CHR", "MAPINFO", "STRAND", "ON_HM450",
               "VARIANT_START", "VARIANT_END", "VARIANT_ID", 
               "REF", "ALT", "AF", "MAX AF")

#--- Load EPIC CpG probe coordinates
probe_annot <- read.csv(gzfile(manifest), skip=7, header=TRUE)
idxs_cg <- grep("^cg", probe_annot$IlmnID)
probe_annot <- probe_annot[idxs_cg, c("IlmnID", "Infinium_Design_Type", "CHR", "MAPINFO", "Strand_FR", "Next_Base")]
rownames(probe_annot) <- probe_annot$IlmnID
# Create relevant ranges: SBE, CPG, REST
type1_idxs <- probe_annot$Infinium_Design_Type=="I"
type2_idxs <- probe_annot$Infinium_Design_Type=="II"
sbe <- with(probe_annot[type1_idxs, ], GRanges(seqnames=CHR, 
                                               ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO-1, MAPINFO+2), width=1),
                                               strand=ifelse(Strand_FR=="F", "+", "-"),
                                               name=IlmnID))
cpg <- with(probe_annot, GRanges(seqnames=CHR, 
                                 ranges=IRanges(start=MAPINFO, width=2), 
                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                 name=IlmnID))
rest1 <- with(probe_annot[type1_idxs, ], GRanges(seqnames=CHR, 
                                                 ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO+2, MAPINFO-5), width=5),
                                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                                 name=IlmnID))
rest2 <- with(probe_annot[type2_idxs, ], GRanges(seqnames=CHR, 
                                                 ranges=IRanges(start=ifelse(Strand_FR=="F", MAPINFO+2, MAPINFO-5), width=5),
                                                 strand=ifelse(Strand_FR=="F", "+", "-"),
                                                 name=IlmnID))
rest <- c(rest1, rest2)

#--- Load variant data
#load(file.path(report.dir, "epic_variants_101_b2.Rdata"))


## SBE Variants
#--- Find overlaps with SBE
ov <- as.matrix(findOverlaps(sbe, variant_regions, select="all"))
rv_sbe <- cbind(probe_annot[sbe$name[ov[,1]], ], variant_table_new[ov[,2], ])
#rv_sbe$type <- unlist(sapply(rv_sbe$type, paste, collapse=";"))
rv_sbe$af_max <- sapply(rv_sbe$af, max)
rv_sbe$af <- unlist(sapply(rv_sbe$af, paste, collapse=";"))
rv_sbe <- rv_sbe[, keep_cols]
colnames(rv_sbe) <- name_cols
write.csv(rv_sbe, file=file.path(report.dir, "epic_SNPs_b2_99_sbe.csv"),
          row.names=FALSE, quote=FALSE)
length(unique(rv_sbe[,1]))


## CpG Variants
#--- Find overlaps with CpG
ov <- as.matrix(findOverlaps(cpg, variant_regions, select="all"))
rv_cpg <- cbind(probe_annot[cpg$name[ov[,1]], ], variant_table_new[ov[,2], ])
#rv_cpg$type <- unlist(sapply(rv_cpg$type, paste, collapse=";"))
rv_cpg$af_max <- sapply(rv_cpg$af, max)
rv_cpg$af <- unlist(sapply(rv_cpg$af, paste, collapse=";"))
rv_cpg <- rv_cpg[, keep_cols]
colnames(rv_cpg) <- name_cols
write.csv(rv_cpg, file=file.path(report.dir, "epic_SNPs_b2_99_cpg.csv"),
          row.names=FALSE, quote=FALSE)
length(unique(rv_cpg[,1]))


## Probe Body Variants
#--- Find overlaps with rest of the probe
ov <- as.matrix(findOverlaps(rest, variant_regions, select="all"))
rv_rest <- cbind(probe_annot[rest$name[ov[,1]], ], variant_table_new[ov[,2], ])
#rv_rest$type <- unlist(sapply(rv_rest$type, paste, collapse=";"))
rv_rest$af_max <- sapply(rv_rest$af, max)
rv_rest$af <- unlist(sapply(rv_rest$af, paste, collapse=";"))
rv_rest <- rv_rest[, keep_cols]
colnames(rv_rest) <- name_cols
write.csv(rv_rest, file=file.path(report.dir, "epic_SNPs_b2_99_body_5bp.csv"),
          row.names=FALSE, col.names=TRUE, quote=FALSE)
length(unique(rv_rest[,1]))


#-------------------------------------------------------------------------------

### 3. MERGE CpG, SBE & BODY PROBES INTO ONE "BLACKLIST" ###

# probes affected by SNPs, load own lists
epic.variants1 <- read.csv("/scratch/chrY_varaints_in_PCa/Jenna/snp_parsing/batch2/reports/epic_SNPs_b2_99_sbe.csv", head = T)
epic.variants2 <- read.csv("/scratch/chrY_varaints_in_PCa/Jenna/snp_parsing/batch2/reports/epic_SNPs_b2_99_cpg.csv", head = T)
epic.variants3 <- read.csv("/scratch/chrY_varaints_in_PCa/Jenna/snp_parsing/batch2/reports/epic_SNPs_b2_99_body_5bp.csv", head = T)

# combine 3 files into 1, include PROBE, VAR_START and VAR_END columns
#epic.snp.probes <- c(as.character(epic.variants1$PROBE), as.character(epic.variants2$PROBE), as.character(epic.variants3$PROBE))
epic.snps.1 = subset(epic.variants1, select = c("PROBE", "VARIANT_START", "VARIANT_END"))
epic.snps.2 = subset(epic.variants2, select = c("PROBE", "VARIANT_START", "VARIANT_END"))
epic.snps.3 = subset(epic.variants3, select = c("PROBE", "VARIANT_START", "VARIANT_END"))
library(dplyr)
epic.snp.probes = bind_rows(epic.snps.1, epic.snps.2, epic.snps.3)

#epic.snp.probes <- unique(epic.snp.probes)
# final list of unique probes, remove duplicates by a single column (Probe IDs)
epic.snp.probes <- epic.snp.probes %>% distinct(PROBE, .keep_all = TRUE)
#epic.snp.probes <- epic.snp.probes[!duplicated(epic.snp.probes$PROBE), ] # another method

write.csv(epic.snp.probes, file=file.path(report.dir, "epic_SNPs_b2_99_all.csv"),
          row.names=T, quote=F)