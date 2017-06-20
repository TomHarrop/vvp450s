#!/usr/bin/env Rscript

library(GenomicRanges)
library(GenomicFeatures)
library(data.table)

# load the maker gff
vvul_gff_file = "data/maker_filtered.gff3"
vvul_gff3 <- rtracklayer::import.gff3(vvul_gff_file)

# convert to txdb
txdb <- makeTxDbFromGRanges(vvul_gff3)

# dom_hits headers from manually reading the file
headers <- c("target_name",
             "accession",
             "tlen",
             "query_name",
             "accession",
             "qlen",
             "E-value",
             "score",
             "bias",
             "#",
             "of",
             "c-Evalue",
             "i-Evalue",
             "score",
             "bias",
             "hmm_from",
             "hmm_to",
             "ali_from",
             "ali_to",
             "env_from",
             "env_to",
             "acc",
             "description_of_target")
# get the data, need to parse headers separately. awesome
dom_hits <- fread("grep '^[^#]' output/hmm/dom_hits.tab")
dom_hits[, V23_new := paste(V23, V24, V25, V26)]
dom_hits[, c("V23", "V24", "V25", "V26") := NULL]
names(dom_hits) <- headers

# add a hit_id column to merge later
dom_hits[, hit_id := paste0("hit", 1:dim(dom_hits)[1])]

# convert the protein coordinates to cds coordinates (i.e. multiply by three)
hits_cds <- dom_hits[, .(
    ali_from_cds = ali_from * 3,
    ali_to_cds = ali_to * 3),
    by = .(hit_id, target_name)]

# convert the hits to GenomicRanges
hits_gr <- makeGRangesFromDataFrame(hits_cds,
                                    keep.extra.columns = TRUE,
                                    ignore.strand = TRUE,
                                    seqnames.field = "target_name",
                                    start.field = "ali_from_cds",
                                    end.field = "ali_to_cds")

hits_gr[seqnames(hits_gr) == "Vvul006892-mRNA-2"]

# map the exon coordinates to genomic coordinates with
# GenomicFeatures::mapFromTranscripts
# group transcripts by...
cds <- cdsBy(txdb, by = "tx", use.names = TRUE)

mapped_hits <- mapFromTranscripts(hits_gr,
                                  transcripts = cds)
 
# need to subtract intron regions here, maybe overlap with exons?

# convert mapped_hits to data.table
mapped_hits_dt <- as.data.table(mapped_hits)

# get hit_id from hits_gr
GetHitIDFromHitNumber <- function(hit_number, hits_gr_object) {
    as.character((hits_gr_object[hit_number])$hit_id)
}
mapped_hits_dt[, hit_id := GetHitIDFromHitNumber(xHits, hits_gr)]

# add source, scores etc
hits_with_name <- merge(
    mapped_hits_dt,
    dom_hits[, .(hit_id,
                 target_name,
                 query_name,
                 `c-Evalue`,
                 `i-Evalue`,
                 score,
                 bias)],
    by = "hit_id")

# convert back to genomic ranges
hits_with_name[, c("hit_id", "xHits", "transcriptsHits") := NULL]
setnames(hits_with_name, "target_name", "Name")
hits_with_name[, type := "match"]
hits_with_name[, source := "HMMER 3.1b2"]
output_gr <- makeGRangesFromDataFrame(hits_with_name,
                         keep.extra.columns = TRUE,
                         seqnames.field = "seqnames",
                         start.field = "start",
                         end.field = "end",
                         strand.field="strand")

# write as gff
rtracklayer::export.gff3(output_gr, "output/parsed_results/hmmer.gff3")

