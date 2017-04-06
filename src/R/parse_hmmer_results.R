library(data.table)


# get the data, need to parse headers separately. awesome
fread("grep '^[^#]' output/hmm/seq_hits.tab")
