library(sleuth)
setwd("~/Downloads/cuffdiff2_data_kallisto_results")
sample_id <- dir(file.path("results"))
kal_dirs <- sapply(sample_id, function(id) file.path("results", id, "kallisto"))
s2c <- read.table(file.path("hiseq_info.txt"), header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

## Can instead aggregate on the gene level, rather than transcript level:
# so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
#  aggregation_column = 'ens_gene')

so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionscramble')

## Or can alternatively do a likelihood ratio test:
# so <- sleuth_fit(so, ~1, 'reduced')
# so <- sleuth_lrt(so, 'reduced', 'full')
# results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

results_table <- sleuth_results(so, test='conditionscramble', test_type = 'wald')
results_table[grep("HOXA1",results_table[,13]),]

sleuth_live(so)
