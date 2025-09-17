dim(dds_ASTRO_filter)
sum(keep_ASTRO)
length(keep_ASTRO)

keep <- edgeR::filterByExpr(counts_ASTRO, group = colData_ASTRO$pathology_class)
sum(keep)
length(keep)

identical(names(keep),rownames(counts_ASTRO))
