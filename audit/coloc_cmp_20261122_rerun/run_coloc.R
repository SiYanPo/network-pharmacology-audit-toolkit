# 由 np_coloc_compare.py 生成：在装有 R + coloc 的机器上运行
# install.packages("coloc")
library(coloc)
d <- read.delim("C:/Users/huawei/AppData/Local/Temp/np_coloc_run_te4w8mxz/loci.tsv", stringsAsFactors = FALSE)
res <- lapply(split(d, d$locus), function(x) {
  e <- list(snp = x$snp, beta = x$beta_x, varbeta = x$se_x^2, type = "quant",
            N = 1000, sdY = 1)
  o <- list(snp = x$snp, beta = x$beta_y, varbeta = x$se_y^2, type = "quant",
            N = 5000, sdY = 1)
  tryCatch(coloc.abf(e, o)$summary, error = function(err) NULL)
})
out <- do.call(rbind, lapply(names(res), function(k) {
  if (is.null(res[[k]])) return(data.frame(locus = k, PP.H4 = NA))
  data.frame(locus = k, PP.H0 = res[[k]][["PP.H0.abf"]], PP.H1 = res[[k]][["PP.H1.abf"]],
             PP.H2 = res[[k]][["PP.H2.abf"]], PP.H3 = res[[k]][["PP.H3.abf"]],
             PP.H4 = res[[k]][["PP.H4.abf"]])
}))
write.csv(out, "C:/Users/huawei/AppData/Local/Temp/np_coloc_run_te4w8mxz/r_coloc_pp.csv", row.names = FALSE)
cat("R coloc done ->", "C:/Users/huawei/AppData/Local/Temp/np_coloc_run_te4w8mxz/r_coloc_pp.csv", "\n")
