args <- commandArgs(TRUE); code <- args[1]; python <- args[2]; out <- args[3]
suppressPackageStartupMessages(library(Seurat))
test_lib <- Sys.getenv("DEGAS_TEST_LIBRARY", "")
if (nzchar(test_lib)) .libPaths(c(test_lib, .libPaths()))
library(DEGASv2)
stopifnot(file.exists(system.file('python','scanpy_markers.py',package='DEGASv2')))
set.seed(42)
x <- matrix(rpois(80*180,1),80,180,dimnames=list(paste0('g',1:80),paste0('c',1:180)))
x[1:10,1:60] <- x[1:10,1:60]+8
x[11:20,61:120] <- x[11:20,61:120]+8
x[21:30,121:180] <- x[21:30,121:180]+8
x[80,] <- 0
obj <- NormalizeData(CreateSeuratObject(counts=Matrix::Matrix(x,sparse=TRUE)),verbose=FALSE)
Idents(obj) <- rep(c('0','1','2'),each=60)
a <- FindAllMarkers(obj,only.pos=FALSE,min.pct=.25,logfc.threshold=.25,verbose=FALSE)
b <- FindAllMarkersScanpy(obj,directory=out,python=python)
z <- merge(a,b,by=c('gene','cluster'),suffixes=c('.r','.py'))
stopifnot(nrow(z)==nrow(a),nrow(z)==nrow(b),
  max(abs(z$avg_log2FC.r-z$avg_log2FC.py))<1e-12,
  max(abs(z$pct.1.r-z$pct.1.py))<1e-12,
  max(abs(z$pct.2.r-z$pct.2.py))<1e-12,
  any(b$avg_log2FC<0),all(b$p_val_adj==pmin(b$p_val*nrow(obj),1)),
  !('g80' %in% b$gene))
# P-values need not be identical: record numerical discrepancies explicitly.
write.csv(z,file.path(out,'comparison.csv'),row.names=FALSE)
jsonlite::write_json(list(status='passed',rows=nrow(z),
  max_absolute_p_difference=max(abs(z$p_val.r-z$p_val.py)),
  same_top_20=identical(as.character(head(a[order(a$p_val_adj,-a$avg_log2FC),]$gene,20)),
                      as.character(head(b[order(b$p_val_adj,-b$avg_log2FC),]$gene,20)))),
  file.path(out,'test.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
