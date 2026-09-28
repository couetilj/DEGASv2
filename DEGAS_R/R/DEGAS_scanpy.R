#' Find cluster markers with Scanpy Wilcoxon
#'
#' Uses Seurat fold changes and detection fractions with Scanpy p-values.
#' Requires a separate Python environment; see inst/python/requirements-scanpy.txt.
#' Clustering and the native global marker-row selection remain unchanged.
#' @param object A normalized Seurat object with cluster identities.
#' @param only.pos Return only positive markers.
#' @param min.pct Minimum detection fraction in either group.
#' @param logfc.threshold Minimum absolute log2 fold change.
#' @param return.thresh Raw p-value threshold.
#' @param directory Directory for exchange files and marker output. Files are retained.
#' @param python Python executable with the pinned Scanpy environment.
#' @param script Path to the installed Scanpy marker script.
#' @return A data frame with Seurat-compatible marker columns.
#' @export
# Preserve Seurat normalization, clustering, FoldChange, and DEGAS ranking.
FindAllMarkersScanpy <- function(object, only.pos=FALSE, min.pct=.25,
                                logfc.threshold=.25, return.thresh=.01,
                                directory = tempfile("degas-scanpy-"),
                                python = Sys.getenv("DEGAS_SCANPY_PYTHON", "python3"),
                                script = system.file("python", "scanpy_markers.py", package="DEGASv2")) {
  if (!nzchar(script) || !file.exists(script)) stop("Scanpy marker script not found; install DEGASv2 or supply script")
  if (!nzchar(Sys.which(python)) && !file.exists(python)) stop("Python executable not found: ", python)
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  x <- methods::as(SeuratObject::GetAssayData(object, assay='RNA', layer='data'), 'dgCMatrix')
  ids <- SeuratObject::Idents(object)
  stopifnot(identical(colnames(x),names(ids)))
  writeLines(rownames(x),file.path(directory,'genes.txt'))
  write.csv(data.frame(cell_id=colnames(x),cluster=as.character(ids)),
            file.path(directory,'cells.csv'),row.names=FALSE)
  writeBin(x@x,file.path(directory,'x.f64'),size=8,endian='little')
  writeBin(x@i,file.path(directory,'i.i32'),size=4,endian='little')
  writeBin(x@p,file.path(directory,'p.i32'),size=4,endian='little')
  fc <- lapply(levels(ids),function(k) {
    z <- Seurat::FoldChange(object,ident.1=k,assay='RNA',slot='data',base=2)
    z$gene <- rownames(z); z$cluster <- k; z
  })
  write.csv(do.call(rbind,fc),gzfile(file.path(directory,'fold_changes.csv.gz')),row.names=FALSE)
  jsonlite::write_json(list(min_pct=min.pct,logfc_threshold=logfc.threshold,
    return_thresh=return.thresh,only_pos=only.pos),file.path(directory,'parameters.json'),auto_unbox=TRUE)
  status <- system2(python,c(shQuote(script),shQuote(directory)))
  if (status != 0L) stop('Scanpy marker calculation failed: ',status)
  read.csv(file.path(directory,'markers.csv'),stringsAsFactors=FALSE,
           colClasses=c(cluster='character',gene='character'),check.names=FALSE)
}
