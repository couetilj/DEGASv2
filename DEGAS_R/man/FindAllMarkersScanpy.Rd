\name{FindAllMarkersScanpy}
\alias{FindAllMarkersScanpy}
\title{Find cluster markers with Scanpy Wilcoxon}
\usage{
FindAllMarkersScanpy(object, only.pos = FALSE, min.pct = .25,
  logfc.threshold = .25, return.thresh = .01,
  directory = tempfile("degas-scanpy-"),
  python = Sys.getenv("DEGAS_SCANPY_PYTHON", "python3"),
  script = system.file("python", "scanpy_markers.py", package = "DEGASv2"))
}
\arguments{
\item{object}{Normalized Seurat object with cluster identities.}
\item{only.pos}{Return only positive markers.}
\item{min.pct}{Minimum detection fraction in either group.}
\item{logfc.threshold}{Minimum absolute log2 fold change.}
\item{return.thresh}{Raw p-value threshold.}
\item{directory}{Exchange/output directory. Files are retained for inspection.}
\item{python}{Python executable with the pinned Scanpy environment.}
\item{script}{Path to the installed Python marker script.}
}
\description{Uses Seurat fold changes and detection fractions with Scanpy's
asymptotic tie-corrected Wilcoxon tests against the rest. Bonferroni correction
uses all assay genes. P-values need not be identical to Seurat. Install the
optional Python dependencies from python/requirements-scanpy.txt in the package.}
\value{A data frame with Seurat-compatible marker columns.}
