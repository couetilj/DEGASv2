## Diagnostic Evidence GAuge of Single cells (DEGAS) version 2

![](figures/DEGASv2_fig1-2.png)

**Package development by:**

**Ziyu Liu**

**Travis S. Johnson (https://github.com/tsteelejohnson91)**

**Sihong Li (https://github.com/alanli97)**

**Jiahui Liu (https://github.com/ElliotLiu1997)**

## Installation

### R
```R
install.packages("devtools")
devtools::install_github("ElliotLiu1997/DEGASv2", subdir = "DEGAS_R")
```

&#9888;&#65039; Note: The folder `DEGAS_R_old` contains archived code for internal reference only.  
Do **NOT** install from that folder.

### Python
```bash
pip install git+https://github.com/ElliotLiu1997/DEGASv2.git#subdirectory=DEGAS_python
```

## **Prerequisites**

**Python packages**

pytorch

numpy

pandas

scipy

tqdm

scikit-learn

## Patient-level out-of-fold diagnostics

The Python bagging entry point can optionally hold out complete patient groups
(for example, source studies) and record calibration diagnostics. Legacy calls
are unchanged. Enable the new path with `patient_oof=True` and pass a group and
identifier for every patient:

```python
options["patient_oof"] = True
options["calibration_penalty"] = 1.0
options["calibration_bins"] = 5

cell_results = DEGAS_python.bagging_all_results(
    options,
    patient_expression,
    patient_labels,
    cell_expression,
    sc_lab_mat=cell_labels,
    pat_groups=patient_study,
    pat_ids=patient_ids,
)
```

Each submodel is trained without the complete held-out group. DEGAS writes raw
submodel predictions, one bagged out-of-fold score per patient, cross-fitted
ridge-logistic probabilities, discrimination and Brier metrics, reliability
bins, fold-specific calibration parameters, and a deployment calibrator. Apply
the deployment calibrator only after averaging a matching ensemble. Cellular
DEGAS scores are not patient disease probabilities and are not calibrated by
this procedure.

**R**

reticulate

Seurat

ggplot2

DESeq2
