# Customer Analytics at ShopCo

Individual project, Modern Computational Methods (Antoine Foucart).

- `data/project_data_Y0PV0T.csv`: personal data file (560 customers)
- `analysis.R`: single script reproducing every number and figure (Tasks 1 to 6)
- `figures/`: figures written by the script
- `output.txt`: full console output of the script
- `report/report.md` → `report/report.pdf`: the report (management summary + one section per task)

Run from the repository root:

```
Rscript analysis.R > output.txt
cd report && pandoc report.md -o report.pdf --pdf-engine=pdflatex
```

Packages: mice, survival, cluster, ca, pROC (plus splines, shipped with R).
