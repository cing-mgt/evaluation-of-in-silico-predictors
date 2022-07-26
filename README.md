# Evaluation of in silico predictors on short nucleotide variants in HBA1, HBA2 and HBB associated with haemoglobinopathies

This is a public repo with input files and all scripts and files produced for 

>Tamana, S., Xenophontos, M., Minaidou, A., Stephanou, C., Harteveld, C. L., Bento, C., Traeger-Synodinos, J., Fylaktou, I., Yasin, N. M., Abdul Hamid, F. S., Esa, E., Halim-Fikri, H., Zilfalil, B. A., Kakouri, A. C., ClinGen Hemoglobinopathy VCEP, Kleanthous, M., & Kountouris, P. (2022). Evaluation of in silico predictors on short nucleotide variants in HBA1, HBA2 and HBB associated with haemoglobinopathies. eLife, 11. https://doi.org/10.7554/eLife.79713

The initial step of the analyis (binary classification) is performed in `evaluate-performance.Rmd`, while the evaluation fo tools using different pathogenic and benign thresholds is performed in `refine-tool-thresholds.Rmd`. Both these scripts make use of custom functions defined in `source/metrics.r`. N.B.: We also provide `metrics-for-R3.6.r` to be able to run the scripts for an older version of R and the `epiR` package.

All the figures in the paper, along with the files containing the formatted and/or transformed data for plotting (found under `plots/plot-source-files`) are produced in `make-plots.Rmd`
