# Updates

- add sQTL cis permutation mudule.
__commit: 1c8c3be549842cd3ab5b11d19764872f96b43827__
__file: l3_permutation.cpp; function: permu_sqtl__
add new source file l3_permutation.cpp and related cis_learn_beta.cpp.  
add --permu and --permu-times flags.
example:
osca --sqtl --permu --permu-time 1000 --bfile filename --befile finename --qcover filename --bed filename --maf 0.01 --call 0.85 --cis

For snp of a gene, the snp have smallest p value is best_snp(nominal).

Then shuffle phenotype-sample data at sample dimension.
then calculate p value of snp belong this gene, and keep the smallest p, this do a permutation. Every permutation will generate a smallest p.

Afer permutation, using function learnBetaParameters(from QTLtools ) and permuted data to ajust the beta_ml1 and beta_ml2.Use function  pbeta(from R) and ajdusted beta_ml1, beta_ml2 to recalculate p_best(nominal).

The result file contain following fields:
beta_ml1: ajusted beta_ml1 using permuted data.
beta_ml2: ajusted beta_ml2 using permuted data.
pemp: number of permuted data which less equal than p_best(before ajusted) devied by permutation times. 
pbml: p_best after recalcute.


- add --trans and trans-wind to ssQTL
__commit: 9471c0bef0f41fdf5e4a667b98a65f309e1fda52__
__file: l3_vqtl.cpp; function: ssQTL__
add trans funtion to ssQTL when using trans in ssQTL, snp not at same chromosome of probe/gene or out of trans rang is collected for later analysis. each probe contain multi transcripts, and union the snp of each transcript.


- add --trans-meta flag to function meta
__commit: 442dbf80d53d441f1ba180a5dd1b822c5ce53619__
__commit: db572d42b9a9bdaca4f8c36026d1c94201c5ca95__

__file: l3_smr.cpp; function: meta__
If use this flag, when extract snp from multi besd files, 
will choose snp by trans method. and this flag will cauuse the output besd in DENSE format too.


- add --trans and --trans-wind to sQTL
__commit: 5b2494472374fd9f887c5ff10088fd9879e9db46__
__file: l3_vqtl.cpp; function: sQTL__
add --trans and --trans-wind to function sQTL, when --trans flag
was used, the genotype will collect by trans condition: located in different 
chromsome or located at same chromsome but out of trans_window. The default value
of --trans-wind is 5000kbp.


- add --use-top to function sQTL and ssQTL
__commit: 932644e2d070f5ea5afcee9c21221a61ae683647__
__commit: 9d8edd820615845dfdc68e2dee189cedc77a3a81__
__commit: 26863223cb1d877afd25a11156c9cea89246d9ac__
__commit: 4b786bbc452adaad6ca9240926a2f1fdbf894eb7__
__file: l3_vqtl.cpp; function: sQTL; function: ssQTL__


- modification
__commit: ca4c782d020402ad7838f1038096d39890677566__
__file: l3_vqtl.cpp; funciton: ssQTL__
Using `SelfAdjointEigenSolver<MatrixXd> es(corr_dev, EigenvaluesOnly);` instead of `SelfAdjointEigenSolver<MatrixXd> es(corr_dev);` to improve efficency.

- bug fix 
__commit: ca4c782d020402ad7838f1038096d39890677566__
__file: l3_vqtl.cpp; function: ssQTL__
Filter `cor_null` if the row contain value 1 except diagnal. 

- modification 
__commit: 8235a21739ac85946929935f8af6174ee01e8006__
__file: l3_vqtl.cpp; function: sQTL__
Using `SelfAdjointEigenSolver<MatrixXd> es; es.compute(corr_dev, EigenvaluesOnly);` instead of `SelfAdjointEigenSolver<MatrixXd> es(corr_dev);` to improve efficency.

- bug fix
__commit: e26ddebbc4beca2772708be13b93a08064e97932__
__commit: b796a7c01e01319fdbf9f90359c30033d41d98dc__
__commit: b1e516c8c893a7eba3855050fb22c3dd14859485__
__file: l3_vqtl.cpp; function: sQTL__
Check `cor_null` metrix, and if a row contain value of none diagnal is eqtal to 1. remove this row and correspoding column. Remove those rows form `trpv`.