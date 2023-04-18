# ENGEP
R package supporting the paper **"ENGEP: An ensemble learning tool for spatially
unmeasured genes expression prediction"**. 

ENGEP combins multiple prediction results from different reference datasets and prediction methods using a weighted average ensemble strategy to predict expression levels of spatially unmeasured genes. ENGEP mainly includes two steps: (i) generating multiple base results using k-nearest-neighbor (k-NN) regression. Different reference datasets, similarity measures, and numbers of neighbors (k) are used for this step. (ii) Combining these base results into a consensus result using a weighted average ensemble strategy. In this step, weights are assigned to different reference datasets to take into account their predictive power. 
![image](https://github.com/Zhangxf-ccnu/ENGEP/blob/master/docs/Figure1.jpg)

## Installation

 ``` buildoutcfg
 install.packages("devtools")
 devtools::install_github("Zhangxf-ccnu/ENGEP")
 ```
Note that package ‘propr’ was removed from the [CRAN](https://cran.r-project.org/web/packages/propr/index.html) repository, we advise that the user should install [propr](https://github.com/tpq/propr) before install ENGEP. 
 ``` buildoutcfg
 devtools::install_github("tpq/propr")
 ```
 
## Tutorial

A tutorial with examples of the usage of `ENGEP` is available at: [ENGEP-examples.html](https://github.com/Zhangxf-ccnu/ENGEP-examples)

## Contact

Please do not hesitate to contact Miss **Yang Shi-Tong** [styang@mails.ccnu.edu.cn](styang@mails.ccnu.edu.cn) or Dr. **Xiao-Fei Zhang** [zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.
