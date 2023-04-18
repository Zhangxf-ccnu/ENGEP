#' @title gene_dataliat
#' @description Function to randomly partition large reference into small reference.
#'     equal-sized sub-datasets.
#' @param spa_counts gene expression matrix of spatial dataset (gene by cell).
#' @param ref_list each element in the list is a gene expression matrix (gene by cell) of
#'     an original single reference dataset.
#' @param pre_genes an array contains names of genes to be predicted, if is NULL, ENGEP will
#'     predict the intersection of unique genes of each references. If you let pre_genes = NULL,
#'     we suggest you to use reference datasets with high variable genes.
#'
#' @return a list contains equal-sized reference sub-datasets with common genes and predicted genes,
#'     spatial dataset with common genes, and unique genes in references.
#'
gene_dataliat <- function(spa_counts,ref_list,pre_genes){
  ref_common = list()
  ref_pre = list()
  common_genes = rownames(as.matrix(spa_counts))
  a=1
  for (i in 1:length(ref_list)){
    common_genes = intersect(common_genes,rownames(ref_list[[i]]))
  }
  unique_genes_i = setdiff(rownames(ref_list[[1]]), common_genes)
  for (i in 1:length(ref_list)){
    unique_genes_ii = setdiff(rownames(ref_list[[i]]), common_genes)
    unique_genes = intersect (unique_genes_i,unique_genes_ii)
    unique_genes_i = unique_genes
  }
  if (is.null(pre_genes)){
    pre_genes = unique_genes
  }

  print(paste("Predict", length(pre_genes), "genes", sep = ""))

  for (i in 1:length(ref_list)){
    n = ncol(as.matrix(ref_list[[i]]))
    K = ceiling(n/8000)
    sample_id = sample(rep(1:K, length.out=n))
    for(j in 1:K){
      sample_use = sample_id==j
      ref_common[[a]] = ref_list[[i]][,sample_use][common_genes,]
      ref_pre[[a]] = ref_list[[i]][,sample_use][pre_genes,]
      a=a+1
    }
  }
  return(list(ref_c = ref_common,
              ref_p = ref_pre,
              qur_c = spa_counts[common_genes,],
              uni_genes = unique_genes))
}

#' @title engep_predict
#' @description Function to get expression levels of unmeasured genes predicted by ENGEP.
#' @param spa_counts gene expression matrix of spatial dataset (gene by cell).
#' @param ref_list each element in the list is a gene expression matrix (gene by cell) of an
#'     original single reference dataset.
#' @param pre_genes an array contains names of genes to be predicted, if is NULL, ENGEP will
#'     predict the intersection of unique genes of each references.
#' @param nCpus number of (maximum) cores to use for parallel execution, default to 6.
#' @param simi_list a vector indicates the ten similarity measures that are used, note that the
#'     similarity measures used should be chosen from the ten similarity measures mentioned in the
#'     document. Default is c("pearson",  "spearman","cosine", "jaccard","weighted_rank",
#'     "manhattan","canberra","euclidean", "phi_s","rho_p")
#' @param parallel a logical value to indicate if the process should be run parallelly in multiple threads,
#'     default to TURE.
#' @param k_list a list contains different values of $k$ (number of neighbors in knn),
#'     default is (20,30,40,50).
#' @param get_baes a logical value to indicate whethe to return base results or not,
#'     default is FAISE.
#' @param get_weight a logical value to indicate whethe to return weights of different references,
#'     default is FAISE.
#' @export
#' @return the predicted expression levels of unmeasured genes.
#'

engep_predict <- function(spa_counts,ref_list,pre_genes,nCpus=6,simi_list= NULL,
                          parallel=TRUE,k_list = NULL,get_baes=FALSE,get_weight=FALSE){

  if(is.null(simi_list)){
    simi_list = c("pearson",  "spearman","cosine","jaccard","weighted_rank","manhattan",
                  "canberra","euclidean","phi_s","rho_p")
  }

  if (length(setdiff(simi_list, c("pearson", "spearman", "cosine", "jaccard", "weighted_rank",
                                  "manhattan", "canberra", "euclidean", "phi_s","rho_p"))) > 0) {
    stop("Error: Similarity measures should be chosen from the ten measures mentioned in the document")
  }

  if(is.null(k_list)){
    k_list = c(20,30,40,50)
  }

  message("Partition large reference")
  data_list = gene_dataliat(spa_counts,ref_list,pre_genes)
  rm(ref_list)
  gc()

  if (is.null(pre_genes)){
    pre_genes = data_list$uni_genes
  }

  k_r = length(data_list$ref_c)
  if (parallel==TRUE){
    message("Run ENGEP parallelly to get base results")
    result_single = parallel::mclapply(1:k_r, imp_new, data_list$qur_c,data_list$ref_c,simi_list,k_list,data_list$ref_p,mc.cores=6)
  }else{
    message("Run ENGEP in one thread to get base results")
    result_single =lapply(1:k_r, imp_new, data_list$qur_c,data_list$ref_c,simi_list,k_list,data_list$ref_p)
  }
  message("Combine base results")
  engep_result = ensemble_result(result_single,get_weight)
  if (get_baes == TRUE){
    return(list("ensemble"=engep_result,"base result"=result_single))
  }else {
    return(engep_result)
  }
}

#' @title imp_new
#' @description Function to generate base results by using different reference datasets, different
#'     similarity measures and different values of k.
#' @param i each time one reference expression matrix in the list is used to predict the expression.
#' @param spcom expression matrix of spatial data with common genes.
#' @param sccomlist a list contains expression matrices of different reference datasets with common genes.
#' @param similist  a vector indicates the ten similarity measure that are used.
#' @param k.list a list contains different values of $k$ (number of neighbors in knn), default is (20,30,40,50).
#' @param sc_implist a list contains expression matrices of different reference datasets with genes
#'     to be predicted.
#'
#' @return the predicted expression of unmeasured genes and the R-squared sore between the predicted expression
#'     and the measured expression of training genes representing predictive power of reference.
#' @export
#'
#' @examples
imp_new <- function(i,spcom,sccomlist,similist,k.list,sc_implist){
  sccom = sccomlist[[i]]
  sc_imp = sc_implist[[i]]
  sccom <- Seurat::CreateSeuratObject(sccom,project = "sc_use")
  sccom <- Seurat::NormalizeData(object = sccom,normalization.method = "LogNormalize",verbose = FALSE)

  spcom <- Seurat::CreateSeuratObject(counts = spcom, project = 'sp_use')
  spcom <- Seurat::NormalizeData(object = spcom,normalization.method = "LogNormalize",verbose = FALSE)

  commongene = rownames(spcom@assays$RNA@data)
  imp_gene = matrix(0,nrow = ncol(spcom@assays$RNA@data),ncol = length(rownames(sc_imp)))
  colnames(imp_gene) <- rownames(sc_imp)
  rownames(imp_gene) <- colnames(spcom@assays$RNA@data)

  imp_gene_t = matrix(0,nrow = ncol(spcom@assays$RNA@data),ncol = length(commongene))
  colnames(imp_gene_t) <- rownames(commongene)
  rownames(imp_gene_t) <- colnames(spcom@assays$RNA@data)

  gene.train = commongene

  for (s in 1:length(similist)){
    tt = compute_simi(spcom[gene.train,],sccom[gene.train,],similarity = similist[s])

    for (k in 1:length(k.list)){
      topKNN_idx <- apply(tt, 1, function(x)
        colnames(tt)[order(x, decreasing = TRUE)][seq_len(k.list[k])])

      topKNN_idx.vec = as.vector(topKNN_idx)
      col.idx = match(topKNN_idx.vec,colnames(tt))
      row.idx = rep(1:dim(tt)[1],each = k.list[k])

      ad_matrix <- Matrix::sparseMatrix(i = row.idx, j = col.idx, x = 1,
                                        dims = c(nrow(tt), ncol(tt)),
                                        dimnames = list(rownames(tt), colnames(tt)))

      w_matrix = ad_matrix * tt
      w_matrix = as.matrix(w_matrix)
      w_matrix = pmax(w_matrix,0)

      w_scale = w_matrix /  (rowSums(w_matrix)+1e-16 )
      rm(w_matrix)
      gc()
      imp_genes = w_scale %*% t(as.matrix(sc_imp))
      imp_train = w_scale %*% t(as.matrix(sccom@assays$RNA@data[gene.train,]))

      imp_gene = imp_gene +  imp_genes
      imp_gene_t = imp_gene_t + imp_train

    }
  }
  imp_gene = imp_gene/(length(k.list)*length(similist))
  imp_gene_t = imp_gene_t/(length(k.list)*length(similist))
  r2_single = caret::R2(as.vector(imp_gene_t), as.vector(t(as.matrix(spcom@assays$RNA@data[gene.train,]))))
  return(list("exp"=imp_gene,"r2"=r2_single))

}

#' @title ensemble_result
#' @description Function to combine base results by a weighted average ensemble method.
#' @param result_single a list contains base results predicted by different reference datasets, different
#'     similarity measures and different values of k.
#' @param get_weight a logical value to indicate whethe to return weights of different references,
#'     default is FAISE.
#'
#' @return ensembel result.
#'
#' @examples
ensemble_result <- function(result_single,get_weight=FALSE){
  r2_vec <- rep(0,length(result_single))
  for(i in 1:length(result_single)){
    r2_vec[i]=result_single[[i]]$r2
  }

  r2_vec <-weight_r2(r2_vec)
  weight_vec_k = r2_vec/sum(r2_vec)

  fianl_re <- matrix(0,nrow = dim(result_single[[1]]$exp)[1],ncol = dim(result_single[[1]]$exp)[2])

  for (i in 1:length(result_single)){
    fianl_re <- fianl_re + result_single[[i]]$exp * weight_vec_k[i]
  }
  if(get_weight==TRUE){
    return(list("exp"=fianl_re,"weight"=weight_vec_k))
  }else{
    return(fianl_re)
  }
}

#' @title weight_r2
#' @description Function to convert the predictive power to weight, with a range of values between 0.1 and 0.9.
#' @param r2_vec R-squared sore which presents the predictive power of reference data.
#'
#' @return weight converted from R-squared sore.
#'
#' @examples
weight_r2<-function(r2_vec){
  a=(0.9-0.1)/(max(r2_vec)-min(r2_vec))
  b=((0.9+0.1)-a*(max(r2_vec)+min(r2_vec)))/2
  rec_n=a*r2_vec +b
  return (rec_n)
}

#' @title compute_simi
#' @description Function to caculate similarity matrix between spatial and reference data.
#'     Distance measures (such as Manhattan distance, Canberra distance, Euclidean distance,
#'     and phi_{s} proportionality measure) are transformed into similarity scores using
#'     $s = 1/(1+d)$. The codes are modified from the calculateSimilarity function scclassify.
#' @param spcom common gene expression matrix of spatial data.
#' @param sccom common gene expression matrix of reference data.
#' @param similarity a vector indicates the ten similarity measure that are used.
#' @return similarity matrix between spatial and reference data.
#'
#' @examples
compute_simi <- function(spcom,sccom,similarity) {
  similarity <- match.arg(similarity, c("pearson",  "spearman","cosine",
                                        "jaccard","weighted_rank",
                                        "manhattan","canberra","euclidean",
                                        "phi_s","rho_p"))


  if (similarity == "pearson") {
    corMat <- proxyC::simil(Matrix::t(spcom@assays$RNA@data),
                            Matrix::t(sccom@assays$RNA@data),
                            method = "correlation")

    corMat[is.na(corMat) | is.infinite(corMat)] <- -1

  } else if (similarity == "spearman") {
    corMat <- stats::cor(as.matrix(spcom@assays$RNA@data),
                         as.matrix(sccom@assays$RNA@data), method = "spearman")
    corMat[is.na(corMat) | is.infinite(corMat)] <- -1

  }else if (similarity == "jaccard") {

    corMat <- proxyC::simil(Matrix::t(spcom@assays$RNA@data),
                            Matrix::t(sccom@assays$RNA@data),
                            method = "jaccard")
    corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat,na.rm=TRUE)


  } else if (similarity == "weighted_rank") {
    corMat <- wtd_rank2(as.matrix(spcom@assays$RNA@data),
                        as.matrix(sccom@assays$RNA@data))
    corMat[is.na(corMat) | is.infinite(corMat)] <- -1

  } else if (similarity == "manhattan") {

    dist_m <- as.matrix(proxy::dist(t(as.matrix(spcom@assays$RNA@data)),
                                    t(as.matrix(sccom@assays$RNA@data)),
                                    method = "Manhattan"))

    dist_m[is.na(dist_m) | is.infinite(dist_m)] <- max(dist_m,na.rm=TRUE)

    corMat <- 1/(1+dist_m)

  } else if (similarity == "cosine") {
    corMat <- proxyC::simil(Matrix::t(spcom@assays$RNA@data),
                            Matrix::t(sccom@assays$RNA@data),
                            method = "cosine")
    corMat[is.na(corMat) | is.infinite(corMat)] <- -1

  } else if (similarity == "phi_s") {
    mat <- cbind(as.matrix(spcom@assays$RNA@counts),
                 as.matrix(sccom@assays$RNA@counts))
    dis =  propr::phis(mat, select = colnames(mat))@matrix
    dist_m <- dis[1:dim(spcom@assays$RNA@counts)[2],(dim(spcom@assays$RNA@counts)[2]+1):dim(mat)[2]]
    dist_m[is.na(dist_m) | is.infinite(dist_m)] <- max(dist_m,na.rm=TRUE)
    corMat <- 1/(1+dist_m)

  } else if (similarity == "rho_p") {
    mat <- cbind(as.matrix(spcom@assays$RNA@counts),
                 as.matrix(sccom@assays$RNA@counts))
    cor = propr::perb(mat, select = colnames(mat))@matrix
    corMat <- cor[1:dim(spcom@assays$RNA@counts)[2],(dim(spcom@assays$RNA@counts)[2]+1):dim(mat)[2]]
    corMat[is.na(corMat) | is.infinite(corMat)] <- -1


  } else if (similarity == "canberra") {
    dist_m <- as.matrix(proxy::dist(t(as.matrix(spcom@assays$RNA@data)),
                                    t(as.matrix(sccom@assays$RNA@data)),
                                    method = "canberra"))

    dist_m[is.na(dist_m) | is.infinite(dist_m)] <- max(dist_m,na.rm=TRUE)
    corMat <- 1/(1+dist_m)

  } else if (similarity == "euclidean") {

    dist_m <- as.matrix(proxy::dist(t(as.matrix(spcom@assays$RNA@data)),
                                    t(as.matrix(sccom@assays$RNA@data)),
                                    method = "euclidean"))
    dist_m[is.na(dist_m) | is.infinite(dist_m)] <- max(dist_m,na.rm=TRUE)
    corMat <- 1/(1+dist_m)

  }

  return(corMat)
}

wtd_rank2 <- function(mat1, mat2) {

  ranks1 <- apply(mat1, 2, function(x) rank(-x, ties.method = "average"))
  n1 <- nrow(mat1)
  reciprocals1 <- 1 / seq_len(n1)
  savage1 <- vapply(seq_len(n1),
                    function(i) sum(reciprocals1[i:n1]),
                    numeric(1L))
  savages1 <- ranks1
  savages1[] <- savage1[ranks1]


  ranks2 <- apply(mat2, 2, function(x) rank(-x, ties.method = "average"))
  n2 <- nrow(mat2)
  reciprocals2 <- 1 / seq_len(n2)
  savage2 <- vapply(seq_len(n2),
                    function(i) sum(reciprocals2[i:n2]),
                    numeric(1L))
  savages2 <- ranks2
  savages2[] <- savage2[ranks2]
  cor <- stats::cor(savages1, savages2)
  return(cor)
}

