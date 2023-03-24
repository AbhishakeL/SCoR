#' Main function for ScoR
#'
#' Creates a \code{list} output which includes an igraph object,
#'         a module count, r threshold value for SparCC correlation,
#'         r threshold for distance correlation and
#'         a final correlation table in long format.
#'
#' @param otu_tab (Required). File path for the OTU table. Remove `#` from the first cell.
#'
#' @param sparCor (Required). File path for SparCC correlation matrix. Remove `#` from the first cell.
#'
#' @param sparP (Required). File path for SparCC p-value matrix. Remove `#` from the first cell.
#'
#' @param zc. Fraction of zero-counts allowed for each pair of OTUs. Default = 0.1.
#'            OTUs which have more zeros than zc will be excluded.
#'
#' @return A \code{list} object.
#' @import dplyr
#' @importFrom energy dcor dcov.test
#' @importFrom igraph graph_from_data_frame E delete.vertices cluster_leading_eigen set_vertex_attr degree
#' @importFrom RMThreshold rm.get.threshold
#' @importFrom reshape2 melt
#' @importFrom corrr colpair_map
#' @export scor
#' @author Abhishake Lahiri \email{abhishake.lahiri@gmail.com}
#' @examples
#' data(enterotype)
#' sample_data(enterotype)
#' EThealthy <- FSprep(enterotype, ClinicalStatus==healthy, prev = 0.2)
#' write.table(EThealthy, sep ="\t", file = "Healthy.txt", quote = FALSE,row.names = FALSE)
#' #RUN FastSpar to generate SparCC correlation and p-value.
#' sr <- scor("Healthy.txt","H_cor.txt","H_p.txt",0.2)
#'

scor <- function(otu_tab,sparCor,sparP,zc=0.2){
  s2ctrl <- read.table(otu_tab, header = TRUE, sep = "\t",
                       row.names = 1, check.names = FALSE)
  s2ctrlSc <- read.table(sparCor, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
  s2ctrlSp <- read.table(sparP, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
  ###Sparcc thresholding####
  s2ctrlScsign <- s2ctrlSc
  thSpar <- rm.get.threshold(as.matrix(abs(s2ctrlSc)),
                             save.fit = FALSE)
  s2ctrlSc[abs(s2ctrlSc) < thSpar$sse.chosen] <- 0.00
  s2ctrlSp.yes <- s2ctrlSp < 0.05
  spr.fi <- s2ctrlSc*s2ctrlSp.yes

  ###dcorr thresholding####
  s2ctrlDcorr <- colpair_map(as.data.frame(t(s2ctrl)), dcor,.diagonal = 1)
  s2ctrlDcorr <- s2ctrlDcorr[-c(1)]
  rownames(s2ctrlDcorr) <- colnames(s2ctrlDcorr)
  thDcor <- rm.get.threshold(as.matrix(s2ctrlDcorr),
                             save.fit = FALSE)
  s2ctrlDcorr[s2ctrlDcorr < thDcor$sse.chosen] <- 0.00
  s2ctrlDp <- colpair_map(as.data.frame(t(s2ctrl)), distPtest,.diagonal = 1)
  s2ctrlDp <- s2ctrlDp[-c(1)]
  rownames(s2ctrlDp) <- colnames(s2ctrlDp)
  s2ctrlDp.yes <- s2ctrlDp < 0.05
  s2ctrlD.fi <- s2ctrlDcorr * s2ctrlDp.yes
  rownames(s2ctrlD.fi) <- colnames(s2ctrlD.fi)

  ###Zero-pairs####
  s2ctrlZc <- colpair_map(as.data.frame(t(s2ctrl)), zeroPair,.diagonal = NA)
  s2ctrlZc <- s2ctrlZc[-c(1)]
  rownames(s2ctrlZc) <- colnames(s2ctrlZc)

  ###Traingles####
  df_s = (get_upper_tri(as.matrix(spr.fi))) %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 3))
  df_d = (get_upper_tri(as.matrix(s2ctrlD.fi))) %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 3))
  df_z = (get_upper_tri(as.matrix(s2ctrlZc))) %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    mutate(presence = ifelse(value<zc*ncol(s2ctrl),1,0))
  df_sd <- dplyr::full_join(df_s,df_d,by=c("Var1","Var2"))
  colnames(df_sd) <- c("Var1","Var2","sparcc","dcorr")

  df_sdz <- dplyr::full_join(df_sd,df_z,by=c("Var1","Var2"))
  df_sSign = (get_upper_tri(as.matrix(s2ctrlScsign))) %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    mutate(Sign = sign(value)) %>%
    select(-value)
  df_sdzS <- dplyr::full_join(df_sdz,df_sSign,by=c("Var1","Var2"))
  df_fi <- df_sdzS %>%
    filter(presence == 1) %>%
    mutate(R = ifelse(abs(sparcc)>dcorr,abs(sparcc),dcorr)) %>%
    mutate(Corl = R*Sign) %>%
    select("Var1","Var2","Corl")

  ###Construct Network####
  com.gr <- graph_from_data_frame(df_fi, directed = FALSE)
  E(com.gr)$weight <- df_fi$Corl
  com.up <- delete.vertices(com.gr, degree(com.gr) < 1)
  clp <- cluster_leading_eigen(com.up,
                               weights = abs(df_fi$Corl),
                               options = list(maxiter = 30000))
  com.up <- set_vertex_attr(com.up, "module",
                            value = clp$membership)

  out <- list(network = com.up,
              module = clp,
              Spar.th = thSpar,
              Dcorr.th = thDcor,
              FinalCor = df_fi)
  return(out)

}

zeroPair <- function(x,y){
  zero = sum(ifelse((ifelse(x>0,1,0) == ifelse(y>0,1,0)),0,1))
  return(zero)
}

distPtest <- function(x,y){
  ptest <- dcov.test(x, y, R=1000)
  return(ptest$p.value)
}

get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)] = NA
  diag(cormat) = NA
  return(cormat)
}
