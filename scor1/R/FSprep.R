#' A function to prepare data for FastSpar
#'
#' Creates a data.frame from the original \code{phyloseq} object
#' with the relevant conditions satisfied.
#'
#' @param physeq (Required). A \code{phyloseq} object containing abundance table,
#'        taxonomic assignment, sample data and/or phylogenetic tree if available.
#'
#' @param condition (Required). A condition that is compatible with \code{phyloseq} object.
#'        A metadata header == metadata state. E.g. \code{Health==healthy}
#' @param prev. Prevalence of an OTU as a fraction of the total number of samples. Default = 0.1
#' @return A \code{data.frame} object with `#OTU ID` inserted in the first cell.
#' @import phyloseq
#' @importFrom tibble rownames_to_column
#' @export FSprep
#' @author Abhishake Lahiri \email{abhishake.lahiri@gmail.com}
#' @examples
#' data(GlobalPatterns)
#' sample_data(GlobalPatterns)
#' GPfeces <- FSprep(GlobalPatterns, SampleType==feces,prev=0.1)
#' write.table(GPfeces, sep ="\t", file = "Can.txt", quote = FALSE,row.names = FALSE)
FSprep <- function(physeq,condition,prev=0.1){

  d2sH <-subset_samples(physeq, condition)

  mdt = fast_melt(d2sH)
  prevdt = mdt[, list(Prevalence = sum(count>0),
                      TotalCounts = sum(count)),
               by = taxaID]
  mytaxa = prevdt[(Prevalence >=prev*nsamples(d2sH)),1]

  d2sHf1 <- prune_taxa(mytaxa$taxaID,d2sH)
  d2sHf1
  outH <- as.data.frame(otu_table(d2sHf1)) %>%
    rownames_to_column("OTU_ID")
  return(outH)
}


