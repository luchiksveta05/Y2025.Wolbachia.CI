#' Build fasta file for cif genes
#' @description Build fasta file for cif genes
#' @param i row in the pseudo bed file
#' @param gene name of the gene(s) as indicated in the pseudo bed file
#'
#' @return A vector of positions
#'
#' @export
#'
#' @examples
#'

# Define the infection status function
extract_gene_fasta <- function(i, gene, ref) {

    contig <- gene$contig3[i]
    start  <- gene$start[i]
    end    <- gene$end[i]
    strand <- gene$strand[i]

    # Find the corresponding reference contig
    idx <- which(names(ref) == contig)

    if (length(idx) != 1) {
      stop("Could not uniquely find contig: ", contig)
    }

    seq <- XVector::subseq(
      ref[[idx]],
      start = start,
      end = end
    )

    # Reverse complement if gene is on minus strand
    if (strand == "-") {
      seq <- Biostrings::reverseComplement(seq)
    }

    seq
  }



