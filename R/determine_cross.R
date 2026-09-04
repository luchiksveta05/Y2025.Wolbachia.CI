#' ## Generate a consensus sequence from one Forward and one reverse-complemented Reverse read
#' @description Generate consensus sequence for forward and reverse sequenced by Sanger
#' @param forward - forward sequence
#' @param reverse_rc - reverse complement sequence
#'
#' @return A vector with infected or uninfected status according given information for the each line
#'
#' @export
#'
#' @examples
#'

## Generate a consensus sequence from one Forward and one reverse-complemented Reverse read
make_consensus <- function(forward, reverse_rc) {

  ## Align the overlapping parts of the two reads
  alignment <- Biostrings::pairwiseAlignment(forward, reverse_rc, type = "overlap", gapOpening = 5, gapExtension = 2)

  ## Extract only the aligned regions, excluding non-overlapping sequence ends
  aligned_forward <- as.character(Biostrings::alignedPattern(alignment))
  aligned_reverse <- as.character(Biostrings::alignedSubject(alignment))

  ## Return an empty result if no overlap is detected
  if (nchar(aligned_forward) == 0 || nchar(aligned_reverse) == 0) return(list(consensus = NA_character_, consensus_length = 0, overlap = 0, identity = NA_real_, mismatches = NA_integer_, gaps = NA_integer_))

  ## Split the aligned sequences into individual positions
  forward_bases <- strsplit(aligned_forward, "")[[1]]
  reverse_bases <- strsplit(aligned_reverse, "")[[1]]

  ## Identify positions covered by nucleotides in both reads
  comparable_positions <- forward_bases != "-" & reverse_bases != "-"

  ## Retain identical nucleotides and replace mismatches or gaps with N
  consensus_bases <- ifelse(forward_bases == "-" | reverse_bases == "-", "N", ifelse(forward_bases == reverse_bases, forward_bases, "N"))

  ## Join all consensus positions and remove ambiguous positions from both ends
  consensus <- paste0(consensus_bases, collapse = "")
  consensus <- sub("^N+", "", consensus)
  consensus <- sub("N+$", "", consensus)

  ## Calculate alignment-quality statistics
  overlap <- sum(comparable_positions)
  mismatches <- sum(forward_bases != reverse_bases & comparable_positions)
  gaps <- sum(forward_bases == "-" | reverse_bases == "-")
  identity <- if (overlap > 0) sum(forward_bases[comparable_positions] == reverse_bases[comparable_positions]) / overlap else NA_real_

  ## Return the consensus sequence and quality-control statistics
  list(consensus = consensus, consensus_length = nchar(consensus), overlap = overlap, identity = identity, mismatches = mismatches, gaps = gaps)
}
