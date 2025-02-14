get_pval <- function(chainA,
  chainB,
  overlap_freq,
  sample_size,
  n_repeats) {
# Probabilities
lenA <- length(chainA)
lenB <- length(chainB)
total_chain <- c(chainA, chainB)
total_len <- lenA + lenB

pA <- overlap_freq
pB <- 1 - pA
probs <- c(rep(pA / lenA, lenA), rep(pB / lenB, lenB))


# Draw sample indices for all replicates in one go
mix_idx <- matrix(
sample.int(
total_len,
size = sample_size * n_repeats,
prob = probs,
replace = TRUE
),
nrow = n_repeats
)

ref_idx <- matrix(sample.int(lenA, size = sample_size * n_repeats, replace = TRUE),
 nrow = n_repeats)

# # Use mapply to compute p-values for each replicate
# p_values <- mapply(
#   FUN = function(mix, ref) {
#     mixed_chain <- total_chain[mix]
#     reference_chain <- chainA[ref]
#     suppressWarnings(
#       epitree::compare_chains(reference_chain, mixed_chain)[, "Pr(>F)"][1]
#     )
#   },
#   mix = asplit(mix_idx, 1),
#   ref = asplit(ref_idx, 1)
# )

# Use vapply to compute p-values for each replicate
p_values <- vapply(seq_len(n_repeats), function(i) {
mixed_chain <- total_chain[mix_idx[i, ]]
reference_chain <- chainA[ref_idx[i, ]]
suppressWarnings(
epitree::compare_chains(reference_chain, mixed_chain)[, "Pr(>F)"][1]
)
}, numeric(1))

return(p_values)
}