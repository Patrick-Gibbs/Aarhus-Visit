
library(susieR)
args <- commandArgs(trailingOnly = TRUE)
input_sum <- args[1]
input_ld_matrix <- args[2]
output_probs <- args[3]
n <- as.integer(args[4])

# Optional parameters
L <- if(length(args) >= 5) as.integer(args[5]) else NULL
max_iter <- if(length(args) >= 6) as.integer(args[6]) else NULL
refine <- if(length(args) >= 7) as.logical(args[7]) else NULL

gwas_sum_stats <- read.csv(input_sum, sep=' ')
# read space sperated matrix for R
R = as.matrix(read.table(input_ld_matrix, sep=' '))
# need to put in variance covariance matrix.

# Create a list of arguments for susie_rss
args_list <- list(z=gwas_sum_stats$Z, R=R, estimate_residual_variance = TRUE, n = n)

# Add optional parameters to the list if they are not NULL
if (!is.null(L)) args_list$L <- L
if (!is.null(max_iter)) args_list$max_iter <- max_iter
if (!is.null(refine)) args_list$refine <- refine

# Call susie_rss with the arguments
result <- do.call(susie_rss, args_list)

# save the results$pip as a csv such that the header is Probability
output_probs <- paste(output_probs, ".probs", sep="")
write.csv(data.frame(Probability=result$pip), output_probs, row.names=FALSE, quote=FALSE)
