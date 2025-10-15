#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    if (!requireNamespace("getopt", quietly = TRUE)) {
        stop("Package 'getopt' is required. Install with install.packages('getopt').")
    }
})

# --- Utility: escape regex metacharacters in a literal label ---
escape_regex <- function(x) gsub("([][(){}.+*?^$|\\\\])", "\\\\\\1", x)

# --- Mid-root by splitting the outgroup edge (string-based, minimal change) ---
# Assumptions:
# - Input contains exactly one Newick tree.
# - Outgroup appears at the top clade level as "(<OUTGROUP>:<len>)".
# - We remove that subtree, split its branch length in half, and rebuild a new root.
mid_root_iqtree <- function(tree_string, outgroup_label) {
    # Normalize: collapse to one line and strip all whitespace
    tree_string <- gsub("\\s+", "", paste(tree_string, collapse = ""))
    
    # Pattern for "(OUTGROUP:<number>)"
    # <number> supports integer/decimal/scientific notation
    num_pat <- "([0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
    pat <- paste0("\\(", escape_regex(outgroup_label), ":", num_pat, "\\)")
    
    # Locate the outgroup subtree
    all_m <- gregexpr(pat, tree_string, perl = TRUE)
    mm <- regmatches(tree_string, all_m)[[1]]
    if (length(mm) == 0) {
        stop(sprintf("Outgroup '%s' not found as a subtree '(LABEL:length)'.", outgroup_label))
    }
    if (length(mm) > 1) {
        stop(sprintf(
            "Multiple outgroup matches found (%d). Ensure the label is unique. Examples: %s",
            length(mm), paste(head(mm, 3), collapse = " | ")
        ))
    }
    extracted_string <- mm[1]
    
    # Extract numeric length via capture group
    m1 <- regexec(pat, tree_string, perl = TRUE)
    cap <- regmatches(tree_string, m1)[[1]]
    out_len <- as.numeric(cap[2])
    if (!is.finite(out_len)) stop("Failed to parse the outgroup branch length.")
    half_len <- out_len / 2
    
    # Remove the outgroup safely: try ",<match>" then "<match>,"
    modified <- sub(paste0(",", extracted_string), "", tree_string, fixed = TRUE)
    if (identical(modified, tree_string)) {
        modified <- sub(paste0(extracted_string, ","), "", tree_string, fixed = TRUE)
    }
    if (identical(modified, tree_string)) {
        stop("Could not remove the outgroup subtree; it may not be a direct sibling at this level.")
    }
    
    # Strip trailing semicolon, rebuild balanced root, and return
    modified <- sub(";$", "", modified)
    final <- sprintf("(%s:%.10g,%s:%.10g);", modified, half_len, outgroup_label, half_len)
    return(final)
}

# --- CLI (getopt) ---
usage <- function() {
    cat("
Mid-root a IQ-TREE treefile by splitting the outgroup edge (string-based).

Usage:
  Rscript midroot_iqtree.R --tree <file> --outgroup <label> --out <file>

Options:
  -t, --tree       Input Newick tree file (single tree)
  -g, --outgroup   Outgroup tip label (exact match)
  -o, --out        Output tree file
  -h, --help       Show this help and exit
")
}

spec <- matrix(c(
    "help",     "h", 0, "logical",
    "tree",     "t", 1, "character",
    "outgroup", "g", 1, "character",
    "out",      "o", 1, "character"
), byrow = TRUE, ncol = 4)

opt <- getopt::getopt(spec)
if (!is.null(opt$help)) { usage(); quit(status = 0) }
if (is.null(opt$tree) || is.null(opt$outgroup) || is.null(opt$out)) {
    usage(); quit(status = 1)
}

# --- Read, process, write ---
txt <- readLines(opt$tree, warn = FALSE)
if (length(txt) == 0) stop("Failed to read the tree file.")
res <- mid_root_iqtree(txt, opt$outgroup)
cat(res, file = opt$out)
