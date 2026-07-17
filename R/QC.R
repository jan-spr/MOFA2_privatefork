#' @importFrom stringi stri_enc_mark
.quality_control <- function(object, verbose = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Check views names
  if (verbose == TRUE) message("Checking views names...")
  stopifnot(!is.null(views_names(object)))
  stopifnot(!duplicated(views_names(object)))
  if (any(grepl("/", views_names(object)))) {
    stop("Some of the views names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    views_names(object) <- gsub(\"/\", \"-\", views_names(object))")
  }
  
  # Check groups names
  if (verbose == TRUE) message("Checking groups names...")
  if (any(grepl("/", groups_names(object)))) {
    stop("Some of the groups names contain `/` symbol, which is not supported.
    This can be fixed e.g. with:
    groups_names(object) <- gsub(\"/\", \"-\", groups_names(object))")
  }
  stopifnot(!is.null(groups_names(object)))
  stopifnot(!duplicated(groups_names(object)))
  
  # Check samples names
  if (verbose == TRUE) message("Checking samples names...")
  stopifnot(!is.null(samples_names(object)))
  stopifnot(!duplicated(unlist(samples_names(object))))
  enc <- stringi::stri_enc_mark(unlist(samples_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(samples_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following samples names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check features names
  if (verbose == TRUE) message("Checking features names...")
  stopifnot(!is.null(features_names(object)))
  stopifnot(!duplicated(unlist(features_names(object))))
  enc <- stringi::stri_enc_mark(unlist(features_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(features_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following features names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check dimensionalities in the input data
  if (verbose == TRUE) message("Checking dimensions...")
  N <- object@dimensions$N
  D <- object@dimensions$D
  for (i in views_names(object)) {
    for (j in groups_names(object)) {
      stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
      stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
      stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
      stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
    }
  }
  
  # Check that there are no features with complete missing values (across all groups)
  if (object@status == "untrained" || object@data_options[["loaded"]]) {
      if (verbose == TRUE) message("Checking there are no features with complete missing values...")
      for (i in views_names(object)) {
        if (!(is(object@data[[i]][[1]], "dgCMatrix") || is(object@data[[i]][[1]], "dgTMatrix"))) {
          tmp <- as.data.frame(sapply(object@data[[i]], function(x) rowMeans(is.na(x)), simplify = TRUE))
          if (any(unlist(apply(tmp, 1, function(x) mean(x==1)))==1))
            warning("You have features which do not contain a single observation in any group, consider removing them...")
        }
      }
    }
    
  # check dimensionalities of sample_covariates 
  if (verbose == TRUE) message("Checking sample covariates...")
  if(.hasSlot(object, "covariates") && !is.null(object@covariates)){
    stopifnot(ncol(object@covariates) == sum(object@dimensions$N))
    stopifnot(nrow(object@covariates) == object@dimensions$C)
    stopifnot(all(unlist(samples_names(object)) == colnames(object@covariates)))
  }
  
  # Sanity checks that are exclusive for an untrained model  
  if (object@status == "untrained") {
    
    # Check features names
    if (verbose == TRUE) message("Checking features names...")
    tmp <- lapply(object@data, function(x) unique(lapply(x,rownames)))
    for (x in tmp) stopifnot(length(x)==1)
    for (x in tmp) if (any(duplicated(x[[1]]))) stop("There are duplicated features names within the same view. Please rename")
    all_names <- unname(unlist(tmp))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (i in names(object@data)) {
      for (j in names(object@data[[i]])) {
        tmp <- which(rownames(object@data[[i]][[j]]) %in% duplicated_names)
        if (length(tmp)>0) {
          rownames(object@data[[i]][[j]])[tmp] <- paste(rownames(object@data[[i]][[j]])[tmp], i, sep="_")
        }
      }
    }
    
  # Sanity checks that are exclusive for a trained model  
  } else if (object@status == "trained") {
    # Check expectations
    if (verbose == TRUE) message("Checking expectations...")
    stopifnot(all(c("W", "Z") %in% names(object@expectations)))
    # if(.hasSlot(object, "covariates") && !is.null(object@covariates)) stopifnot("Sigma" %in% names(object@expectations))
    stopifnot(all(sapply(object@expectations$W, is.matrix)))
    stopifnot(all(sapply(object@expectations$Z, is.matrix)))
    
    # Check for intercept factors
    if (object@data_options[["loaded"]]) { 
      if (verbose == TRUE) message("Checking for intercept factors...")
      if (!is.null(object@data)) {
        factors <- do.call("rbind",get_factors(object))
        r <- suppressWarnings( t(do.call('rbind', lapply(object@data, function(x) 
          abs(cor(colMeans(do.call("cbind",x),na.rm=TRUE),factors, use="pairwise.complete.obs"))
        ))) )
        colnames(r) <- views_names(object)
        intercept_factors <- which(rowSums(r>0.75, na.rm=TRUE)>0)
        if (length(intercept_factors)) {
            # report factors with views they are flagged in
            flagged_views <- vapply(intercept_factors, function(i)
              paste(colnames(r)[which(r[i,]>0.75)], collapse=", "), character(1))
            factor_list <- paste(sprintf("%s (%s)", intercept_factors, flagged_views), collapse=", ")
            warning(sprintf("Factor(s) %s are strongly correlated with the per-sample averages of the indicated omics (See 'plot_factor_mean_cor()' for more detail).\nSuch factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.",factor_list))
        }
      }
    }
  
    # Check for correlated factors
    if (verbose == TRUE) message("Checking for highly correlated factors...")
    Z <- do.call("rbind",get_factors(object))
    noise <- matrix(rnorm(n=length(Z), mean=0, sd=1e-10), nrow(Z), ncol(Z))
    # suppress warnings locally instead of toggling globally
    tmp <- suppressWarnings(cor(Z+noise)); diag(tmp) <- NA
    if (max(tmp,na.rm=TRUE)>0.5) {
      warning("The model contains highly correlated factors (see `plot_factor_cor(MOFAobject)`). \nWe recommend that you train the model with less factors and that you let it train for a longer time.\n")
    }
  
  }
  
  return(object)  
}

#' @title Heatmap of factor correlation with mean feature values
#' @name plot_factor_mean_cor
#' @description
#' Plots (absolute) Pearson correlation between every factor and the per-sample mean feature
#' values, with correlations above 0.75 highlighted (red outline).
#' @param object a trained \code{\linkS4class{MOFA}} object.
#' @param split_by_groups logical. If \code{FALSE} (default), groups are pooled
#'   into a single heatmap; if \code{TRUE}, one panel per group.
#' @details
#' Highlighted factors often capture differences in the total 'levels' between
#' samples. Depending on the modality in question, this could be valid biological structure,
#' or a normalisation artifact of the data, e.g. from incomplete preprocessing.
#'
#' To decide whether these factors are of interest, it may help to take a closer look
#' at what signal the highlighted view-factor combinations capture -
#'  e.g. with \code{\link{plot_factor}}, \code{\link{plot_top_weights}}
#' or \code{\link{plot_factors_vs_cov}} (in the case with covariates).
#'
#' Correlations are computed with \code{use = "pairwise.complete.obs"} to
#' tolerate missing values. Pooled and per group scores are computed differently -
#' across all samples at once or split within each group.
#' @return
#' A \code{\link[ggplot2]{ggplot}} object: a heatmap of factors (rows) by views
#' (columns), faceted by group when \code{split_by_groups = TRUE}.
#' @seealso \code{\link{plot_factor}}, \code{\link{plot_weights}},
#'   \code{\link{plot_top_weights}}, \code{\link{plot_factors_vs_cov}},
#'   \code{\link{plot_variance_explained}}, \code{\link{load_model}}
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#'
#' # Correlation of each factor with the mean feature values per view
#' plot_factor_mean_cor(model)
#'
#' # Same, but computed within each group and shown as one panel per group
#' plot_factor_mean_cor(model, split_by_groups = TRUE)
#' @export
plot_factor_mean_cor <- function(object, split_by_groups = FALSE) {
  # sanity checks
  if (!object@data_options[["loaded"]]) {
    stop("No training data loaded in this object. This can be fixed by reloading with:
      load_model(file, load_data = TRUE)")
  }
  if (is.null(object@data)) {
    stop("The 'data' slot is empty, so correlations with the mean feature values
    cannot be computed.")
  }

  .mean_cor <- function(data_by_view, factors) {
    r <- suppressWarnings(
      do.call("rbind", lapply(data_by_view, function(x) {
        sample_means <- colMeans(x, na.rm = TRUE)
        abs(cor(sample_means, factors, use = "pairwise.complete.obs"))
      }
    )))
    return(t(r))
  }
  .label_dims <- function(r) {
    dimnames(r) <- list(factor = factors_names(object), view = views_names(object))
    return(r)
  }
  
  # Compute for intercept factors
  Z <- get_factors(object)
  if (split_by_groups) {
    # one correlation matrix per group, melted and stacked with a `group` column
    df <- do.call("rbind", lapply(groups_names(object), function(g) {
      data_g <- lapply(object@data, function(x) x[[g]])  # this group's matrix per view
      r <- .label_dims(.mean_cor(data_g, Z[[g]]))
      cbind(
        as.data.frame.table(r, responseName = "correlation"),
        group = g
      )
    }))
    df$group <- factor(df$group, levels = groups_names(object))
  } else {
    # all groups concatenated (samples pooled across groups)
    Z_comb <- do.call("rbind", Z)
    data_comb <- lapply(object@data, function(x) do.call("cbind", x))
    r <- .label_dims(.mean_cor(data_comb,Z_comb))
    df <- as.data.frame.table(r, responseName = "correlation")
  }
  
  # factors flagged in at least one view (in any group) get a bold axis label.
  df$flag <- !is.na(df$correlation) & df$correlation > 0.75
  factor_levels <- levels(df$factor)
  is_flagged <- factor_levels %in% unique(as.character(df$factor[df$flag]))
  axis_labels <- Map(
    function(nm, flagged) if (flagged) bquote(bold(.(nm))) else bquote(plain(.(nm))),
    factor_levels, is_flagged
  )

  # Heatmap
  p <- ggplot(df, aes(view, factor, fill = correlation)) +
    geom_tile(color = "grey90") +
    geom_tile(data = subset(df, flag), fill = NA, color = "red", linewidth = 1) +
    geom_text(aes(label = ifelse(is.na(correlation), "", sprintf("%.2f", correlation)),
                  fontface = ifelse(flag, "bold", "plain")), size = 3) +
    scale_fill_gradient(name = "|cor|", low = "white", high = "steelblue",
                        limits = c(0, 1), na.value = "grey95") +
    scale_y_discrete(limits = factor_levels, labels = axis_labels) +
    labs(x = NULL, y = NULL,
          title = "Factor correlation with per-sample mean") +
    theme_minimal()

  if (split_by_groups) p <- p + facet_wrap(~group, nrow = 1)
  p
}