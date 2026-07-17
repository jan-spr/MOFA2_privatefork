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
        intercept_factors <- which(rowSums(r>0.75)>0)
        if (length(intercept_factors)) {
            warning(sprintf("Factor(s) %s are strongly correlated with the average values of features for at least one of your omics.\nSuch factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps. See 'plot_factor_intercept()' for a detailed breakdown.",paste(intercept_factors,collapse=", ")))
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
#' @name plot_factor_intercept
#' @description
#' Diagnostic plot for detecting \emph{intercept factors}, that
#' capture ovearally shift in the signal rather than
#' biological structure. The absolute Pearson correlation
#' between every factor and the per-sample mean feature values is computed
#' and displayed as a heatmap.
#' Factors above a correlation of 0.75 in at least one view are highlighted
#' (red outline, bold value), as these often reflect differences in the total
#' 'levels' between samples, sometimes caused by incomplete normalisation during
#' preprocessing.
#' @param object a trained \code{\linkS4class{MOFA}} object.
#' @details
#' Correlations are computed with \code{use = "pairwise.complete.obs"} to
#' tolerate missing values.
#' @return
#' A \code{\link[ggplot2]{ggplot}} object: a heatmap of factors (rows) by views
#' (columns). If the training data are not loaded
#' (\code{object@data_options$loaded} is \code{FALSE}) the function returns
#' \code{NULL} invisibly, as the plot cannot be computed.
#' @seealso \code{\link{run_mofa}}, \code{\link{load_model}}
#' @examples
#' \dontrun{
#' model <- load_model("model.hdf5", load_data = TRUE)
#' plot_factor_intercept(model)
#' }
#' @export
plot_factor_intercept <- function(object) {
  # sanity checks
  stopifnot(object@data_options[["loaded"]])
  if (is.null(object@data)) {stop("Data is Null")}
  
  # Compute for intercept factors
  factors <- do.call("rbind", get_factors(object))
  r <- suppressWarnings(t(do.call("rbind", lapply(object@data, function(x) {
    abs(cor(colMeans(do.call("cbind", x), na.rm = TRUE), factors, use = "pairwise.complete.obs"))
  }))))
  dimnames(r) <- list(factor = factors_names(object), view = views_names(object))
  # Heatmap
  df <- as.data.frame.table(r, responseName = "correlation")
  df$flag <- !is.na(df$correlation) & df$correlation > 0.75
  ggplot(df, aes(view, factor, fill = correlation)) +
    geom_tile(color = "grey90") +
    geom_tile(data = subset(df, flag), fill = NA, color = "red", linewidth = 1) +
    geom_text(aes(label = ifelse(is.na(correlation), "", sprintf("%.2f", correlation)),
                  fontface = ifelse(flag, "bold", "plain")), size = 3) +
    scale_fill_gradient(name = "|cor|", low = "white", high = "steelblue",
                        limits = c(0, 1), na.value = "grey95") +
    scale_y_discrete(limits = levels(df$factor)) +
    labs(x = NULL, y = NULL,
          title = "Factor correlation with mean feature values") +
    theme_minimal()
}