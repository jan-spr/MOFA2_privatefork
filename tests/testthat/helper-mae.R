# Shared fixtures for the create_mofa / mofa2 tests.
# testthat sources helper-*.R automatically before running the tests.

# Build the miniACC-derived MultiAssayExperiment used by several tests.
# Selects four experiments, makes feature names unique per experiment (so the
# duplicated-feature renaming does not kick in) and adds a "log1p" assay to the
# two RNA-based experiments.
# Assumes MultiAssayExperiment and SummarizedExperiment are attached, which the
# calling tests guarantee via skip_if_not_installed() + library().
make_test_mae <- function() {
	data(miniACC, package = "MultiAssayExperiment", envir = environment())
	miniACC <- MultiAssayExperiment::intersectColumns(miniACC)
	mae_sub <- miniACC
	MultiAssayExperiment::experiments(mae_sub) <- MultiAssayExperiment::experiments(miniACC)[
		c("RNASeq2GeneNorm", "RPPAArray", "Mutations", "miRNASeqGene")
	]

	# Make feature names unique per experiment
	rownames(mae_sub[["RNASeq2GeneNorm"]]) <- paste0(rownames(mae_sub[["RNASeq2GeneNorm"]]), "_1")
	rownames(mae_sub[["RPPAArray"]])       <- paste0(rownames(mae_sub[["RPPAArray"]]), "_2")
	rownames(mae_sub[["Mutations"]])       <- paste0(rownames(mae_sub[["Mutations"]]), "_3")
	rownames(mae_sub[["miRNASeqGene"]])    <- paste0(rownames(mae_sub[["miRNASeqGene"]]), "_4")

	# Add a log1p assay to the count-based experiments
	for (nm in c("RNASeq2GeneNorm", "miRNASeqGene")) {
		se <- mae_sub[[nm]]
		SummarizedExperiment::assay(se, "log1p") <- log1p(SummarizedExperiment::assay(se, "exprs"))
		mae_sub[[nm]] <- se
	}

	mae_sub
}
