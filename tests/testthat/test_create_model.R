library(MOFA2)

test_that("a model can be created from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	expect_warning(create_mofa(list("view1" = m)))  # no feature names provided
	rownames(m) <- paste("feature", seq_len(nrow(m)), paste = "", sep = "")
	expect_s4_class(create_mofa(list("view1" = m)), "MOFA")
	expect_error(create_mofa(m))
})

test_that("a model can be created from a list of sparse matrices", {
	skip_if_not_installed("Matrix")
	library(Matrix)
	# Generate a sparse matrix
	m <- matrix(rnorm(100 * 5), ncol = 5) %*% t(matrix(rnorm(5 * 50), ncol = 5))
	m[sample(1:nrow(m), 100, replace = TRUE), sample(1:ncol(m), 100, replace = TRUE)] <- 0
	m <- Matrix(m, sparse = TRUE)
	# Set feature names
	rownames(m) <- paste("feature_", seq_len(nrow(m)), paste = "", sep = "")
	# Set sample names
	colnames(m) <- paste("sample_", seq_len(ncol(m)), paste = "", sep = "")
	# Test if a sparse matrix can be imported to the MOFA
	expect_s4_class(create_mofa(list("view1" = m)), "MOFA")
})

test_that("a model can be created from a Seurat object", {
	skip_if_not_installed("Seurat")
	skip_if_not_installed("SeuratObject")
	library(Seurat)
	library(Matrix)
	m <- as(readMM('matrix.mtx'),'dgCMatrix')
	genes <- read.delim('genes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,2]
	cells <- read.delim('barcodes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
	colnames(m) <- cells
	rownames(m) <- genes
	srt <- SeuratObject::CreateSeuratObject(m)
	# only for testing purpose, should use scale.data
	expect_s4_class(create_mofa(srt, features = genes, layer = "counts"), "MOFA")
})

test_that("a list of matrices per view is split correctly into a nested list of matrices according to samples groups", {
	n_groups <- 3
	# Create view 1
	m <- as.matrix(read.csv('matrix.csv'))
	rownames(m) <- paste("feature", seq_len(nrow(m)), paste = "", sep = "")
	colnames(m) <- paste("sample", seq_len(ncol(m)), paste = "", sep = "")
	# Add second view
	m2 <- m[1:(nrow(m)/3),]
	rownames(m2) <- paste("view2", rownames(m2), sep = "_")
	# Define multiple groups
	samples_groups <- sample(x = paste0("group", 1:n_groups), replace = TRUE, size = ncol(m))
	# Split the data
	data_split <- .split_data_into_groups(list("view1" = m, "view2" = m2), samples_groups)
	# Check group assignments
	for (g in 1:n_groups) {
		g_name <- paste0("group", g)
		expect_equal(colnames(data_split[[1]][[g_name]]), colnames(m)[which(samples_groups == g_name)])
		expect_equal(colnames(data_split[[2]][[g_name]]), colnames(m)[which(samples_groups == g_name)])
	}
})


test_that("a model can be created from a MultiAssayExperiment Object", {
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)

	# Import and preprocess the miniACC fixture (see helper-mae.R)
	mae_sub <- make_test_mae()

	# create MOFA model
	# "" drops the Mutations experiment, leaving 3 views
	model <- create_mofa(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		extract_metadata = TRUE
	)

	# do checks
	# class check
	expect_s4_class(model, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(model)$M,3) # right number of views

	# right views retained (Mutations dropped via "")
	expect_equal(
		views_names(model),
		c("RNASeq2GeneNorm", "RPPAArray", "miRNASeqGene")
	)


    # right data matrix (colnames differ: the sampleMap maps assay barcodes
    # to primary sample names)
    expect_equal(
        get_data(model, views = c("RPPAArray"))$RPPAArray$group1,
        assay(mae_sub[["RPPAArray"]]),
        ignore_attr = TRUE
    )

	# extract_metadata should propagate colData onto the MOFA object
	meta <- samples_metadata(model)
	cd <- as.data.frame(colData(mae_sub))
	metadata_col <- colnames(cd)[1]
	expect_true(all(c("sample", "group", metadata_col) %in% colnames(meta)))
	expect_equal(meta[[metadata_col]], cd[meta$sample, metadata_col])

	#check also with HintikkaXOData

})

test_that("a model can be created from a MultiAssayExperiment Object - HintikkaXOData data", {
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)
	library(mia)

	#check also with HintikkaXOData
	data("HintikkaXOData")

	# Prepare data 
	mae <- HintikkaXOData
	altExp(mae[[1]], "asd") <- mae[[1]]

	MOFAobject <- create_mofa(
		mae,
		alt_experiments = list("asd", "main", "main"),
		assays = list("counts",NULL,"signals")
	)

	# do checks
	# class check
	expect_s4_class(MOFAobject, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(MOFAobject)$M,2) # right number of views

	MOFAobject <- create_mofa(
		mae,
		experiments = c("microbiota","biomarkers"),
		alt_experiments = list("asd", "main"),
		assays = list("counts","signals")
	)

	# do checks
	# class check
	expect_s4_class(MOFAobject, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(MOFAobject)$M,2) # right number of views
})

test_that("a model can be created from a SingleCellExperiment Object", {
	skip_if_not_installed("SingleCellExperiment")
	skip_if_not_installed("MOFAdata")
	skip_if_not_installed("data.table")
	library(SingleCellExperiment)
	library(data.table)
	
	# Import CLL data
	utils::data("CLL_data", package = "MOFAdata")
	CLL_metadata <- data.table::fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
	IGHV <- CLL_metadata$IGHV 
	IGHV_filled <- ifelse(is.na(IGHV), 'NA', IGHV)
	CLL_metadata$IGHV_filled <- IGHV_filled

	# Create SCE Object
	sce <- SingleCellExperiment(assays = list(mRNA = CLL_data$mRNA), colData = CLL_metadata)
	altExps(sce) <- list(
		Drugs = SummarizedExperiment(list(expr = CLL_data$Drugs)),
		Methylation = SummarizedExperiment(list(expr = CLL_data$Methylation)),
		Mutations = SummarizedExperiment(list(expr = CLL_data$Mutations))
	)

	# create MOFA model
	MOFAobject <- create_mofa_from_SingleCellExperiment(
		sce,
		alt_experiments = c("Main","Drugs", "Methylation", "Mutations"),
		assays = c("mRNA","expr","expr","expr"),
		groups = "IGHV_filled",
		extract_metadata = TRUE
	)

	# do checks
	# class check
	expect_s4_class(MOFAobject, "MOFA")

	# extract_metadata should propagate the SCE colData onto the MOFA object.
	# Samples are reordered by group, so match on the sample column.
	meta <- samples_metadata(MOFAobject)
	cd <- as.data.frame(colData(sce))
	expect_true(all(c("age", "IGHV_filled") %in% colnames(meta)))
	expect_equal(meta$IGHV_filled, cd[meta$sample, "IGHV_filled"])
	expect_equal(meta$age, cd[meta$sample, "age"])

	
    # right data matrix - with groups
    expect_equal(
        get_data(MOFAobject, views = c("Main"), groups = c("0"))$Main$`0`,
        CLL_data$mRNA[, CLL_metadata$IGHV == 0 & !is.na(CLL_metadata$IGHV)]
    )
    expect_equal(
        get_data(MOFAobject, views = c("Main"), groups = c("1"))$Main$`1`,
        CLL_data$mRNA[, CLL_metadata$IGHV == 1 & !is.na(CLL_metadata$IGHV)]
    )

    # check dimensions
    expect_equal(get_dimensions(MOFAobject)$G, 3) # right number of groups

    # right data matrix - without groups
    # create MOFA model
    MOFAobject <- create_mofa_from_SingleCellExperiment(
        sce,
        alt_experiments = c("Main", "Drugs", "Methylation", "Mutations"),
        assays = c("mRNA", "expr", "expr", "expr")
    )
    expect_equal(
        get_data(MOFAobject, views = c("Main"))$Main$group1,
        CLL_data$mRNA
    )
})

test_that("MAE experiments can be selected by numeric index", {
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)

	mae_sub <- make_test_mae()  # 4 experiments (see helper-mae.R)

	# Select experiments 1, 2 and 4 by index (drops Mutations at index 3).
	# assays are matched positionally to the *selected* experiments.
	model <- create_mofa(
		mae_sub,
		experiments = c(1, 2, 4),
		assays = list("log1p", "exprs", "log1p")
	)

	expect_s4_class(model, "MOFA")
	expect_equal(get_dimensions(model)$M, 3)
	expect_equal(
		views_names(model),
		c("RNASeq2GeneNorm", "RPPAArray", "miRNASeqGene")
	)


    # selecting by index and by name should give the same object
    model_by_name <- create_mofa(
        mae_sub,
        experiments = c("RNASeq2GeneNorm", "RPPAArray", "miRNASeqGene"),
        assays = list("log1p", "exprs", "log1p")
    )
    expect_equal(get_data(model), get_data(model_by_name))
})

test_that("SCE alt_experiments can be selected by numeric index", {
	skip_if_not_installed("SingleCellExperiment")
	library(SingleCellExperiment)
	library(SummarizedExperiment)

	set.seed(1)
	main <- matrix(rnorm(20 * 10), nrow = 20,
		dimnames = list(paste0("g1_", 1:20), paste0("s_", 1:10)))
	alt <- matrix(rnorm(15 * 10), nrow = 15,
		dimnames = list(paste0("g2_", 1:15), paste0("s_", 1:10)))
	sce <- SingleCellExperiment(
		assays = list(counts = main),
		altExps = list(alt_data = SummarizedExperiment(list(counts = alt)))
	)

	# reference the altExp by numeric index (1 == "alt_data"), mixed with the
	# main experiment. Numeric entries require the list() form of alt_experiments.
	model <- create_mofa_from_SingleCellExperiment(
		sce,
		alt_experiments = list("main", 1),
		assays = c("counts", "counts")
	)

	expect_s4_class(model, "MOFA")
	expect_equal(get_dimensions(model)$M, 2)
	# numeric index should resolve to the altExp name
	expect_equal(views_names(model), c("Main", "alt_data"))

    # index and name should select the same altExp data
    model_by_name <- create_mofa_from_SingleCellExperiment(
        sce,
        alt_experiments = c("main", "alt_data"),
        assays = c("counts", "counts")
    )
    expect_equal(get_data(model), get_data(model_by_name))
})

test_that("create_mofa() dispatches to the SingleCellExperiment constructor", {
    skip_if_not_installed("SingleCellExperiment")
    library(SingleCellExperiment)
    library(SummarizedExperiment)

    set.seed(1)
    logc <- matrix(rnorm(20 * 10),
        nrow = 20,
        dimnames = list(paste0("g1_", 1:20), paste0("s_", 1:10))
    )
    alt <- matrix(rnorm(15 * 10),
        nrow = 15,
        dimnames = list(paste0("g2_", 1:15), paste0("s_", 1:10))
    )
    sce <- SingleCellExperiment(
        assays = list(logcounts = logc),
        altExps = list(alt_data = SummarizedExperiment(list(counts = alt)))
    )

    # Regression - create_mofa(sce) with no assays argument (with one view names view after assay)
    model <- create_mofa(sce)
    expect_s4_class(model, "MOFA")
    expect_equal(get_dimensions(model)$M, 1)
    expect_equal(views_names(model), "logcounts")
    expect_equal(get_data(model)$logcounts$group1, logc)

    # The dispatcher forwards alt_experiments / assays (via ...) to the constructor
    model2 <- create_mofa(
        sce,
        alt_experiments = c("Main", "alt_data"),
        assays = c("logcounts", "counts")
    )
    expect_s4_class(model2, "MOFA")
    expect_equal(get_dimensions(model2)$M, 2)
    expect_equal(views_names(model2), c("Main", "alt_data"))

    # assays length must match alt_experiments length
    expect_error(
        create_mofa(sce, alt_experiments = c("main", "alt_data"), assays = "logcounts"),
        "must match"
    )
})