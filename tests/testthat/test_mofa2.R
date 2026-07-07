test_that("mofa2 wrapper correctly initializes MOFAobject", {
	#check that manually setting the parameters creates the same model as the mofa2() wrapper

	# 1. get sample data
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)

	# Import and preprocess the miniACC fixture (see helper-mae.R)
	mae_sub <- make_test_mae()

	# 2. create MOFA instances
	MOFA_init <- create_mofa(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		extract_metadata = TRUE
	)

	data_opts <- get_default_data_options(MOFA_init)
	model_opts <- get_default_model_options(MOFA_init)
	model_opts$num_factors <- 4
	train_opts <- get_default_training_options(MOFA_init)
	train_opts$seed <- 42
	train_opts$convergence_mode <- "fast"

	MOFA_prep1 <- prepare_mofa(MOFA_init,
		data_options = data_opts,
		model_options = model_opts,
		training_options = train_opts
	)

	MOFA_prep2 <- mofa2(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		num_factors = 4,
		seed = 42,
		convergence_mode = "fast"
	)
	
	expect_identical(MOFA_prep1,MOFA_prep2)

	# check manual parameter values
	MOFA_prep2 <- mofa2(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		num_factors = 4,
		seed = 1337,
		convergence_mode = "fast"
	)
	expect_equal(MOFA_prep2@training_options$seed,1337)
	expect_equal(MOFA_prep2@training_options$convergence_mode,"fast")
	expect_equal(MOFA_prep2@training_options$maxiter,1000)
	
	# add: SCE with altexperiments?
})

test_that("mofa2 wrapper rejects misspelled arguments", {
    # a typo'd option must error rather than being silently dropped
    skip_if_not_installed("MultiAssayExperiment")
    library(MultiAssayExperiment)
    library(SummarizedExperiment)

    mae_sub <- make_test_mae()

    expect_error(
        mofa2(
            mae_sub,
            assays = c("log1p", "exprs", "", "log1p"),
            num_factrs = 4
        ),
        "Unrecognised argument"
    )
})

test_that("mofa2 dispatches on a list of matrices and matches the manual call", {
	# matrix input needs no extra packages
	set.seed(1)
	mtx <- make_example_data(n_views = 2, n_features = 30, n_samples = 40)$data

	MOFA_wrap <- mofa2(mtx, num_factors = 4, seed = 7)
	expect_s4_class(MOFA_wrap, "MOFA")

	# same object as an explicit create_mofa + prepare_mofa (subsumes the
	# individual status / views / option-slot checks)
	MOFA_init  <- create_mofa(mtx)
	model_opts <- get_default_model_options(MOFA_init)
	model_opts$num_factors <- 4
	train_opts <- get_default_training_options(MOFA_init)
	train_opts$seed <- 7
	MOFA_manual <- prepare_mofa(MOFA_init, model_options = model_opts, training_options = train_opts)

	expect_identical(MOFA_wrap, MOFA_manual)
})

test_that("mofa2 fills stochastic options only when stochastic training is enabled", {
	set.seed(1)
	mtx <- make_example_data(n_views = 2, n_features = 30, n_samples = 40)$data

	MOFA_default <- mofa2(mtx, num_factors = 4)
	expect_length(MOFA_default@stochastic_options, 0)

	# enabling stochastic on small data always warns (N << 10,000) - expected
	MOFA_stoch <- suppressWarnings(
		mofa2(mtx, stochastic = TRUE, learning_rate = 0.75, num_factors = 4)
	)
	expect_gt(length(MOFA_stoch@stochastic_options), 0)
	expect_equal(MOFA_stoch@stochastic_options$learning_rate, 0.75)
})

test_that("mofa2 fills mefisto options only when covariates are provided", {
	set.seed(1)
	mtx <- make_example_data(n_views = 2, n_features = 30, n_samples = 40)$data

	# without covariates the mefisto slot stays empty
	expect_length(mofa2(mtx, num_factors = 4)@mefisto_options, 0)

	# a single continuous covariate, named to match the samples
	covs <- matrix(seq_len(40), nrow = 1,
	               dimnames = list("time", colnames(mtx[[1]])))

	MOFA_mefisto <- mofa2(mtx, covariates = covs, scale_cov = TRUE, num_factors = 4)
	expect_gt(length(MOFA_mefisto@covariates), 0)          # set_covariates ran
	expect_gt(length(MOFA_mefisto@mefisto_options), 0)     # options were stored
	expect_true(MOFA_mefisto@mefisto_options$scale_cov)    # ... routed to the mefisto group
})