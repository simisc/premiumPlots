#' Summarise the posterior distribution of the variable selection
#'     parameter \code{rho}. Modified version of
#'     \code{\link[PReMiuM]{summariseVarSelectRho}}, allowing plotting
#'     functions to access \code{rho} summaries.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @return A \code{\link[tibble]{tibble}} containing the mean, median,
#'     lower (0.05) and upper (0.95) quantiles for each covariate,
#'     and a ranking from the covariate with the highest value of
#'     \code{rho} to the lowest.
summariseVarSelectRhoFromRiskProfObj <- function (riskProfObj) {

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    rhoFileName <-
        file.path(directoryPath, paste(fileStem, "_rho.txt", sep = ""))
    rhoMat <- matrix(scan(rhoFileName, what = double(), quiet = T),
                     ncol = nCovariates,
                     byrow = T)
    firstLine <- ifelse(reportBurnIn, nBurn / nFilter + 2, 1)
    lastLine <-
        (nSweeps + ifelse(reportBurnIn, nBurn + 1, 0)) / nFilter
    rhoMat <- rhoMat[firstLine:lastLine, ]
    rhoMean <- apply(rhoMat, 2, mean)
    rhoMedian <- apply(rhoMat, 2, median)
    rhoLowerCI <- apply(rhoMat, 2, quantile, 0.05)
    rhoUpperCI <- apply(rhoMat, 2, quantile, 0.95)

    tibble::tibble(
        var = covNames,
        rhoMean = rhoMean,
        rhoMedian = rhoMedian,
        rhoLowerCI = rhoLowerCI,
        rhoUpperCI = rhoUpperCI
    ) %>%
        dplyr::mutate(rhoRank = rank(-rhoMean))
}

#' Plots the (discrete) covariate profiles, facetted by cluster.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed. If
#'     a single integer, this is intepreted as the 'top N' covariates,
#'     as ranked by their value of \code{rho} (requires the model to
#'     have been fitted with variable selection).
#' @param rhoMinimum As an alternative to \code{whichCovariates}, plot
#'     plot all covariates whose value of \code{rho} is greater than
#'     some threshold (requires the model to have been fitted with
#'     variable selection).
#' @param useProfileStar o be set equal to TRUE only if a variable
#'     selection procedure has been run. The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#' @param covariate_levels Vector of integer values taken by the
#'     (discrete) covariates, specifying the order of the facets.
#' @param covariate_labels Vector of strings giving labels for each
#'     level of the covariate, in the same order as \code{covariate_levels}.
plotProfilesByCluster <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoMinimum = NULL,
                                   useProfileStar = TRUE,
                                   covariate_levels = NULL,
                                   covariate_labels = NULL) {
    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoMinimum = rhoMinimum,
        useProfileStar = useProfileStar
    )

    if (is.null(rhoMinimum)) {
        rhoMinimum <- min(profileDF$rhoMean)
    }

    if (!is.null(covariate_levels) & !is.null(covariate_labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = covariate_levels,
                                            labels = covariate_labels))
    }

    covtab <- profileDF %>%
        dplyr::group_by(cluster, category, covname, fillColor, rhoMean, rhoRank) %>%
        dplyr::summarise(prop = mean(est)) %>%
        dplyr::mutate(covrho = sprintf("%s (%i%%)", covname, round(100 * rhoMean)))

    expected_proportions <-
        tapply(profileDF$mean, profileDF$category, mean)
    expected_proportions <-
        cumsum(expected_proportions)[-length(expected_proportions)]

    ggplot2::ggplot(covtab,
                    ggplot2::aes(
                        x = reorder(covrho, rhoRank),
                        y = prop,
                        fill = factor(category),
                        alpha = fillColor != "avg"
                    )) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::facet_grid(cluster ~ .) +
        ggplot2::geom_hline(yintercept = expected_proportions, linetype = "dashed") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )) +
        ggplot2::labs(
            x = sprintf("Covariates (rho >= %.i%%)", round(100 * rhoMinimum)),
            y = "Proportion (by cluster)",
            title = "Covariate profiles"
        ) +
        ggplot2::scale_fill_discrete(name = "Marker\nprevalence\ncategory") +
        ggplot2::scale_alpha_discrete(name = "Proportion\ndifferent\nfrom mean",
                                      range = c(0.25, 1))
}

#' Plots the (discrete) covariate profiles, facetted by covariate,
#'     similar to the layout in \code{\link[PReMiuM]{plotRiskProfile}}.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed. If
#'     a single integer, this is intepreted as the 'top N' covariates,
#'     as ranked by their value of \code{rho} (requires the model to
#'     have been fitted with variable selection).
#' @param rhoMinimum As an alternative to \code{whichCovariates}, plot
#'     plot all covariates whose value of \code{rho} is greater than
#'     some threshold (requires the model to have been fitted with
#'     variable selection).
#' @param useProfileStar o be set equal to TRUE only if a variable
#'     selection procedure has been run. The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#' @param covariate_levels Vector of integer values taken by the
#'     (discrete) covariates, specifying the order of the facets.
#' @param covariate_labels Vector of strings giving labels for each
#'     level of the covariate, in the same order as \code{covariate_levels}.
plotCovariateProfiles <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoMinimum = NULL,
                                   useProfileStar = TRUE,
                                   covariate_levels = NULL,
                                   covariate_labels = NULL) {
    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoMinimum = rhoMinimum,
        useProfileStar = useProfileStar
    )

    length(unique(profileDF$covname)) <= 10 ||
        stop("Don't plot so many covariates at once.")

    if (!is.null(covariate_levels) & !is.null(covariate_labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = covariate_levels,
                                            labels = covariate_labels))
    }

    profileDF <- profileDF %>%
        dplyr::mutate(covrho = sprintf("%s (%i%%)", covname, round(100 * rhoMean)))

    cols <- c(high = "#CC0033",
              low = "#0066CC",
              avg = "#33CC66")

    ggplot2::ggplot(profileDF, ggplot2::aes(x = factor(cluster), col = factor(fillColor))) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean), linetype = "dashed") +
        ggplot2::geom_boxplot(ggplot2::aes(y = est, fill = factor(fillColor)),
                              col = "black",
                              outlier.size = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = lower), size = 1.5, shape = 4) +
        ggplot2::geom_point(ggplot2::aes(y = upper), size = 1.5, shape = 4) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::theme(legend.position = "none",
                       strip.text = ggplot2::element_text(size = 6)) +
        ggplot2::labs(x = "Cluster", title = "Covariate profiles", y = "Probability") +
        ggplot2::facet_grid(factor(category) ~ reorder(covrho, rhoRank))
}

#' @export
tabulateCovariateProfiles <- function (riskProfObj,
                                       whichCovariates = NULL,
                                       rhoMinimum = NULL,
                                       useProfileStar = TRUE) {
    # This now does most of the work for plotCovariateProfiles
    # Use same output for new covariate profiles plot
    # Vectorisation of inner loop made this 3x faster...

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    xModel == "Discrete" ||
        stop(
            "Only implemented for discrete covariates.\n
            Use plotRiskProfile() for Normal- or Mixed-covariate models."
        )

    if (useProfileStar) {
        profile <- profileStar
    }

    rhotab <- summariseVarSelectRhoFromRiskProfObj(riskProfObj)

    if (is.null(whichCovariates) & !is.null(rhoMinimum)) {
        # If specify minimum rho value instead of which covariates
        rhotab <- filter(rhotab, rhoMean >= rhoMinimum)
        whichCovariates <- rhotab$var
    }

    if (length(whichCovariates) == 1) {
        # If whichCovariates is a single integer, treat it as the 'top N'
        rhotab <- arrange(rhotab, rhoRank)
        whichCovariates <- rhotab$var[1:whichCovariates]
    }

    if (!is.null(whichCovariates)) {
        if (!is.numeric(whichCovariates)) {
            whichCovariates <- match(whichCovariates, covNames)
        }
        covNames <- covNames[whichCovariates]
        nCovariates <- length(whichCovariates)
        profile <- profile[, , whichCovariates, ]
        nCategories <- nCategories[whichCovariates]
    }

    # clusters numbered in consistent order (same in plotResponse and plotClusterSizes)
    orderStat <- apply(risk, 2, median)
    meanSortIndex <- order(orderStat, decreasing = F)
    clusterSizes <- clusterSizes[meanSortIndex]
    profile <- profile[, meanSortIndex, , ]

    profDFlist = list()

    for (j in 1:nCovariates) {
        miniDFlist = list() # still growing lists, not ideal...

        for (k in 1:nCategories[j]) {
            probMat <- profile[, , j, k]
            probMeans <- apply(probMat, 2, mean)

            probMean <-
                sum(probMeans * clusterSizes) / sum(clusterSizes)
            probLower <- apply(probMat, 2, quantile, 0.05)
            probUpper <- apply(probMat, 2, quantile, 0.95)

            clusterDF <- tibble::tibble(
                cluster = 1:nClusters,
                category = k - 1,
                mean = probMean,
                lower = probLower,
                upper = probUpper
            ) %>%
                dplyr::mutate(
                    fillColor = ifelse(
                        lower > mean,
                        "high",
                        ifelse(upper < mean, "low", "avg")
                    ),
                    fillColor = as.character(fillColor)
                )

            profileDF <-
                tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                               est = c(probMat)) %>%
                dplyr::left_join(clusterDF, by = "cluster")

            miniDFlist[[k]] <- profileDF
        }

        profileDF <- dplyr::bind_rows(miniDFlist)
        profDFlist[[covNames[j]]] <- profileDF
    }

    dplyr::bind_rows(profDFlist, .id = "covname") %>%
        dplyr::left_join(rhotab, by = c("covname" = "var"))
}

#' Plots the (categorical) response profiles, similar to the layout
#'     of the 'risk' plot in \code{\link[PReMiuM]{plotRiskProfile}}.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param response_levels Vector of integer values taken by the
#'     (categorical) response, specifying the order of the facets.
#' @param response_labels Vector of strings giving labels for each
#'     level of the response, in the same order as \code{response_levels}.
plotResponse <- function (riskProfObj,
                          response_levels = NULL,
                          response_labels = NULL) {
    # Restructured like tabulateCovariateProfile()
    # 33x times faster than before
    # Corrected ordering of clusters along x-axis!

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    if (!includeResponse) {
        stop("if !includeResponse, why are you trying to plotResponse()?")
    }

    if (yModel != "Categorical") {
        stop("Only for categorical response, use plotRiskProfile() for others...")
    }

    orderStat <- apply(risk, 2, median)
    meanSortIndex <- order(orderStat, decreasing = F)
    clusterSizes <- clusterSizes[meanSortIndex]
    risk <- risk[, meanSortIndex, ]

    miniDFlist = list() # still growing lists, not ideal...

    for (k in 1:nCategoriesY) {
        probMat <- risk[, , k] # 2D matrix for categ. k
        probMeans <- apply(probMat, 2, mean, trim = 0.005)

        probMean <-
            sum(probMeans * clusterSizes) / sum(clusterSizes)
        probLower <- apply(probMat, 2, quantile, 0.05)
        probUpper <- apply(probMat, 2, quantile, 0.95)

        clusterDF <- tibble::tibble(
            cluster = 1:nClusters,
            category = k - 1,
            mean = probMean,
            lower = probLower,
            upper = probUpper
        ) %>%
            dplyr::mutate(
                fillColor = ifelse(lower > mean, "high",
                                   ifelse(upper < mean, "low", "avg")),
                fillColor = as.character(fillColor)
            )

        profileDF <-
            tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                           est = c(probMat)) %>%
            dplyr::left_join(clusterDF, by = "cluster")

        miniDFlist[[k]] <- profileDF
    }

    profileDF <- dplyr::bind_rows(miniDFlist)

    if (!is.null(response_levels) & !is.null(response_labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = response_levels,
                                            labels = response_labels))
    }

    cols <- c(high = "#CC0033",
              low = "#0066CC",
              avg = "#33CC66")

    ggplot2::ggplot(profileDF, ggplot2::aes(x = factor(cluster), col = factor(fillColor))) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean), linetype = "dashed") +
        ggplot2::geom_boxplot(ggplot2::aes(y = est, fill = factor(fillColor)),
                              col = "black",
                              outlier.size = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = lower), size = 1.5, shape = 4) +
        ggplot2::geom_point(ggplot2::aes(y = upper), size = 1.5, shape = 4) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Cluster", title = "Response risk", y = "Probability") +
        ggplot2::facet_grid(. ~ factor(category))
}

#' Plots the cluster sizes
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotClusterSizes <- function (...) {
    # riskprofs <- list(risk.profile.object, ...)
    # arg <- deparse(substitute(risk.profile.object))
    # dots <- substitute(list(...))[-1]
    # names <- c(arg, sapply(dots, deparse))
    # names(riskprofs) <- names

    riskprofs <- list(...)
    length(riskprofs) >= 1 ||
        stop("Supply at least one risk profile object.")
    # dots <- substitute(list(...))[-1]
    # names <- sapply(dots, deparse)
    # names(riskprofs) <- names

    data <- lapply(riskprofs, function(r) {
        orderStat <- apply(r$risk, 2, median)
        meanSortIndex <- order(orderStat, decreasing = F)

        nClusters <- r$riskProfClusObj$nClusters
        clusterSizes <-
            r$riskProfClusObj$clusterSizes[meanSortIndex]

        tibble::tibble(cluster = 1:nClusters,
                       clusterSize = clusterSizes)

    })
    data <- dplyr::bind_rows(data, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(
        x = factor(cluster),
        y = clusterSize,
        col = model
    )) +
        ggplot2::geom_point(size = 3,
                            position = ggplot2::position_dodge(width = .5)) +
        ggplot2::labs(title = "Cluster size",
                      y = "Number of subjects",
                      x = "Cluster") +
        ggplot2::expand_limits(y = 0)
}

#' Plots the similarity matrix, in the original order and sorted by
#'     (hierarchical) cluster to illustrate cluster structure.
#' @export
#' @param ... Object(s) of type \code{disSimObj}, output of
#'     \code{\link[PReMiuM]{calcDissimilarityMatrix}}.
plotSimilarityMatrix <- function(...) {
    # dissims <- list(dissim.object, ...)
    # arg <- deparse(substitute(dissim.object))
    # dots <- substitute(list(...))[-1]
    # names <- c(arg, sapply(dots, deparse))
    # names(dissims) <- names

    dissims <- list(...)
    length(dissims) >= 1 ||
        stop("Supply at least one dissimilarity object.")
    # dots <- substitute(list(...))[-1]
    # names <- sapply(dots, deparse)
    # names(dissims) <- names

    data <- lapply(dissims, function(d) {
        n <- d$disSimRunInfoObj$nSubjects

        dsmat <- matrix(0, nrow = n, ncol = n)
        dsmat[lower.tri(dsmat)] <- d$disSimMat
        dsmat[upper.tri(dsmat)] <- t(dsmat)[upper.tri(dsmat)]

        dsmat <- 1 - vec2mat(d$disSimMat, nrow = n)

        ddist <- as.dist(1 - dsmat)
        hc <- hclust(ddist)
        dsmat2 <- dsmat[hc$order, hc$order]

        dsmat[lower.tri(dsmat)] <- NA
        dsmat2[lower.tri(dsmat2)] <- NA

        m <- tibble::as_tibble(reshape2::melt(dsmat, na.rm = TRUE))
        m2 <-
            tibble::as_tibble(reshape2::melt(dsmat2, na.rm = TRUE))
        m$type <- "Subjects in original order"
        m2$type <- "Subjects reordered by hclust"
        m3 <- dplyr::bind_rows(m, m2)
    })
    data <- dplyr::bind_rows(data, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(model ~ type) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme(
            line = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        ) +
        ggplot2::scale_fill_gradient(
            name = "",
            low = "white",
            high = "black",
            limits = c(0, 1)
        ) +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Similarity matrix")
}

#' Plots the cumulative distribution of the variable selection
#'     parameter \code{rho}. Can be interpreted as the proportion of
#'     covariates that would be dropped if only those above a given
#'     threshold of \code{rho} were included.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotVarSelectRho <- function(...) {
    riskprofs <- list(...)
    length(riskprofs) >= 1 ||
        stop("Supply at least one risk profile object.")
    ## To do: add condition that names models if names not supplied:
    # dots <- substitute(list(...))[-1]
    # names <- sapply(dots, deparse)
    # names(models) <- names

    data <- lapply(riskprofs, function(m) {
        rho <- summariseVarSelectRhoFromRiskProfObj(m)$rhoMean
        ecd <- ecdf(rho)(rho) # ecdf(rho) returns a function!
        tibble::tibble(rho = rho, ecd = ecd)
    })
    data <- dplyr::bind_rows(data, .id = "model")
    ggplot2::ggplot(data, ggplot2::aes(x = rho, y = ecd, col = model)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "ECDF of rho",
                      x = "rho (mean for each variable)",
                      y = "proportion of variables deselected")
}

#' Assemble MCMC samples from \code{link[PReMiuM]{profRegr}} into
#'     a \code{\link[coda]{mcmc.list}} object with multiple chains.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A \code{\link[coda]{mcmc.list}}.
codaFromPremium <- function(global.parameter, ...) {
    # models <- list(premium.model, ...)
    # arg <- deparse(substitute(premium.model))
    # dots <- substitute(list(...))[-1]
    # names <- c(arg, sapply(dots, deparse))
    # names(models) <- names

    models <- list(...)
    length(models) >= 1 ||
        stop("Supply at least one premium model.")
    # dots <- substitute(list(...))[-1]
    # names <- sapply(dots, deparse)
    # names(models) <- names

    data <- lapply(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_", global.parameter, ".txt", sep = ""))
        parData <- coda::as.mcmc(read.table(parFileName)$V1)
    })
    coda::mcmc.list(data)
}

#' Get hyperparameters used to fit a PReMiuM model.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A matrix with a named column for each model and a named
#'     row for each hyperparameter.
getHyperparams <- function(...) {
    # Must be a cleaner way to do this...

    models <- list(...)
    length(models) >= 1 ||
        stop("Supply at least one premium model.")
    # dots <- substitute(list(...))[-1]
    # names <- sapply(dots, deparse)
    # names(models) <- names

    data <- sapply(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_hyper.txt", sep = ""))
        hypstr <- read.table(parFileName)$V1
        hypstr <- stringr::str_split(hypstr, "=")
        hyptab <- t(sapply(hypstr, function (h) {
            h
        }))
        hyp <- as.numeric(hyptab[, 2])
        names(hyp) <- hyptab[, 1]
        hyp
    })
    data
}

#' Render a standardised HTML report on PReMiuM model(s) using
#'     markdown, with covariate/response profiles and MCMC diagnostics.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @param filename File name (.html) to save markdown report.
#' @param rmd.template An alternate Rmd template to populate using
#'     the models in \code{...}.
#' @return An HTML file written to \code{filename}.
renderPremiumReport <-
    function(...,
             filename = "premium_report.html",
             rmd.template = NULL) {
        premium.models <- list(...)
        length(premium.models) >= 1 ||
            stop("Supply at least one premium model.")
        if (is.null(rmd.template)) {
            rmd.template <-
                system.file("rmd/premium_report.Rmd", package = "premiumPlots")
        }
        knitr.root.directory <- getwd()
        filepath <- file.path(knitr.root.directory, filename)
        rmarkdown::render(rmd.template, output_file = filepath)
    }
