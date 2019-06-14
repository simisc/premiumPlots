#' Pipe
#'
#' Imported from \code{\link[magrittr]{pipe}}
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

#' Summarise \code{rho}
#'
#' Summarise the posterior distribution of the variable selection
#'     parameter \code{rho}. Modified version of
#'     \code{\link[PReMiuM]{summariseVarSelectRho}}, allowing plotting
#'     functions to access \code{rho} summaries.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param rho_file Path to the rho.txt file written by \code{\link[PReMiuM]{profRegr}}.
#'     If NULL (the default), path is taken from the \code{ClusObjRunInfoObj} element
#'     of \code{riskProfObj}.
#' @return A \code{\link[tibble]{tibble}} containing the mean, median,
#'     lower (0.05) and upper (0.95) quantiles for each covariate,
#'     and a ranking from the covariate with the highest value of
#'     \code{rho} to the lowest.
tabulateVarSelectRho <- function (riskProfObj, rho_file = NULL) {

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    if (is.null(rho_file)) {
        # Allows rho_file to be specified in call (e.g. if WD different to the one in which profRegr was called)
        # Propagate this option to functions that call tabulateVarSelectRho()
        rho_file <- file.path(directoryPath, paste(fileStem, "_rho.txt", sep = ""))
    }

    rhoMat <- matrix(scan(rho_file, what = double(), quiet = T),
                     ncol = nCovariates,
                     byrow = T)
    firstLine <- ifelse(reportBurnIn, nBurn / nFilter + 2, 1)
    lastLine <-
        (nSweeps + ifelse(reportBurnIn, nBurn + 1, 0)) / nFilter
    rhoMat <- rhoMat[firstLine:lastLine, ]
    rhoMean <- apply(rhoMat, 2, mean)
    rhoMedian <- apply(rhoMat, 2, stats::median)
    rhoLowerCI <- apply(rhoMat, 2, stats::quantile, 0.05)
    rhoUpperCI <- apply(rhoMat, 2, stats::quantile, 0.95)

    tibble::tibble(
        covname = covNames,
        rhoMean = rhoMean,
        rhoMedian = rhoMedian,
        rhoLowerCI = rhoLowerCI,
        rhoUpperCI = rhoUpperCI
        ) %>%
        dplyr::mutate(rhoRank = rank(-rhoMean)) %>%
        dplyr::arrange(rhoRank)
}

#' Plot covariate profiles by cluster
#'
#' Plots the (discrete) covariate profiles, facetted vertically by cluster and
#'     (optionally) horizontally by some grouping (\code{split}) of the covariates.
#'     Covariate categories that are more common within a cluster than they are
#'     overall are shaded dark; those that are not are shaded light.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed.
#' @param rhoOrder Should covariates be ordered in descending order of \code{rho}?
#'     If \code{TRUE}, the default, then indices given in \code{whichCovariates}
#'     will be interpreted as indices of the reordered covariates.
#'     Ignored silently if variable selection was not used.
#' @param useProfileStar The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#'     Ignored silently if variable selection was not used.
#' @param covariate_info Optional list of details about the covariates,
#'     with (some of) the following named elements:
#'     \describe{
#'         \item{\code{title}:}{String describing the covariate, default
#'         to "Covariate category".}
#'         \item{\code{levels}:}{Vector of integer values taken by the
#'         (discrete) covariates, specifying the order of the facets.}
#'         \item{\code{labels}:}{Vector of strings giving labels for each
#'         level of the covariate, in the same order as \code{levels}.}
#'         \item{\code{split}:}{A character vector specifying the covariate
#'         grouping for horizontal facetting. The first element appears
#'         in a subset of covariates to be plotted in a separate horizontal
#'         facet. The second and third elements, if present, are facet labels
#'         for covariates that do and do not meet this criterion, respectively,
#'         e.g. \code{c("_citr", "Citrullinated", "Not citrullinated")}.}
#'     }
plotProfilesByCluster <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoOrder = TRUE,
                                   useProfileStar = TRUE,
                                   covariate_info = list(title = "Covariate\ncategory",
                                                         levels = NULL,
                                                         labels = NULL,
                                                         split = NULL)
                                   ) {

    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoOrder = rhoOrder,
        useProfileStar = useProfileStar
    )

    if (!is.null(covariate_info$levels) & !is.null(covariate_info$labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = covariate_info$levels,
                                            labels = covariate_info$labels))
    }

    empirical_proportions <- profileDF %>%
        dplyr::group_by(category, covname) %>%
        dplyr::summarise(x = mean(mean)) %>%
        dplyr::group_by(covname) %>%
        dplyr::arrange(dplyr::desc(category)) %>%
        dplyr::mutate(emp_propn = cumsum(x)) %>%
        dplyr::filter(emp_propn < max(emp_propn)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-x)

    covtab <- profileDF %>%
        dplyr::left_join(empirical_proportions, by = c("covname", "category")) %>%
        dplyr::group_by(cluster, category, covname, covOrder, fillColor, emp_propn) %>%
        dplyr::summarise(prop = mean(est)) %>%
        dplyr::ungroup()

    if (!is.null(covariate_info$split)) {
        if(length(covariate_info$split) < 3) {
            covariate_info$split <- rep(covariate_info$split, 3 - length(covariate_info$split))
            covariate_info$split <- c(covariate_info$split, sprintf("not %s", covariate_info$split[2]))
        }
        covtab <- covtab %>%
            dplyr::mutate(type = stringr::str_detect(covname, covariate_info$split[1]),
                          type = dplyr::recode(
                              as.numeric(type), `1` = covariate_info$split[2], `0` = covariate_info$split[3]
                          ))
        facetting_layer <- list(
            ggplot2::facet_grid(cluster ~ type, scales = "free_x", space = "free_x")
        )
    } else {
        facetting_layer <- list(ggplot2::facet_grid(cluster ~ .))
    }

    ggplot2::ggplot(covtab,
                    ggplot2::aes(
                        x = stats::reorder(covname, covOrder),
                        y = prop,
                        fill = factor(category),
                        alpha = fillColor == "high"
                    )) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::geom_point(ggplot2::aes(y = emp_propn, group = category),
                            col = "black", fill = "white", alpha = 1, shape = 18, na.rm = TRUE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )) +
        ggplot2::labs(
            x = "Covariate",
            y = "Proportion (by cluster)",
            title = "Covariate profiles",
            subtitle = "Black pips: empirical proportions in each covariate category.\nDark fill: category more prevalent within cluster than overall.") +
        ggplot2::scale_fill_discrete(name = covariate_info$title) +
        ggplot2::scale_alpha_discrete(guide = FALSE,
                                      range = c(0.25, 1)) +
        facetting_layer
}

#' Plot covariate profiles by covariate
#'
#' Plots the (discrete) covariate profiles, facetted by covariate,
#'     similar to the layout in \code{\link[PReMiuM]{plotRiskProfile}}.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed.
#' @param rhoOrder Should covariates be ordered in descending order of \code{rho}?
#'     If \code{TRUE}, the default, then indices given in \code{whichCovariates}
#'     will be interpreted as indices of the reordered covariates.
#'     Ignored silently if variable selection was not used.
#' @param useProfileStar The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#'     Ignored silently if variable selection was not used.
#' @param covariate_levels Vector of integer values taken by the
#'     (discrete) covariates, specifying the order of the facets.
#' @param covariate_labels Vector of strings giving labels for each
#'     level of the covariate, in the same order as \code{covariate_levels}.
plotCovariateProfiles <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoOrder = TRUE,
                                   useProfileStar = TRUE,
                                   covariate_levels = NULL,
                                   covariate_labels = NULL) {
    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoOrder = rhoOrder,
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
        ggplot2::facet_grid(factor(category) ~ stats::reorder(covname, covOrder))
}

#' Make table of covariate profiles with info for plotting
#'
#' Used by \code{\link{plotProfilesByCluster}} and \code{\link{plotCovariateProfiles}}
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed.
#' @param rhoOrder Should covariates be ordered in descending order of \code{rho}?
#'     If \code{TRUE}, the default, then indices given in \code{whichCovariates}
#'     will be interpreted as indices of the reordered covariates.
#'     Ignored silently if variable selection was not used.
#' @param useProfileStar The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#'     Ignored silently if variable selection was not used.
#' @param rho_file Path to the rho.txt file written by \code{\link[PReMiuM]{profRegr}}.
#'     If NULL (the default), path is taken from the \code{ClusObjRunInfoObj} element
#'     of \code{riskProfObj}.
tabulateCovariateProfiles <- function (riskProfObj,
                                       whichCovariates = NULL,
                                       rhoOrder = TRUE,
                                       useProfileStar = TRUE,
                                       rho_file = NULL) {
    # This now does most of the work for plotCovariateProfiles
    # Use same output for new covariate profiles plot

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
        stop("Only implemented for discrete covariates.\nUse plotRiskProfile() for Normal- or Mixed-covariate models.")

    if (varSelect) {

        rhotab <- tabulateVarSelectRho(riskProfObj, rho_file = rho_file)

        if (useProfileStar) {
            profile <- profileStar
        }

        if(rhoOrder) {

            if (is.numeric(whichCovariates)) {
                whichCovariates <- rhotab %>%
                    dplyr::filter(rhoRank %in% whichCovariates) %>%
                    "$"(covname)
            } else if (is.character(whichCovariates)) {
                whichCovariates <- rhotab %>%
                    dplyr::filter(covname %in% whichCovariates) %>%
                    "$"(covname)
            } else {
                whichCovariates <- rhotab %>%
                    "$"(covname)
            }

        }

    }

    if(!is.null(whichCovariates)) {

        if (is.character(whichCovariates)) {
            whichCovariates <- match(whichCovariates, covNames)
        }

        covNames <- covNames[whichCovariates]
        nCovariates <- length(whichCovariates)
        profile <- profile[, , whichCovariates, , drop = FALSE]
        nCategories <- nCategories[whichCovariates]

    }

    orderStat <- apply(risk, 2, stats::median)
    sortIndex <- order(clusterSizes, orderStat, decreasing = T)
    clusterSizes <- clusterSizes[sortIndex]
    profile <- profile[, sortIndex, , , drop = FALSE]

    covOrder <- tibble::tibble(
        covname = covNames,
        covOrder = 1:nCovariates
    )

    purrr::map_dfr(1:nCovariates, function(j) {

        purrr::map_dfr(1:nCategories[j], function(k) {

            probMat <- profile[, , j, k]
            probMeans <- apply(probMat, 2, mean)

            probMean <- sum(probMeans * clusterSizes) / sum(clusterSizes)
            probLower <- apply(probMat, 2, stats::quantile, 0.05)
            probUpper <- apply(probMat, 2, stats::quantile, 0.95)

            clusterDF <- tibble::tibble(
                cluster = 1:nClusters,
                category = k - 1,
                mean = probMean,
                lower = probLower,
                upper = probUpper
            )

            tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                           est = c(probMat)) %>%
                dplyr::left_join(clusterDF, by = "cluster")
            }) %>%
            dplyr::mutate(covname = covNames[j])
        }) %>%
        dplyr::left_join(covOrder, by = "covname") %>%
        dplyr::mutate(
            fillColor = dplyr::case_when(
                lower > mean ~ "high",
                upper < mean ~ "low",
                TRUE ~ "avg"))
}

#' Plot response profiles
#'
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

    orderStat <- apply(risk, 2, stats::median)
    sortIndex <- order(clusterSizes, orderStat, decreasing = T)
    clusterSizes <- clusterSizes[sortIndex]
    risk <- risk[, sortIndex, ]

    profileDF <- purrr::map_dfr(1:nCategoriesY, function(k) {

        probMat <- risk[, , k] # 2D matrix for categ. k
        probMeans <- apply(probMat, 2, mean, trim = 0.005)

        probMean <- sum(probMeans * clusterSizes) / sum(clusterSizes)
        probLower <- apply(probMat, 2, stats::quantile, 0.05)
        probUpper <- apply(probMat, 2, stats::quantile, 0.95)

        clusterDF <- tibble::tibble(
            cluster = 1:nClusters,
            category = k - 1,
            mean = probMean,
            lower = probLower,
            upper = probUpper
        )

        tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                       est = c(probMat)) %>%
            dplyr::left_join(clusterDF, by = "cluster")
        }) %>%
        dplyr::mutate(
            fillColor = dplyr::case_when(
                lower > mean ~ "high",
                upper < mean ~ "low",
                TRUE ~ "avg"))

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
        ggplot2::labs(x = "Cluster", title = "Response profile", y = "Probability") +
        ggplot2::facet_grid(. ~ factor(category))
}

#' Plot cluster sizes
#'
#' Plots the cluster sizes from one or more PReMiuM models
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotClusterSizes <- function (...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- purrr::map_dfr(riskprofs, function(r) {

        nClusters <- r$riskProfClusObj$nClusters
        clusterSizes <- r$riskProfClusObj$clusterSizes

        orderStat <- apply(r$risk, 2, stats::median)
        sortIndex <- order(clusterSizes, orderStat, decreasing = T)
        clusterSizes <- clusterSizes[sortIndex]

        tibble::tibble(cluster = 1:nClusters,
                       clusterSize = clusterSizes)

        }, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(x = factor(cluster), y = clusterSize, col = model)) +
        ggplot2::geom_point(size = 3,
                            position = ggplot2::position_dodge(width = .5)) +
        ggplot2::labs(title = "Cluster size",
                      y = "Number of subjects",
                      x = "Cluster") +
        ggplot2::expand_limits(y = 0) +
        ggplot2::scale_colour_discrete(guide = length(riskprofs) > 1)
}

#' Plot similarity matrix
#'
#' Plots the similarity matrix, in the original order and sorted by
#'     (hierarchical) cluster to illustrate cluster structure.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotSimilarityMatrix <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- purrr::map_dfr(riskprofs, function(d) {

        disSimMat <- d$riskProfClusObj$clusObjDisSimMat

        n <- d$riskProfClusObj$clusObjRunInfoObj$nSubjects

        dsmat <- matrix(0, nrow = n, ncol = n)
        dsmat[lower.tri(dsmat)] <- disSimMat
        dsmat[upper.tri(dsmat)] <- t(dsmat)[upper.tri(dsmat)]

        dsmat <- 1 - vec2mat(disSimMat, nrow = n)

        ddist <- stats::as.dist(1 - dsmat)
        hc <- stats::hclust(ddist)
        dsmat2 <- dsmat[hc$order, hc$order]

        dsmat[lower.tri(dsmat)] <- NA
        dsmat2[lower.tri(dsmat2)] <- NA

        m <- tibble::as_tibble(reshape2::melt(dsmat, na.rm = TRUE))
        m2 <- tibble::as_tibble(reshape2::melt(dsmat2, na.rm = TRUE))
        m$type <- "Subjects in original order"
        m2$type <- "Subjects reordered by hclust"
        m3 <- dplyr::bind_rows(m, m2)
        }, .id = "model")

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

#' Plot cumulative distribution of \code{rho}
#'
#' Plots the cumulative distribution of the variable selection
#'     parameter \code{rho}. Can be interpreted as the proportion of
#'     covariates that would be dropped if only those above a given
#'     threshold of \code{rho} were included.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotVarSelectRho <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- purrr::map_dfr(riskprofs, function(m) {
        rho <- tabulateVarSelectRho(m)$rhoMean # Add option for rho_file location
        ecd <- stats::ecdf(rho)(rho) # ecdf(rho) returns a function!
        tibble::tibble(rho = rho, ecd = ecd)
        }, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(x = rho, y = ecd, col = model)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "ECDF of rho",
                      x = "rho (mean for each covariate)",
                      y = "cumulative proportion of covariates")
}

#' Plot distribution of \code{rho} per covariate
#'
#' Plots the mean, median and CI of the variable selection
#'     parameter \code{rho} for each covariate.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotRhoDistributions <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- purrr::map_dfr(riskprofs, tabulateVarSelectRho, .id = "model") %>%
        dplyr::rename(Mean = rhoMean, Median = rhoMedian) %>%
        tidyr::gather(centre, value, Mean, Median)

    ggplot2::ggplot(data, ggplot2::aes(x = factor(covname), col = model, group = model)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = rhoLowerCI, ymax = rhoUpperCI), position = ggplot2::position_dodge(width = 0.25)) +
        ggplot2::geom_point(ggplot2::aes(y = value, shape = centre), position = ggplot2::position_dodge(width = 0.25)) +
        ggplot2::labs(x = "Covariate", y = "rho") +
        ggplot2::scale_shape_manual(name = NULL,
                                    values = c(Mean = 16, Median = 4))
}

#' Make a coda object PReMiuM samples
#'
#' Assemble MCMC samples from \code{\link[PReMiuM]{profRegr}} into
#'     a \code{\link[coda]{mcmc.list}} object with multiple chains.
#'     Useful for convergence diagnositics.
#' @export
#' @param global_parameter Parameter to convert
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A \code{\link[coda]{mcmc.list}}.
codaFromPremium <- function(global_parameter, ...) {

    models <- list(...)
    if (is.null(names(models))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(models) <- names
    }

    purrr::map(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_", global_parameter, ".txt", sep = ""))
        parData <- coda::as.mcmc(utils::read.table(parFileName)$V1)
        }) %>%
        coda::mcmc.list()
}

#' Get PReMiuM hyperparameters
#'
#' Get hyperparameters used to fit a PReMiuM model.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A matrix with a named column for each model and a named
#'     row for each hyperparameter.
getHyperparams <- function(...) {

    # Must be a cleaner way to do this... hyperparameters are only
    # written to file if they differ from the defaults...

    models <- list(...)
    if (is.null(names(models))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(models) <- names
    }

    data <- sapply(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_hyper.txt", sep = ""))
        hypstr <- utils::read.table(parFileName)$V1
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

#' Render PReMiuM report
#'
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
        if (is.null(names(premium.models))) {
            # Name models if ... arguments are not supplied with names
            dots <- substitute(list(...))[-1]
            names <- sapply(dots, deparse)
            names(premium.models) <- names
        }

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

# to appease R CMD check
utils::globalVariables(c("Var1", "Var2", "category", "centre", "clusObjRunInfoObj", "cluster", "clusterSize",
                         "covNames", "covOrder", "covname", "directoryPath", "ecd", "emp_propn", "est", "fileStem",
                         "fillColor", "includeResponse", "lower", "model", "nBurn", "nCategoriesY",
                         "nClusters", "nCovariates", "nFilter", "nSweeps", "profileStar", "prop", "reportBurnIn",
                         "rho", "rhoLowerCI", "rhoMean", "rhoMedian", "rhoRank", "rhoUpperCI", "risk",
                         "riskProfClusObj", "type", "upper", "value", "varSelect",
                         "vec2mat", "x", "xModel", "yModel", "Mean", "Median"))
