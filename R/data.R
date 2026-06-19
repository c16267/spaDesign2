# =============================================================================
# data.R — documentation for the example datasets shipped under data/
# Provenance: all objects derive from the human dorsolateral prefrontal cortex
# (DLPFC) 10x Visium data of Maynard et al. (2021), obtained via the
# spatialLIBD package (Pardo et al., 2022). Genes and spots are subsampled so
# that the examples run quickly.
# =============================================================================


#' Mini DLPFC Pilot Data List
#'
#' A small, downsampled list of human dorsolateral prefrontal cortex (DLPFC)
#' spatial transcriptomics pilot samples, in the input format expected by
#' \code{\link{createSpaDesign2Object}}.
#'
#' @format A named \code{list} of pilot samples (DLPFC sections, e.g.
#'   \code{"151507"}--\code{"151509"}). Each element is itself a \code{list}
#'   with:
#'   \describe{
#'     \item{\code{counts}}{A genes \eqn{\times} spots count matrix (reduced
#'       gene panel).}
#'     \item{\code{coords}}{A \code{data.frame} of spot-level metadata with
#'       columns \code{x}, \code{y} (spatial coordinates) and \code{domain}
#'       (annotated spatial domain, e.g. \code{"Layer1"}--\code{"Layer6"} and
#'       \code{"WM"}).}
#'     \item{\code{group}}{Experimental group label (control \code{= 0},
#'       case \code{= 1}).}
#'     \item{\code{sample_id}}{Character sample identifier.}
#'   }
#' @keywords datasets
#' @seealso \code{\link{createSpaDesign2Object}}
#' @source Derived from the human DLPFC 10x Visium data of Maynard et al.
#'   (2021), \emph{Nature Neuroscience} 24:425--436, obtained via the
#'   \pkg{spatialLIBD} package (Pardo et al., 2022), \emph{BMC Genomics} 23:434.
"mini_pilot_data_list"


#' Mini spaDesign2 Object (Created, Unfitted)
#'
#' A \code{\linkS4class{spaDesign2}} object built from
#' \code{\link{mini_pilot_data_list}} via \code{\link{createSpaDesign2Object}}.
#' The \code{pilot_data} slot is populated, but no parameters have been
#' estimated yet, so \code{params_composition}, \code{params_geometry},
#' \code{params_expression}, \code{synthetic_data}, and \code{testing_result}
#' are still empty. It is the natural starting point for feature selection and
#' parameter estimation.
#'
#' @format An object of class \code{\linkS4class{spaDesign2}} with the
#'   \code{pilot_data} slot populated and the remaining slots empty.
#' @keywords datasets
#' @seealso \code{\link{createSpaDesign2Object}}, \code{\link{featureSelection}},
#'   \code{\link{mini_obj_features}}, \code{\link{mini_obj_fitted}}
#' @source See \code{\link{mini_pilot_data_list}}.
"mini_spaDesign2_obj"


#' Mini spaDesign2 Object (After Feature Selection)
#'
#' A \code{\linkS4class{spaDesign2}} object obtained from
#' \code{\link{mini_spaDesign2_obj}} after running \code{\link{featureSelection}}
#' and \code{\link{featureSelectionStable}}. Relative to the unfitted object,
#' its \code{params_expression} slot now carries the selected gene sets
#' (\code{top_genes}, \code{stable_genes}) together with the supporting
#' \code{stable_gene_stats} table; the generative parameters (composition,
#' geometry, expression) have not yet been fitted.
#'
#' @format An object of class \code{\linkS4class{spaDesign2}}. The
#'   \code{params_expression} slot contains \code{top_genes} (domain-informative
#'   markers), \code{stable_genes} (empirical-null / spike candidates), and
#'   \code{stable_gene_stats} (a \code{data.frame} of per-gene stability
#'   statistics).
#' @keywords datasets
#' @seealso \code{\link{featureSelection}}, \code{\link{featureSelectionStable}},
#'   \code{\link{makeCustomGeneSets}}, \code{\link{mini_obj_fitted}}
#' @source See \code{\link{mini_pilot_data_list}}.
"mini_obj_features"


#' Mini spaDesign2 Object (Fully Fitted)
#'
#' A fully fitted \code{\linkS4class{spaDesign2}} object obtained from
#' \code{\link{mini_obj_features}} after estimating all three generative layers
#' with \code{\link{estimateCompositionParams_v2}},
#' \code{\link{estimateGeometryParams_v2}}, and
#' \code{\link{estimateExpressionParams_v2}}. It is ready to drive synthetic-data
#' generation (\code{\link{simulateSpaDesign2}}) and power evaluation.
#'
#' @format An object of class \code{\linkS4class{spaDesign2}} with all parameter
#'   slots populated:
#'   \describe{
#'     \item{\code{params_composition}}{Baseline-anchored Binomial--logit fit,
#'       including \code{beta_binomial} with intercept \code{beta0} and effect
#'       \code{beta1}.}
#'     \item{\code{params_geometry}}{Fitted Fisher--Gaussian kernel mixture
#'       (FGKMM) geometry parameters.}
#'     \item{\code{params_expression}}{Gene sets plus the fitted expression
#'       parameters, e.g. domain means (\code{mu_grand}), spatial covariance
#'       (\code{theta}), and between-sample variance (\code{sigma_bio}).}
#'   }
#' @keywords datasets
#' @seealso \code{\link{estimateExpressionParams_v2}},
#'   \code{\link{simulateSpaDesign2}}, \code{\link{pBANSKY}}
#' @source See \code{\link{mini_pilot_data_list}}.
"mini_obj_fitted"


#' Mini Custom Gene Sets (WM vs. Layer6)
#'
#' The gene-set object returned by \code{\link{makeCustomGeneSets}} for the
#' Target/Reference pair \code{"WM"} (target) versus \code{"Layer6"} (reference),
#' built from the mini DLPFC pilot data. It supplies the labeling, empirical-null,
#' and signal-injection gene sets used by the simulation and testing endpoints.
#'
#' @format A named \code{list} with elements:
#'   \describe{
#'     \item{\code{G_svg}}{Character vector. Domain-informative + stable genes
#'       (labeling set).}
#'     \item{\code{G_null}}{Character vector. Stable / empirical-null genes.}
#'     \item{\code{G_spike}}{Character vector. Top differentially expressed genes
#'       used for signal injection.}
#'     \item{\code{TR_stats}}{A \code{data.frame} with columns \code{gene},
#'       \code{mean_T}, \code{mean_R}, \code{TR} (= \code{mean_T - mean_R}), and
#'       \code{abs_TR}.}
#'   }
#' @keywords datasets
#' @seealso \code{\link{makeCustomGeneSets}}, \code{\link{SaLFC}},
#'   \code{\link{simulateSpaDesign2}}
#' @source See \code{\link{mini_pilot_data_list}}.
"mini_custom_genes"
