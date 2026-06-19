# =============================================================================
# LOR: sample-level composition log-odds-ratio endpoint + pilot TREAT threshold.
#
# `%||%` is defined once in spaDesign2-package.R and is assumed available here.
# =============================================================================

#' @title Composition Test (K-driven, Sample-Level): Naive Two-Sample or TREAT
#' @description
#' Computes per-sample \eqn{L_k = \log((y_k + \epsilon)/(b_k + \epsilon))} and
#' compares Case vs Control at the sample level using a Welch t-test or a
#' Wilcoxon rank-sum test, for a single (target, reference) domain pair. Here
#' \eqn{y_k} and \eqn{b_k} are the spot counts of the target and reference
#' domains in sample \eqn{k}.
#'
#' Two null modes, controlled by \code{null_type}:
#' \describe{
#'   \item{\code{"treat"} (default)}{TREAT-style test (Smyth, 2009) centered on
#'     the pilot baseline \eqn{\delta_0 = \beta_1}.
#'     \itemize{
#'       \item \code{tau = 0} (centered null): \eqn{H_0:\ \hat\delta = \delta_0}.
#'       \item \code{tau > 0} (TREAT): \eqn{H_0:\ |\hat\delta - \delta_0| \le \tau}.
#'     }}
#'   \item{\code{"two_sample"}}{Naive two-sample test; the pilot baseline is
#'     ignored and \eqn{\delta_0 = 0}, i.e. \eqn{H_0:\ \bar\delta_1 = \bar\delta_0}.}
#' }
#'
#' The threshold \code{tau} may be a numeric scalar or the string \code{"pilot"},
#' in which case it is estimated from \code{object@pilot_data} via
#' \code{estimate_lor_tau_from_pilot()} for this single (target, reference) pair.
#'
#' @param processed_data List of samples. Each element is a list with
#'   \code{coords} (containing \code{hat_d} or \code{domain}) and \code{meta}
#'   (with \code{group} in (0/1) and optionally \code{sample_id}).
#' @param object A spaDesign2 object with \code{@params_composition} (and
#'   \code{@pilot_data} when \code{tau = "pilot"}).
#' @param target_domain Character. Target domain label.
#' @param reference_domain Character. Reference domain label.
#' @param test_method "t" (Welch, default) or "wilcoxon".
#' @param alternative "greater" (default), "two.sided", or "less".
#' @param null_type "treat" (default) or "two_sample".
#' @param tau Numeric scalar >= 0, or the string \code{"pilot"}. TREAT threshold
#'   (used only when \code{null_type = "treat"}); ignored (forced to 0) for
#'   "two_sample". If \code{"pilot"}, the threshold is estimated on the fly.
#' @param pilot_tau_args Named list of extra arguments forwarded to
#'   \code{estimate_lor_tau_from_pilot()} when \code{tau = "pilot"}
#'   (e.g. \code{list(q = 0.75, min_tau = 0.05)}). \code{object},
#'   \code{target_domain}, \code{reference_domain}, \code{test_method},
#'   \code{alternative}, and \code{eps} are supplied automatically.
#' @param eps Numeric > 0. Continuity correction (default 0.5).
#' @param verbose Logical (default TRUE).
#'
#' @return A named list with elements \code{null_type, target, reference, K0, K1,
#'   delta_hat, delta0, delta_tilde, tau, eps, Z_bar0, Z_bar1, y, m, sample_id,
#'   group, stat, df, p_value, method} (and \code{stat_treat} in TREAT mode).
#'
#' @references Smyth, G.K. (2009). Testing significance relative to a fold-change
#'   threshold is a TREAT. \emph{Bioinformatics}, 25(6), 765-771.
#'
#' @seealso \code{\link{evaluatePowerLOR}}, \code{\link{estimate_lor_tau_from_pilot}}
#'
#' @importFrom stats wilcox.test var pt
#' @export
#'
#' @examples
#' \dontrun{
#' data("mini_obj_fitted",   package = "spaDesign2")
#' data("mini_custom_genes", package = "spaDesign2")
#'
#' ## One-point LOR power calculation (single design point: K = 5,
#' ## case target proportion = 0.7). Expression DE is switched off.
#' sim <- simulateSpaDesign2(mini_obj_fitted, n_sample_per_group = 5,
#'   scenario_settings = list(DE_lfc = 0, target_prop_case = 0.7),
#'   target_domain = "WM", genes_to_simulate = mini_custom_genes$G_svg,
#'   de_genes = character(0), seed_base = 1)
#'
#' syn <- pBANSKY(sim, mini_obj_fitted@pilot_data, lambda = 0.5, do_hungarian = TRUE)
#' syn <- setNames(lapply(syn, function(s) list(coords = s$coords,
#'          meta = data.frame(sample_id = s$sample_id, group = s$group))),
#'        vapply(syn, function(s) s$sample_id, character(1)))
#'
#' res <- LOR(syn, mini_obj_fitted, target_domain = "WM",
#'   reference_domain = "Layer6", null_type = "treat", tau = "pilot")
#'
#' ## single-test rejection at this design point
#' res$p_value <= 0.05
#' }
LOR <- function(processed_data,
                object,
                target_domain    = "WM",
                reference_domain = "Layer6",
                test_method      = c("t", "wilcoxon"),
                alternative      = c("greater", "two.sided", "less"),
                null_type        = c("treat", "two_sample"),
                tau              = 0,
                pilot_tau_args   = list(),
                eps              = 0.5,
                verbose          = TRUE) {

  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)
  null_type   <- match.arg(null_type)

  # -- Input validation --------------------------------------------------------
  if (!is.list(processed_data) || length(processed_data) == 0)
    stop("processed_data must be a non-empty list.")
  if (is.null(object@params_composition) || length(object@params_composition) == 0)
    stop("object@params_composition is empty.")
  if (!is.finite(eps) || eps <= 0) stop("eps must be > 0.")

  # -- Optional: estimate tau from the pilot on the fly (single pair) ----------
  # The estimator runs LOR on pilot with a numeric tau = 0, so no recursion.
  if (is.character(tau) && length(tau) == 1L && tolower(tau) == "pilot") {
    est_defaults <- list(
      object = object, target_domain = target_domain,
      reference_domain = reference_domain, test_method = test_method,
      alternative = alternative, eps = eps, verbose = verbose,
      test_fn = sys.function()
    )
    est_args <- utils::modifyList(est_defaults, pilot_tau_args)
    tau <- do.call(estimate_lor_tau_from_pilot, est_args)
  }

  if (!is.finite(tau) || tau < 0) stop("tau must be >= 0 (or the string \"pilot\").")
  if (null_type == "two_sample" && tau > 0) {
    warning("tau > 0 is ignored when null_type = 'two_sample'. Setting tau = 0.")
    tau <- 0
  }

  # -- Pilot baseline delta0 ---------------------------------------------------
  # treat      : delta0 = pilot beta1 (null center)
  # two_sample : delta0 = 0           (H0: bar1 = bar0)
  if (null_type == "two_sample") {
    delta0 <- 0
    if (verbose) message(">>> [two_sample] delta0 = 0. Pilot baseline ignored.")
  } else {
    delta0 <- object@params_composition$beta_binomial$beta1
    if (!is.finite(delta0)) {
      warning("delta0 (beta1) is non-finite; defaulting to 0.")
      delta0 <- 0
    }
    if (verbose) message(sprintf(
      ">>> [treat] delta0 = %.4f, tau = %.4f  (%s)",
      delta0, tau, if (tau == 0) "centered null" else "TREAT threshold"))
  }

  # -- Per-sample L_k (vectorized) ---------------------------------------------
  nm <- names(processed_data)
  if (is.null(nm) || any(nm == "")) nm <- paste0("s", seq_along(processed_data))

  stats_mat <- vapply(seq_along(processed_data), function(i) {
    samp   <- processed_data[[i]]
    na_out <- c(L_k = NA_real_, grp = NA_real_, m_k = NA_real_, y_k = NA_real_)

    if (is.null(samp$coords)) return(na_out)
    preds <- as.character(samp$coords$hat_d %||% samp$coords$domain)
    if (is.null(preds)) return(na_out)

    g <- NA_real_
    if (!is.null(samp$meta) && "group" %in% names(samp$meta))
      g <- as.numeric(samp$meta$group)
    if (!is.finite(g) || !(g %in% c(0, 1))) return(na_out)

    y_k <- sum(preds == target_domain,    na.rm = TRUE)
    b_k <- sum(preds == reference_domain, na.rm = TRUE)
    m_k <- y_k + b_k
    if (!is.finite(m_k) || m_k <= 0) return(na_out)

    c(L_k = log(y_k + eps) - log(b_k + eps), grp = g, m_k = m_k, y_k = y_k)
  }, numeric(4))

  # -- Filter valid samples ----------------------------------------------------
  valid <- is.finite(stats_mat["L_k", ]) & is.finite(stats_mat["grp", ])
  if (!any(valid)) return(list(p_value = NA_real_, msg = "No valid samples."))

  Z     <- stats_mat["L_k", valid]
  grp   <- stats_mat["grp", valid]
  m_vec <- stats_mat["m_k", valid]
  y_vec <- stats_mat["y_k", valid]

  sid <- vapply(seq_along(processed_data)[valid], function(i) {
    samp <- processed_data[[i]]
    if (!is.null(samp$meta) && "sample_id" %in% names(samp$meta))
      return(as.character(samp$meta$sample_id))
    nm[i]
  }, character(1))

  if (length(unique(grp)) < 2)
    return(list(p_value = NA_real_, msg = "Only one group present."))

  Z0 <- Z[grp == 0]; Z1 <- Z[grp == 1]
  K0 <- length(Z0);  K1 <- length(Z1)
  if (K0 < 1 || K1 < 1)
    return(list(p_value = NA_real_, msg = "Insufficient samples per group."))

  bar0        <- mean(Z0); bar1 <- mean(Z1)
  delta_hat   <- bar1 - bar0
  delta_tilde <- delta_hat - delta0  # two_sample: delta0 = 0 -> delta_tilde = delta_hat

  out <- list(
    null_type   = null_type,
    target      = target_domain,    reference   = reference_domain,
    K0          = K0,               K1          = K1,
    delta_hat   = delta_hat,        delta0      = delta0,
    delta_tilde = delta_tilde,      tau         = tau,
    eps         = eps,
    Z_bar0      = bar0,             Z_bar1      = bar1,
    y           = y_vec,            m           = m_vec,
    sample_id   = sid,              group       = grp
  )

  # ===========================================================================
  # (i) Welch t-test
  # ===========================================================================
  if (test_method == "t") {
    s0 <- if (K0 >= 2) stats::var(Z0) else NA_real_
    s1 <- if (K1 >= 2) stats::var(Z1) else NA_real_

    if (!is.finite(s0) || !is.finite(s1)) {
      out$stat <- out$df <- out$p_value <- NA_real_
      out$method <- "welch_t"
      out$msg    <- "Welch t-test requires >= 2 samples per group."
      return(out)
    }

    se <- sqrt(s1 / K1 + s0 / K0)
    df <- (s1/K1 + s0/K0)^2 / ((s1/K1)^2/(K1 - 1) + (s0/K0)^2/(K0 - 1))

    if (!is.finite(se) || se <= 0 || !is.finite(df) || df <= 0) {
      out$stat <- out$df <- out$p_value <- NA_real_
      out$method <- "welch_t"
      out$msg    <- "Non-finite SE or df."
      return(out)
    }

    # t_stat = delta_tilde / se  (treat: (delta_hat-delta0)/se; two_sample: delta_hat/se)
    t_stat <- delta_tilde / se

    if (null_type == "treat" && tau > 0) {
      if (alternative == "two.sided") {
        t_shift    <- (abs(delta_tilde) - tau) / se
        p_val      <- 2 * stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_two_sided"
      } else if (alternative == "greater") {
        t_shift    <- (delta_tilde - tau) / se
        p_val      <- stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_greater"
      } else {
        t_shift    <- (-delta_tilde - tau) / se
        p_val      <- stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_less"
      }
      out$stat_treat <- t_shift
      out$p_value    <- min(1, p_val)
    } else {
      if (alternative == "two.sided") {
        p_val <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
      } else if (alternative == "greater") {
        p_val <- stats::pt(t_stat, df = df, lower.tail = FALSE)
      } else {
        p_val <- stats::pt(t_stat, df = df, lower.tail = TRUE)
      }
      out$p_value <- p_val
      out$method  <- if (null_type == "two_sample")
        paste0("welch_t_two_sample_", alternative)
      else
        paste0("welch_t_treat_centered_", alternative)  # tau = 0
    }

    out$stat <- t_stat
    out$df   <- df

    if (verbose) message(sprintf(
      ">>> Composition K-test (t) [%s]: alt=%s, delta_hat=%.4f, delta0=%.4f, delta_tilde=%.4f, tau=%.4f, p=%.3e",
      null_type, alternative, delta_hat, delta0, delta_tilde, tau, out$p_value))
    return(out)
  }

  # ===========================================================================
  # (ii) Wilcoxon rank-sum
  # ===========================================================================
  if (test_method == "wilcoxon") {
    # mu: two_sample -> 0; treat tau=0 -> delta0; treat tau>0 -> delta0 +/- tau
    mu_base <- delta0

    if (null_type == "treat" && tau > 0) {
      if (alternative == "two.sided") {
        p_one <- if (delta_tilde >= 0) {
          tryCatch(stats::wilcox.test(Z1, Z0, alternative = "greater",
                                      mu = mu_base + tau, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        } else {
          tryCatch(stats::wilcox.test(Z1, Z0, alternative = "less",
                                      mu = mu_base - tau, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        }
        p_val      <- if (is.finite(p_one)) min(1, 2 * p_one) else NA_real_
        out$method <- "wilcoxon_treat_two_sided"
      } else if (alternative == "greater") {
        p_val      <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = "greater",
                                                  mu = mu_base + tau, exact = FALSE)$p.value,
                               error = function(e) NA_real_)
        out$method <- "wilcoxon_treat_greater"
      } else {
        p_val      <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = "less",
                                                  mu = mu_base - tau, exact = FALSE)$p.value,
                               error = function(e) NA_real_)
        out$method <- "wilcoxon_treat_less"
      }
    } else {
      p_val <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = alternative,
                                           mu = mu_base, exact = FALSE)$p.value,
                        error = function(e) NA_real_)
      out$method <- if (null_type == "two_sample")
        paste0("wilcoxon_two_sample_", alternative)
      else
        paste0("wilcoxon_treat_centered_", alternative)  # tau = 0
    }

    out$stat    <- NA_real_
    out$p_value <- p_val

    if (verbose) message(sprintf(
      ">>> Composition K-test (wilcoxon) [%s]: alt=%s, delta_hat=%.4f, delta0=%.4f, delta_tilde=%.4f, tau=%.4f, p=%.3e",
      null_type, alternative, delta_hat, delta0, delta_tilde, tau, out$p_value))
    return(out)
  }

  out
}


#' @title Estimate the Composition TREAT Threshold (tau) from Pilot Data
#' @description
#' Fast plug-in estimator of the scalar TREAT threshold \eqn{\tau} for the
#' sample-level composition test \code{\link{LOR}}, mirroring the SaLFC pilot
#' estimator but on the log-odds-ratio scale. \code{LOR} is run once on the pilot
#' data with \code{tau = 0} (centered null); \eqn{\tau} is then defined
#' direction-aware from the pilot effect and its uncertainty:
#' \deqn{\tau = \max(\mathrm{effect\_ref},\; z_q \cdot \mathrm{se}),}
#' where \eqn{\mathrm{effect\_ref}} depends on \code{alternative}, the effect is
#' the de-baselined \eqn{\tilde\delta = \hat\delta - \delta_0}, and
#' \eqn{z_q = \Phi^{-1}(q)}. Returns one scalar per (target, reference) pair.
#'
#' Normally invoked automatically by \code{LOR(..., tau = "pilot")}; can also be
#' called directly. The pilot run uses the ground-truth \code{coords$domain}
#' (pilot data has no \code{hat_d}).
#'
#' @param object A spaDesign2 object with \code{@pilot_data} and
#'   \code{@params_composition}.
#' @param target_domain,reference_domain Character domain labels.
#' @param test_method "t" (default) or "wilcoxon". Only "t" yields an se for the
#'   plug-in; "wilcoxon" falls back to \code{effect_ref} alone (se term = 0).
#' @param alternative "greater" (default), "less", or "two.sided".
#' @param q Numeric in (0,1). Standard-normal quantile level for \eqn{z_q}.
#' @param min_tau Numeric. Lower bound for \eqn{\tau}.
#' @param cap_effect Numeric or NULL. Optional upper cap on \eqn{\tau}.
#' @param eps Numeric > 0. Continuity correction forwarded to \code{LOR}.
#' @param verbose Logical.
#' @param test_fn Function. The composition test to run on the pilot. If NULL,
#'   resolves \code{LOR} by name. When invoked via \code{LOR(..., tau = "pilot")},
#'   the running test function is passed automatically.
#'
#' @return A numeric scalar \eqn{\tau} (>= \code{min_tau}).
#'
#' @seealso \code{\link{LOR}}
#'
#' @importFrom stats qnorm
#' @export
#'
estimate_lor_tau_from_pilot <- function(object,
                                        target_domain    = "WM",
                                        reference_domain = "Layer6",
                                        test_method      = c("t", "wilcoxon"),
                                        alternative      = c("greater", "less", "two.sided"),
                                        q                = 0.75,
                                        min_tau          = 0,
                                        cap_effect       = NULL,
                                        eps              = 0.5,
                                        verbose          = TRUE,
                                        test_fn          = NULL) {
  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)

  if (is.null(object@pilot_data) || length(object@pilot_data) == 0)
    stop("object@pilot_data is empty.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0,1).")
  if (!is.finite(min_tau) || min_tau < 0) stop("min_tau must be >= 0.")

  if (is.null(test_fn)) {
    if (exists("LOR", mode = "function")) {
      test_fn <- match.fun("LOR")
    } else if (exists("test_composition_K_treat", mode = "function")) {
      test_fn <- match.fun("test_composition_K_treat")
    } else {
      stop("No composition test function found; pass it via `test_fn`.")
    }
  }
  if (!is.function(test_fn)) stop("`test_fn` must be a function.")

  # 1) Build pilot processed_data (ground-truth domains; LOR uses domain via %||%)
  pilot_pd <- lapply(seq_along(object@pilot_data), function(i) {
    s <- object@pilot_data[[i]]
    if (is.null(s$coords) || is.null(s$group)) return(NULL)
    if (!("domain" %in% names(s$coords))) return(NULL)
    sid <- if (!is.null(names(object@pilot_data)) && names(object@pilot_data)[i] != "")
      names(object@pilot_data)[i] else paste0("pilot_", i)
    list(coords = s$coords,
         meta   = data.frame(sample_id = sid, group = as.numeric(s$group),
                             stringsAsFactors = FALSE))
  })
  pilot_pd <- Filter(Negate(is.null), pilot_pd)
  if (length(pilot_pd) == 0) stop("No valid pilot samples with coords/group.")
  names(pilot_pd) <- vapply(pilot_pd, function(x) x$meta$sample_id, character(1))

  # 2) Run the composition test on the pilot with tau = 0 (centered null)
  pilot_res <- tryCatch(
    test_fn(processed_data   = pilot_pd,
            object           = object,
            target_domain    = target_domain,
            reference_domain = reference_domain,
            test_method      = test_method,
            alternative      = alternative,
            null_type        = "treat",
            tau              = 0,
            eps              = eps,
            verbose          = FALSE),
    error = function(e) { message("\n[ERROR in LOR pilot] ", e$message); NULL })

  if (is.null(pilot_res) || is.null(pilot_res$delta_tilde) ||
      !is.finite(pilot_res$delta_tilde)) {
    if (verbose) message(">>> [pilot tau] pilot effect unavailable -> tau = min_tau")
    return(max(min_tau, 0))
  }

  # 3) Direction-aware plug-in tau (scalar)
  delta_tilde <- pilot_res$delta_tilde
  se <- NA_real_
  if (test_method == "t" && !is.null(pilot_res$stat) &&
      is.finite(pilot_res$stat) && pilot_res$stat != 0) {
    se <- abs(delta_tilde / pilot_res$stat)  # t_stat = delta_tilde / se
  }
  if (!is.finite(se) || se < 0) se <- 0

  z_q <- stats::qnorm(q); if (!is.finite(z_q)) z_q <- 0; z_q <- max(z_q, 0)

  effect_ref <- switch(alternative,
                       greater   = max(delta_tilde, 0),
                       less      = max(-delta_tilde, 0),
                       two.sided = abs(delta_tilde))
  if (!is.finite(effect_ref)) effect_ref <- 0

  tau <- max(effect_ref, z_q * se)
  if (!is.finite(tau)) tau <- effect_ref
  if (!is.finite(tau)) tau <- 0
  tau <- max(tau, min_tau)
  if (!is.null(cap_effect) && is.finite(cap_effect) && cap_effect > 0)
    tau <- min(tau, cap_effect)

  if (verbose) message(sprintf(
    ">>> [pilot tau] target=%s ref=%s alt=%s q=%.2f (z_q=%.3f): delta_tilde=%.4f, se=%.4f -> tau=%.4f",
    target_domain, reference_domain, alternative, q, z_q, delta_tilde, se, tau))

  tau
}
