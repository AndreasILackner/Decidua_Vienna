library('glmmTMB')
library('DHARMa')

# load data
d <- read.csv("./Xenium/core_level_for_glmm.csv", stringsAsFactors = FALSE)

# factors and reference levels
d$donor     <- factor(d$donor)
d$condition <- factor(d$condition)

# condition decP corresponds to A, decB to B
d$condition <- relevel(d$condition, ref = "A")

# area must be > 0 for log-offset
d <- subset(d, is.finite(area) & area > 0)
stopifnot(nrow(d) > 0)

# NB-GLMM with donor random intercept and log(area) offset
m_full <- glmmTMB(hubs ~ condition + offset(log(area)) + (1|donor),
                  family = nbinom2(), data = d)
m_red  <- update(m_full, . ~ . - condition)

# asymptotic LRT (for reference)
print(anova(m_red, m_full))

co <- summary(m_full)$coefficients$cond
est <- co["conditionB","Estimate"]; se <- co["conditionB","Std. Error"]
irr    <- exp(est)
irr_ci <- exp(est + qnorm(c(0.025, 0.975))*se)
cat(sprintf("IRR (B vs A): %.3f  95%% CI [%.3f, %.3f]\n", irr, irr_ci[1], irr_ci[2]))

# parametric bootstrap LRT under H0 (reduced model)
pb_lrt <- function(m_red, m_full, B = 2000, seed = 123){
  set.seed(seed)
  lr_obs <- 2*(as.numeric(logLik(m_full)) - as.numeric(logLik(m_red)))
  ysims  <- simulate(m_red, nsim = B, seed = seed)
  exceed <- 0L
  for (b in seq_len(B)) {
    d$hubs <- ysims[[b]]
    m0 <- suppressWarnings(update(m_red,  data = d))
    m1 <- suppressWarnings(update(m_full, data = d))
    lr <- 2*(as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
    if (is.finite(lr) && lr >= lr_obs) exceed <- exceed + 1L
  }
  data.frame(
    lr_obs = lr_obs, B = B, exceed = exceed,
    p_value = (exceed + 1)/(B + 1)  # unbiased finite-sample estimate
  )
}

boot_res <- pb_lrt(m_red, m_full, B = 2000, seed = 123)
print(boot_res)

# effect size: IRR for B vs A (with Wald 95% CI)
co <- summary(m_full)$coefficients$cond
est <- co["conditionB","Estimate"]; se <- co["conditionB","Std. Error"]
irr    <- exp(est)
irr_ci <- exp(est + qnorm(c(0.025, 0.975))*se)
cat(sprintf("IRR (B vs A): %.3f  95%% CI [%.3f, %.3f]\n", irr, irr_ci[1], irr_ci[2]))


