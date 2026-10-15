# =====================================================================
# Customer Analytics at ShopCo - Individual project
# Antoine Foucart - data file: project_data_Y0PV0T.csv
#
# Run from the repository root:  Rscript analysis.R
# Every number and figure of the report is produced by this script.
# Figures are written to the folder figures/
# =====================================================================

library(mice)       # multiple imputation
library(survival)   # Kaplan-Meier, Cox model
library(cluster)    # silhouette
library(ca)         # correspondence analysis
library(pROC)       # ROC curve and AUC
library(splines)    # natural splines (comes with R)

dir.create("figures", showWarnings = FALSE)

d <- read.csv("data/project_data_Y0PV0T.csv", stringsAsFactors = TRUE)
d$education_level <- factor(d$education_level, ordered = TRUE,
                            levels = c("high_school", "bachelor", "graduate"))

cat("Number of customers:", nrow(d), "\n")

sat_items <- c("sat_price", "sat_value", "sat_promotions",
               "sat_staff", "sat_store", "sat_checkout")


# =====================================================================
# TASK 1 - Data quality and missing values
# =====================================================================
cat("\n================ TASK 1 ================\n")

# --- 1a. How much is missing? ---------------------------------------
n_miss   <- colSums(is.na(d))
pct_miss <- round(100 * n_miss / nrow(d), 1)
miss_table <- data.frame(n_missing = n_miss, pct_missing = pct_miss)
print(miss_table[miss_table$n_missing > 0, ])

cat("Complete rows:", sum(complete.cases(d)), "of", nrow(d), "\n")
print(md.pattern(d[, c("education_level", "sat_value", "spend_electronics")],
                 plot = FALSE))

# --- 1b. Is missingness related to observed variables? --------------
# For each incomplete variable, compare customers with and without a value
incomplete <- c("education_level", "sat_value", "spend_electronics")
num_vars   <- c("age", "spend_grocery", "annual_spend", "sat_price", "sat_staff")
cat_vars   <- c("gender", "region", "preferred_channel", "loyalty_program")

miss_tests <- data.frame()
for (v in incomplete) {
  is_miss <- is.na(d[[v]])
  for (x in num_vars) {
    m_miss <- mean(d[[x]][is_miss])
    m_obs  <- mean(d[[x]][!is_miss])
    p <- t.test(d[[x]] ~ is_miss)$p.value
    miss_tests <- rbind(miss_tests,
                        data.frame(incomplete_var = v, compared_with = x,
                                   mean_if_missing = round(m_miss, 1),
                                   mean_if_observed = round(m_obs, 1),
                                   test = "t-test", p_value = round(p, 4)))
  }
  for (x in cat_vars) {
    p <- chisq.test(table(d[[x]], is_miss))$p.value
    miss_tests <- rbind(miss_tests,
                        data.frame(incomplete_var = v, compared_with = x,
                                   mean_if_missing = NA, mean_if_observed = NA,
                                   test = "chi-square", p_value = round(p, 4)))
  }
}
print(miss_tests)

# A logistic regression of the missingness indicator on all observed
# variables summarises the same question in one model per variable
for (v in incomplete) {
  d$miss_ind <- as.numeric(is.na(d[[v]]))
  m <- glm(miss_ind ~ age + gender + region + spend_grocery + annual_spend +
             sat_price + sat_staff + loyalty_program,
           data = d, family = binomial)
  cat("\nMissingness of", v, "- logistic regression\n")
  print(round(summary(m)$coefficients, 3))
}
d$miss_ind <- NULL

# Figure: age of customers with / without spend_electronics
png("figures/t1_missing_age.png", width = 1400, height = 650, res = 150)
par(mfrow = c(1, 3))
boxplot(d$age ~ is.na(d$spend_electronics), names = c("observed", "missing"),
        main = "spend_electronics", ylab = "Age", xlab = "", col = "grey85")
boxplot(d$age ~ is.na(d$sat_value), names = c("observed", "missing"),
        main = "sat_value", ylab = "Age", xlab = "", col = "grey85")
boxplot(d$age ~ is.na(d$education_level), names = c("observed", "missing"),
        main = "education_level", ylab = "Age", xlab = "", col = "grey85")
dev.off()

# --- 1c. Multiple imputation with mice ------------------------------
meth <- make.method(d)
meth["education_level"]   <- "polr"   # ordered factor -> prop. odds model
meth["sat_value"]         <- "pmm"    # 1-7 rating -> pmm keeps valid values
meth["spend_electronics"] <- "pmm"    # skewed, >= 0 -> pmm keeps it realistic
print(meth[meth != ""])

pred <- make.predictorMatrix(d)
pred[, "customer_id"] <- 0            # an identifier must not predict anything

imp <- mice(d, m = 20, maxit = 20, method = meth, predictorMatrix = pred,
            seed = 2026, printFlag = FALSE)
print(imp$loggedEvents)

# Convergence check: trace plots of mean and sd of imputed values
png("figures/t1_convergence.png", width = 1400, height = 750, res = 150)
print(plot(imp, c("sat_value", "spend_electronics"), layout = c(2, 2)))
dev.off()

# Diagnostic: imputed (red) vs observed (blue) distributions
png("figures/t1_density.png", width = 1400, height = 650, res = 150)
print(densityplot(imp, ~ spend_electronics + sat_value))
dev.off()

# Imputed education levels vs observed levels
cat("\nEducation level: observed shares vs imputed shares (all 20 sets)\n")
long_imp <- complete(imp, "long")
edu_imputed <- long_imp$education_level[rep(is.na(d$education_level), 20)]
print(round(prop.table(table(d$education_level)), 3))
print(round(prop.table(table(edu_imputed)), 3))

# --- 1d. Mean of spend_electronics: complete cases vs pooled MI -----
cc_mean <- mean(d$spend_electronics, na.rm = TRUE)
cc_se   <- sd(d$spend_electronics, na.rm = TRUE) / sqrt(sum(!is.na(d$spend_electronics)))
mi_fit  <- with(imp, lm(spend_electronics ~ 1))
mi_pool <- summary(pool(mi_fit), conf.int = TRUE)
cat("\nMean spend_electronics, complete cases:", round(cc_mean, 2),
    " (95% CI", round(cc_mean - 1.96 * cc_se, 2), "-", round(cc_mean + 1.96 * cc_se, 2), ")\n")
cat("Mean spend_electronics, multiple imputation (pooled):\n")
print(mi_pool)
cat("Mean of imputed values only:",
    round(mean(long_imp$spend_electronics[rep(is.na(d$spend_electronics), 20)]), 2), "\n")


# =====================================================================
# TASK 2 - What does "satisfaction" mean?
# =====================================================================
cat("\n================ TASK 2 ================\n")
# Exploratory methods: we use the FIRST completed data set
d1 <- complete(imp, 1)

X <- d1[, sat_items]
print(round(cor(X), 2))

# --- 2a. PCA ---------------------------------------------------------
pca <- prcomp(X, scale. = TRUE)
eig <- pca$sdev^2
var_tab <- data.frame(component = 1:6, eigenvalue = round(eig, 3),
                      pct_var = round(100 * eig / sum(eig), 1),
                      cum_pct = round(100 * cumsum(eig) / sum(eig), 1))
print(var_tab)
print(round(pca$rotation[, 1:2], 3))

png("figures/t2_scree.png", width = 1000, height = 700, res = 150)
plot(eig, type = "b", pch = 19, xlab = "Component", ylab = "Eigenvalue",
     main = "Scree plot - six satisfaction items")
abline(h = 1, lty = 2, col = "red")
dev.off()

# --- 2b. Factor analysis with 2 factors, varimax rotation ------------
fa <- factanal(X, factors = 2, rotation = "varimax", scores = "regression")
print(fa, cutoff = 0.3)

# Factor 1 / factor 2 order can change: identify which one is "price"
L <- fa$loadings[, 1:2]
price_col <- which.max(abs(L["sat_price", ]))
serv_col  <- 3 - price_col
cat("Price/value factor = column", price_col, "\n")

# --- 2c. Factor scores ----------------------------------------------
# Regression (Thomson) scores: F = Z %*% solve(R) %*% L
# Computed with the loadings of data set 1, but applied to EVERY imputed
# data set, so that the imputed sat_value enters the scores of Tasks 4-5
make_scores <- function(data) {
  Z <- scale(data[, sat_items], center = colMeans(X), scale = apply(X, 2, sd))
  S <- Z %*% solve(cor(X)) %*% L
  data.frame(f_price = S[, price_col], f_service = S[, serv_col])
}
sc1 <- make_scores(d1)
cat("Check: correlation with factanal scores:",
    round(cor(sc1$f_price, fa$scores[, price_col]), 4), "\n")
d1$f_price   <- sc1$f_price
d1$f_service <- sc1$f_service

# --- 2d. Canonical correlation --------------------------------------
Y <- d1[, c("spend_grocery", "spend_electronics", "annual_spend")]
cc <- cancor(scale(X), scale(Y))
cat("Canonical correlations:", round(cc$cor, 3), "\n")

# Bartlett's chi-square test that all canonical correlations are 0
n <- nrow(X); p <- ncol(X); q <- ncol(Y)
for (k in 1:length(cc$cor)) {
  lambda <- prod(1 - cc$cor[k:length(cc$cor)]^2)
  chi2 <- -(n - 1 - (p + q + 1) / 2) * log(lambda)
  df <- (p - k + 1) * (q - k + 1)
  cat("Test of correlations", k, "to", length(cc$cor), ": chi2 =", round(chi2, 2),
      " df =", df, " p =", signif(1 - pchisq(chi2, df), 3), "\n")
}

# Canonical variates and their correlations with the original variables
U <- scale(X) %*% cc$xcoef[, 1]
V <- scale(Y) %*% cc$ycoef[, 1]
cat("Structure correlations (first canonical variates)\n")
print(round(cor(X, U), 3))
print(round(cor(Y, V), 3))

png("figures/t2_cancor.png", width = 900, height = 800, res = 150)
plot(U, V, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Satisfaction variate U1", ylab = "Spending variate V1",
     main = paste("First canonical pair, r =", round(cc$cor[1], 2)))
dev.off()


# =====================================================================
# TASK 3 - Channels and customer segments
# =====================================================================
cat("\n================ TASK 3 ================\n")

# --- 3a. Region x channel (both complete, no imputation needed) -----
tab <- table(d$region, d$preferred_channel)
print(tab)
print(round(100 * prop.table(tab, 1), 1))
chi <- chisq.test(tab)
print(chi)
cat("Cramer's V:", round(sqrt(chi$statistic / (sum(tab) * (min(dim(tab)) - 1))), 3), "\n")
print(round(chi$stdres, 2))

ca_res <- ca(tab)
print(summary(ca_res))

png("figures/t3_ca.png", width = 1000, height = 900, res = 150)
plot(ca_res, main = "Correspondence analysis: region x channel")
dev.off()

# --- 3b. Clustering on grocery and electronics spend (data set 1) ---
S <- scale(d1[, c("spend_grocery", "spend_electronics")])

hc <- hclust(dist(S), method = "ward.D2")
png("figures/t3_dendrogram.png", width = 1200, height = 700, res = 150)
plot(hc, labels = FALSE, main = "Ward hierarchical clustering", xlab = "", sub = "")
rect.hclust(hc, k = 3, border = "red")
dev.off()

set.seed(2026)
wss <- numeric(8)
sil <- numeric(8)
for (k in 1:8) {
  km_k <- kmeans(S, centers = k, nstart = 25)
  wss[k] <- km_k$tot.withinss
  if (k > 1) sil[k] <- mean(silhouette(km_k$cluster, dist(S))[, 3])
}
print(data.frame(k = 1:8, within_SS = round(wss, 1), avg_silhouette = round(sil, 3)))

sil_hc <- numeric(8)
for (k in 2:8) sil_hc[k] <- mean(silhouette(cutree(hc, k), dist(S))[, 3])
cat("Silhouette, hierarchical:", round(sil_hc[2:8], 3), "\n")

png("figures/t3_elbow_silhouette.png", width = 1400, height = 650, res = 150)
par(mfrow = c(1, 2))
plot(1:8, wss, type = "b", pch = 19, xlab = "Number of clusters k",
     ylab = "Total within sum of squares", main = "Elbow (k-means)")
plot(2:8, sil[2:8], type = "b", pch = 19, xlab = "Number of clusters k",
     ylab = "Average silhouette width", main = "Silhouette",
     ylim = range(c(sil[2:8], sil_hc[2:8])))
lines(2:8, sil_hc[2:8], type = "b", pch = 17, col = "red")
legend("topright", c("k-means", "Ward"), col = c("black", "red"), pch = c(19, 17))
dev.off()

K <- 3
set.seed(2026)
km <- kmeans(S, centers = K, nstart = 25)
hc_cl <- cutree(hc, K)
print(table(kmeans = km$cluster, hierarchical = hc_cl))

# Adjusted Rand index (written by hand, no extra package)
ari <- function(a, b) {
  t <- table(a, b)
  comb <- function(x) sum(x * (x - 1) / 2)
  index <- comb(t); exp_i <- comb(rowSums(t)) * comb(colSums(t)) / comb(sum(t))
  max_i <- (comb(rowSums(t)) + comb(colSums(t))) / 2
  (index - exp_i) / (max_i - exp_i)
}
cat("Adjusted Rand index k-means vs Ward:", round(ari(km$cluster, hc_cl), 3), "\n")

# --- 3c. Profile the k-means segments -------------------------------
d1$cluster <- km$cluster
profile <- aggregate(cbind(spend_grocery, spend_electronics, annual_spend, age,
                           f_price, f_service, premium_upgrade) ~ cluster,
                     data = d1, FUN = mean)
profile$size <- as.vector(table(d1$cluster))
print(round(profile, 2))

# Name segments from their centres (done automatically so the names do
# not depend on the arbitrary cluster numbers)
centres <- aggregate(cbind(spend_grocery, spend_electronics) ~ cluster, d1, mean)
seg_name <- character(K)
seg_name[which.max(centres$spend_electronics)] <- "Tech Enthusiasts"
seg_name[which.max(centres$spend_grocery)]     <- "Grocery Loyalists"
seg_name[seg_name == ""]                       <- "Balanced Shoppers"
print(cbind(centres, seg_name))
d1$segment <- factor(seg_name[d1$cluster],
                     levels = c("Balanced Shoppers", "Grocery Loyalists", "Tech Enthusiasts"))

print(table(d1$segment, d1$preferred_channel))

png("figures/t3_segments.png", width = 1000, height = 850, res = 150)
cols <- c("grey50", "darkgreen", "steelblue")
plot(d1$spend_grocery, d1$spend_electronics, col = cols[d1$segment], pch = 19,
     xlab = "Monthly grocery spend (EUR)", ylab = "Monthly electronics spend (EUR)",
     main = "Customer segments (k-means, k = 3)")
legend("topright", levels(d1$segment), col = cols, pch = 19, bty = "n")
dev.off()

# Add factor scores (per imputation) and segment (from data set 1)
# to every imputed data set -> new mids object for Tasks 4 and 5
long <- complete(imp, "long", include = TRUE)
long$f_price   <- NA
long$f_service <- NA
for (i in 1:20) {
  rows <- long$.imp == i
  sc <- make_scores(long[rows, ])
  long$f_price[rows]   <- sc$f_price
  long$f_service[rows] <- sc$f_service
}
long$segment <- rep(d1$segment, 21)
long$education <- factor(long$education_level, ordered = FALSE,
                         levels = c("high_school", "bachelor", "graduate"))
imp2 <- as.mids(long)


# =====================================================================
# TASK 4 - What drives annual spend?
# =====================================================================
cat("\n================ TASK 4 ================\n")
# age and annual_spend have no missing values: plain lm on d

m_lin <- lm(annual_spend ~ age, data = d)
print(summary(m_lin))

png("figures/t4_resid_age.png", width = 1400, height = 650, res = 150)
par(mfrow = c(1, 2))
plot(d$age, residuals(m_lin), pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Residual (EUR)", main = "Residuals of linear model")
abline(h = 0, col = "red")
lines(lowess(d$age, residuals(m_lin)), col = "blue", lwd = 2)
plot(d$age, d$annual_spend, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Annual spend (EUR)", main = "Fitted age effects")
age_grid <- data.frame(age = 18:75)
m_quad <- lm(annual_spend ~ age + I(age^2), data = d)
m_ns   <- lm(annual_spend ~ ns(age, df = 3), data = d)
lines(age_grid$age, predict(m_lin, age_grid), col = "red", lwd = 2)
lines(age_grid$age, predict(m_quad, age_grid), col = "darkgreen", lwd = 2)
lines(age_grid$age, predict(m_ns, age_grid), col = "blue", lwd = 2, lty = 2)
legend("bottomleft", c("linear", "quadratic", "spline df=3"),
       col = c("red", "darkgreen", "blue"), lty = c(1, 1, 2), lwd = 2, bty = "n")
dev.off()

print(anova(m_lin, m_quad))
print(anova(m_lin, m_ns))
print(AIC(m_lin, m_quad, m_ns))
print(summary(m_quad)$coefficients)

b <- coef(m_quad)
cat("Peak age (quadratic, -b1/2b2):", round(-b[2] / (2 * b[3]), 1), "\n")
p_ns <- predict(m_ns, age_grid)
cat("Peak age (spline):", age_grid$age[which.max(p_ns)],
    " predicted spend:", round(max(p_ns), 1), "\n")
p_q <- predict(m_quad, age_grid)
cat("Quadratic predictions at 20/45/70:",
    round(predict(m_quad, data.frame(age = c(20, 45, 70))), 1), "\n")

# Full model, fitted on each imputed data set and pooled (Rubin)
fit_full <- with(imp2, lm(annual_spend ~ age + I(age^2) + f_price + f_service + segment))
pool_full <- pool(fit_full)
print(summary(pool_full, conf.int = TRUE))
print(pool.r.squared(fit_full))

# Does the block of segments matter? (pooled Wald test, D1)
fit_noseg <- with(imp2, lm(annual_spend ~ age + I(age^2) + f_price + f_service))
print(D1(fit_full, fit_noseg))
fit_nosat <- with(imp2, lm(annual_spend ~ age + I(age^2) + segment))
print(D1(fit_full, fit_nosat))

cat("Peak age in full model:",
    round(-pool_full$pooled$estimate[2] / (2 * pool_full$pooled$estimate[3]), 1), "\n")


# =====================================================================
# TASK 5 - Who upgrades to premium?
# =====================================================================
cat("\n================ TASK 5 ================\n")
cat("Share upgraded:", round(mean(d$premium_upgrade), 3), "\n")

fit_logit <- with(imp2, glm(premium_upgrade ~ education + spend_grocery + f_price + f_service,
                            family = binomial(link = "logit")))
fit_probit <- with(imp2, glm(premium_upgrade ~ education + spend_grocery + f_price + f_service,
                             family = binomial(link = "probit")))
pool_logit  <- pool(fit_logit)
pool_probit <- pool(fit_probit)

res_logit <- summary(pool_logit, conf.int = TRUE, exponentiate = TRUE)
print(res_logit)
cat("Odds ratio for +10 EUR grocery spend:",
    round(exp(10 * pool_logit$pooled$estimate[4]), 3), "\n")

res_probit <- summary(pool_probit, conf.int = TRUE)
print(res_probit)

compare <- data.frame(term = pool_logit$pooled$term,
                      logit = round(pool_logit$pooled$estimate, 3),
                      probit = round(pool_probit$pooled$estimate, 3))
compare$ratio <- round(compare$logit / compare$probit, 2)
print(compare)

# Predicted probabilities: average over the 20 imputed data sets
p_logit  <- rep(0, nrow(d))
p_probit <- rep(0, nrow(d))
for (i in 1:20) {
  p_logit  <- p_logit  + fitted(fit_logit$analyses[[i]])  / 20
  p_probit <- p_probit + fitted(fit_probit$analyses[[i]]) / 20
}
cat("Correlation of predicted probabilities:", round(cor(p_logit, p_probit), 4), "\n")
cat("Max absolute difference:", round(max(abs(p_logit - p_probit)), 4), "\n")

png("figures/t5_logit_probit.png", width = 800, height = 800, res = 150)
plot(p_logit, p_probit, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "P(upgrade) - logit", ylab = "P(upgrade) - probit",
     main = "Logit vs probit predictions")
abline(0, 1, col = "red")
dev.off()

# ROC curve and AUC (in-sample)
roc_logit  <- roc(d$premium_upgrade, p_logit, quiet = TRUE)
roc_probit <- roc(d$premium_upgrade, p_probit, quiet = TRUE)
cat("AUC logit:", round(auc(roc_logit), 3), " 95% CI:",
    round(ci.auc(roc_logit)[c(1, 3)], 3), "\n")
cat("AUC probit:", round(auc(roc_probit), 3), "\n")

# Threshold: Youden index (maximises sensitivity + specificity)
best <- coords(roc_logit, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
print(best)
thr <- best$threshold[1]

png("figures/t5_roc.png", width = 850, height = 850, res = 150)
plot(roc_logit, main = paste("ROC curve - AUC =", round(auc(roc_logit), 3)),
     legacy.axes = TRUE)
points(best$specificity[1], best$sensitivity[1], pch = 19, col = "red")
dev.off()

for (t in c(0.5, thr)) {
  pred_class <- as.numeric(p_logit >= t)
  cm <- table(predicted = pred_class, actual = d$premium_upgrade)
  cat("\nConfusion matrix at threshold", round(t, 3), "\n")
  print(cm)
  cat("Accuracy:", round(mean(pred_class == d$premium_upgrade), 3),
      " Sensitivity:", round(sum(pred_class == 1 & d$premium_upgrade == 1) / sum(d$premium_upgrade == 1), 3),
      " Specificity:", round(sum(pred_class == 0 & d$premium_upgrade == 0) / sum(d$premium_upgrade == 0), 3),
      " Precision:", round(sum(pred_class == 1 & d$premium_upgrade == 1) / sum(pred_class == 1), 3), "\n")
}


# =====================================================================
# TASK 6 - Did the retention program work?
# =====================================================================
cat("\n================ TASK 6 ================\n")
# All variables used here are complete: no imputation needed

print(table(loyalty = d$loyalty_program, churn = d$churn_event))
print(round(prop.table(table(d$loyalty_program, d$churn_event), 1), 3))
print(aggregate(months_observed ~ loyalty_program, d, mean))
print(table(loyalty = d$loyalty_program,
            censored_before_6 = d$churn_event == 0 & d$months_observed < 6))

# --- 6a. Kaplan-Meier + log-rank ------------------------------------
km_fit <- survfit(Surv(months_observed, churn_event) ~ loyalty_program, data = d)
print(km_fit)
print(summary(km_fit, times = c(3, 6, 9, 11.8)))

png("figures/t6_km.png", width = 1100, height = 800, res = 150)
plot(km_fit, col = c("red", "blue"), lwd = 2, conf.int = TRUE,
     xlab = "Months since start of observation", ylab = "Share still a customer",
     main = "Kaplan-Meier retention curves")
legend("bottomleft", c("Not enrolled", "Enrolled in program"),
       col = c("red", "blue"), lwd = 2, bty = "n")
dev.off()

lr <- survdiff(Surv(months_observed, churn_event) ~ loyalty_program, data = d)
print(lr)

# --- 6b. Cox model --------------------------------------------------
cox <- coxph(Surv(months_observed, churn_event) ~ loyalty_program + age, data = d)
print(summary(cox))
hr <- exp(coef(cox))
cat("Enrolled customers leave at a rate", round(100 * (1 - hr[1]), 1), "% lower\n")
cat("Each extra 10 years of age changes the churn rate by",
    round(100 * (exp(10 * coef(cox)[2]) - 1), 1), "%\n")

zph <- cox.zph(cox)
print(zph)
png("figures/t6_zph.png", width = 1400, height = 650, res = 150)
par(mfrow = c(1, 2))
plot(zph[1], main = "Schoenfeld residuals: loyalty_program")
abline(h = coef(cox)[1], col = "red", lty = 2)
plot(zph[2], main = "Schoenfeld residuals: age")
abline(h = coef(cox)[2], col = "red", lty = 2)
dev.off()

# --- 6c. Discrete time: left within 6 months -------------------------
# Customers censored before month 6 are excluded: we do not know whether
# they would have left before month 6 or not
d6 <- d[!(d$churn_event == 0 & d$months_observed < 6), ]
d6$left6 <- as.numeric(d6$churn_event == 1 & d6$months_observed <= 6)
cat("Customers kept for the 6-month analysis:", nrow(d6), "of", nrow(d), "\n")
print(round(prop.table(table(d6$loyalty_program, d6$left6), 1), 3))

lg6 <- glm(left6 ~ loyalty_program + age, data = d6, family = binomial(link = "logit"))
pb6 <- glm(left6 ~ loyalty_program + age, data = d6, family = binomial(link = "probit"))
print(summary(lg6)$coefficients)
print(summary(pb6)$coefficients)
cat("Odds ratios (logit):\n")
print(round(exp(cbind(OR = coef(lg6), confint.default(lg6))), 3))
cat("Ratio logit/probit coefficients:", round(coef(lg6) / coef(pb6), 2), "\n")

# Predicted 6-month churn probability at the mean age
nd <- data.frame(loyalty_program = c(0, 1), age = mean(d6$age))
cat("Predicted P(leave within 6 months) at mean age, not enrolled / enrolled\n")
cat("  logit :", round(predict(lg6, nd, type = "response"), 3), "\n")
cat("  probit:", round(predict(pb6, nd, type = "response"), 3), "\n")
cat("  Kaplan-Meier 1 - S(6):",
    round(1 - summary(km_fit, times = 6)$surv, 3), "\n")

cat("\nSession info\n")
print(sessionInfo())
