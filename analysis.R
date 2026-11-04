###############################################################################
##  Customer Analytics at ShopCo -- individual project
##  Modern Computational Methods, Fall 2026
##  Antoine Foucart -- data file: project_data_Y0PV0T.csv
##
##  Run it section by section in RStudio: the plots appear in the Plots pane.
##  All numbers and figures of the report come from this script.
##
##  Packages: mice (imputation), survival (Kaplan-Meier, Cox),
##            cluster (silhouette), ca (correspondence analysis),
##            pROC (ROC curve). splines comes with R.
##
##  Sections
##    0. Setup
##    1. Task 1 - Missing values and multiple imputation
##    2. Task 2 - PCA, factor analysis, canonical correlation
##    3. Task 3 - Region x channel (CA) and clustering
##    4. Task 4 - What drives annual spend?
##    5. Task 5 - Who upgrades to premium? (logit / probit)
##    6. Task 6 - Retention program (survival analysis)
###############################################################################


## ---- 0. Setup ---------------------------------------------------------------

library(mice)
library(survival)
library(cluster)
library(ca)
library(pROC)
library(splines)


# prints a visible header in the console (same idea as in the course scripts)
section <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n",
                               strrep("=", 70), "\n", sep = "")

d <- read.csv("data/project_data_Y0PV0T.csv", stringsAsFactors = TRUE)

# education is an ordered variable: high school < bachelor < graduate
d$education_level <- factor(d$education_level, ordered = TRUE,
                            levels = c("high_school", "bachelor", "graduate"))

str(d)
cat("Number of customers:", nrow(d), "\n")

# the six satisfaction items (1-7 scale)
sat_items <- c("sat_price", "sat_value", "sat_promotions",
               "sat_staff", "sat_store", "sat_checkout")


###############################################################################
## ---- 1. TASK 1 - Missing values ------------------------------------------
###############################################################################

section("1a. How much is missing?")

n_miss   <- colSums(is.na(d))
pct_miss <- round(100 * n_miss / nrow(d), 1)
miss_table <- data.frame(n_missing = n_miss, pct_missing = pct_miss)
print(miss_table[miss_table$n_missing > 0, ])

cat("Complete rows:", sum(complete.cases(d)), "of", nrow(d), "\n")

# missing data pattern (1 = observed, 0 = missing)
print(md.pattern(d[, c("education_level", "sat_value", "spend_electronics")],
                 plot = FALSE))


section("1b. Is missingness related to the observed variables?")

# For each incomplete variable, compare customers with a missing value
# and customers with an observed value:
#   - numeric variables     -> t-test on the mean
#   - categorical variables -> chi-square test of independence
# If something is significant, MCAR is not plausible.

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
    new_row <- data.frame(incomplete_var = v, compared_with = x,
                          mean_if_missing = round(m_miss, 1),
                          mean_if_observed = round(m_obs, 1),
                          test = "t-test", p_value = round(p, 4))
    miss_tests <- rbind(miss_tests, new_row)
  }

  for (x in cat_vars) {
    p <- chisq.test(table(d[[x]], is_miss))$p.value
    new_row <- data.frame(incomplete_var = v, compared_with = x,
                          mean_if_missing = NA, mean_if_observed = NA,
                          test = "chi-square", p_value = round(p, 4))
    miss_tests <- rbind(miss_tests, new_row)
  }
}
print(miss_tests)

# Same question with one model per variable: logistic regression of the
# missingness indicator (1 = missing) on all the observed variables
for (v in incomplete) {
  d$miss_ind <- as.numeric(is.na(d[[v]]))
  m <- glm(miss_ind ~ age + gender + region + spend_grocery + annual_spend +
             sat_price + sat_staff + loyalty_program,
           data = d, family = binomial)
  cat("\nMissingness of", v, "- logistic regression\n")
  print(round(summary(m)$coefficients, 3))
}
d$miss_ind <- NULL

# Share of missing spend_electronics by age group and missing education by region
print(round(tapply(is.na(d$spend_electronics),
                   cut(d$age, c(17, 35, 50, 65, 75)), mean), 3))
print(round(tapply(is.na(d$education_level), d$region, mean), 3))

# Figure: age of customers with / without a value
par(mfrow = c(1, 3))
boxplot(d$age ~ is.na(d$spend_electronics), names = c("observed", "missing"),
        main = "spend_electronics", ylab = "Age", xlab = "", col = "grey85")
boxplot(d$age ~ is.na(d$sat_value), names = c("observed", "missing"),
        main = "sat_value", ylab = "Age", xlab = "", col = "grey85")
boxplot(d$age ~ is.na(d$education_level), names = c("observed", "missing"),
        main = "education_level", ylab = "Age", xlab = "", col = "grey85")
par(mfrow = c(1, 1))


section("1c. Multiple imputation with mice")

# One method per incomplete variable, chosen from its type:
#   education_level  : ordered factor   -> "polr" (proportional odds model)
#   sat_value        : rating 1-7       -> "pmm"  (only imputes values that exist)
#   spend_electronics: skewed, >= 0     -> "pmm"
meth <- make.method(d)
meth["education_level"]   <- "polr"
meth["sat_value"]         <- "pmm"
meth["spend_electronics"] <- "pmm"
print(meth[meth != ""])

# all variables are predictors, except the identifier
pred <- make.predictorMatrix(d)
pred[, "customer_id"] <- 0

imp <- mice(d, m = 20, maxit = 20, method = meth, predictorMatrix = pred,
            seed = 2026, printFlag = FALSE)
print(imp$loggedEvents)          # NULL = no problem

# Convergence check: the 20 chains should mix, with no trend
print(plot(imp, c("sat_value", "spend_electronics"), layout = c(2, 2)))

# Observed (blue) vs imputed (red) densities
print(densityplot(imp, ~ spend_electronics + sat_value))

# imp$imp$education_level has one row per missing customer and one column
# per imputed data set: we compare the shares of the imputed values with
# the observed shares
edu_imputed <- unlist(imp$imp$education_level)
cat("\nEducation level: observed shares\n")
print(round(prop.table(table(d$education_level)), 3))
cat("Education level: imputed shares (all 20 data sets)\n")
print(round(prop.table(table(edu_imputed)), 3))


section("1d. Mean of spend_electronics: complete cases vs multiple imputation")

# (a) complete cases
x_obs   <- d$spend_electronics[!is.na(d$spend_electronics)]
cc_mean <- mean(x_obs)
cc_se   <- sd(x_obs) / sqrt(length(x_obs))
cat("Complete cases: mean =", round(cc_mean, 2),
    " 95% CI:", round(cc_mean - 1.96 * cc_se, 2), "-",
    round(cc_mean + 1.96 * cc_se, 2), "\n")

# (b) multiple imputation, pooled with Rubin's rules -- first by hand ...
m <- 20
means <- numeric(m)
vars  <- numeric(m)
for (i in 1:m) {
  di <- complete(imp, i)
  means[i] <- mean(di$spend_electronics)
  vars[i]  <- var(di$spend_electronics) / nrow(di)    # squared standard error
}
Q_bar <- mean(means)              # pooled estimate
U_bar <- mean(vars)               # within-imputation variance
B     <- var(means)               # between-imputation variance
T_var <- U_bar + (1 + 1 / m) * B  # total variance
cat("Rubin's rules by hand: mean =", round(Q_bar, 2),
    " se =", round(sqrt(T_var), 3), "\n")

# ... then with mice: same estimate and standard error
mi_fit  <- with(imp, lm(spend_electronics ~ 1))
mi_pool <- summary(pool(mi_fit), conf.int = TRUE)
print(mi_pool)

cat("Mean of the imputed values only:",
    round(mean(as.matrix(imp$imp$spend_electronics)), 2), "\n")
cat("Correlations with spend_electronics (complete cases):  age",
    round(cor(d$age, d$spend_electronics, use = "complete.obs"), 3),
    "  grocery", round(cor(d$spend_grocery, d$spend_electronics, use = "complete.obs"), 3), "\n")


###############################################################################
## ---- 2. TASK 2 - What does "satisfaction" mean? -------------------------
###############################################################################

# Exploratory methods -> we use the FIRST completed data set
d1 <- complete(imp, 1)
X  <- d1[, sat_items]

section("2a. PCA on the six satisfaction items")

R <- cor(X)
print(round(R, 2))

# PCA by hand: eigen decomposition of the correlation matrix ...
eig <- eigen(R, symmetric = TRUE)
lambda <- eig$values

# ... and with prcomp (scaled data): the eigenvalues are the same
pca <- prcomp(X, scale. = TRUE)
cat("max |eigen(R) - prcomp sdev^2| =", signif(max(abs(lambda - pca$sdev^2)), 3), "\n")

var_tab <- data.frame(component  = 1:6,
                      eigenvalue = round(lambda, 3),
                      pct_var    = round(100 * lambda / sum(lambda), 1),
                      cum_pct    = round(100 * cumsum(lambda) / sum(lambda), 1))
print(var_tab)
cat("Kaiser rule (eigenvalue > 1): keep", sum(lambda > 1), "components\n")
print(round(pca$rotation[, 1:2], 3))      # loadings of PC1 and PC2

plot(lambda, type = "b", pch = 19, xlab = "Component", ylab = "Eigenvalue",
     main = "Scree plot - six satisfaction items")
abline(h = 1, lty = 2, col = "red")


section("2b. Factor analysis: 2 factors, varimax rotation")

fa <- factanal(X, factors = 2, rotation = "varimax", scores = "regression")
print(fa, cutoff = 0.3)

# Which column is the "price" factor? (the order of factors is arbitrary)
L <- fa$loadings[, 1:2]
if (abs(L["sat_price", 1]) > abs(L["sat_price", 2])) {
  price_col <- 1
} else {
  price_col <- 2
}
serv_col <- 3 - price_col
cat("Price & Value factor = column", price_col,
    "   Store Experience factor = column", serv_col, "\n")


section("2c. Factor scores")

# Regression (Thomson) scores: F = Z R^{-1} L, with Z the standardised items.
# We write the formula ourselves so we can apply the SAME rule (loadings and
# correlation matrix of data set 1) to each of the 20 imputed data sets
# in Tasks 4 and 5.
mean_X <- colMeans(X)
sd_X   <- apply(X, 2, sd)

make_scores <- function(data) {
  Z <- scale(data[, sat_items], center = mean_X, scale = sd_X)
  S <- Z %*% solve(R) %*% L
  data.frame(f_price = S[, price_col], f_service = S[, serv_col])
}

sc1 <- make_scores(d1)
# check: our formula gives exactly the factanal scores
cat("Correlation with factanal scores:",
    round(cor(sc1$f_price, fa$scores[, price_col]), 4), "\n")
d1$f_price   <- sc1$f_price
d1$f_service <- sc1$f_service


section("2d. Canonical correlation: satisfaction items vs spending")

Y  <- d1[, c("spend_grocery", "spend_electronics", "annual_spend")]
cc <- cancor(scale(X), scale(Y))
cat("Canonical correlations:", round(cc$cor, 3), "\n")

# Bartlett's sequential test: row k tests "correlations k, k+1, ... are all 0"
n <- nrow(X); p <- ncol(X); q <- ncol(Y)
r <- length(cc$cor)
for (k in 1:r) {
  wilks <- prod(1 - cc$cor[k:r]^2)
  chi2  <- -(n - 1 - (p + q + 1) / 2) * log(wilks)
  df    <- (p - k + 1) * (q - k + 1)
  cat("Test of correlations", k, "to", r, ": chi2 =", round(chi2, 2),
      " df =", df, " p =", signif(1 - pchisq(chi2, df), 3), "\n")
}

# canonical variates and structure correlations (read them like loadings)
U <- scale(X) %*% cc$xcoef[, 1]
V <- scale(Y) %*% cc$ycoef[, 1]
cat("Structure correlations, satisfaction side:\n"); print(round(cor(X, U), 3))
cat("Structure correlations, spending side:\n");     print(round(cor(Y, V), 3))

plot(U, V, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Satisfaction variate U1", ylab = "Spending variate V1",
     main = paste("First canonical pair, r =", round(cc$cor[1], 2)))


###############################################################################
## ---- 3. TASK 3 - Channels and customer segments --------------------------
###############################################################################

section("3a. Region x preferred channel")

# region and channel have no missing values: we use the original data
tab <- table(d$region, d$preferred_channel)
print(tab)
print(round(100 * prop.table(tab, 1), 1))       # row percentages

chi <- chisq.test(tab)
print(chi)
cramer_v <- sqrt(chi$statistic / (sum(tab) * (min(dim(tab)) - 1)))
cat("Cramer's V:", round(cramer_v, 3), "\n")
cat("Standardised residuals (|value| > 2 = cell that departs from independence):\n")
print(round(chi$stdres, 2))


section("3b. Correspondence analysis")

# CA by hand: SVD of the matrix of standardised residuals
P_mat <- tab / sum(tab)                  # relative frequencies
r_mass <- rowSums(P_mat)                 # row masses
c_mass <- colSums(P_mat)                 # column masses
S_res <- diag(1 / sqrt(r_mass)) %*% (P_mat - r_mass %*% t(c_mass)) %*% diag(1 / sqrt(c_mass))
inertias <- svd(S_res)$d^2
cat("Principal inertias by hand:", round(inertias, 6), "\n")
cat("Total inertia =", round(sum(inertias), 6), " = chi2 / n =",
    round(chi$statistic / sum(tab), 6), "\n")

# ... and with the ca package (same inertias, plus contributions)
ca_res <- ca(tab)
print(summary(ca_res))

plot(ca_res, main = "Correspondence analysis: region x channel")


section("3c. Clustering on grocery and electronics spend")

# completed data set 1, both variables standardised so they count equally
S <- scale(d1[, c("spend_grocery", "spend_electronics")])

# (a) Ward hierarchical clustering
hc <- hclust(dist(S), method = "ward.D2")
plot(hc, labels = FALSE, main = "Ward hierarchical clustering", xlab = "", sub = "")
rect.hclust(hc, k = 3, border = "red")

# (b) k-means for k = 1..8: elbow (within SS) and average silhouette
set.seed(2026)
wss <- numeric(8)
sil <- numeric(8)
for (k in 1:8) {
  km_k <- kmeans(S, centers = k, nstart = 25)
  wss[k] <- km_k$tot.withinss
  if (k > 1) {
    sil[k] <- mean(silhouette(km_k$cluster, dist(S))[, 3])
  }
}
print(data.frame(k = 1:8, within_SS = round(wss, 1), avg_silhouette = round(sil, 3)))

sil_hc <- numeric(8)
for (k in 2:8) {
  sil_hc[k] <- mean(silhouette(cutree(hc, k), dist(S))[, 3])
}
cat("Silhouette, Ward:", round(sil_hc[2:8], 3), "\n")

par(mfrow = c(1, 2))
plot(1:8, wss, type = "b", pch = 19, xlab = "Number of clusters k",
     ylab = "Total within sum of squares", main = "Elbow (k-means)")
plot(2:8, sil[2:8], type = "b", pch = 19, xlab = "Number of clusters k",
     ylab = "Average silhouette width", main = "Silhouette",
     ylim = range(c(sil[2:8], sil_hc[2:8])))
lines(2:8, sil_hc[2:8], type = "b", pch = 17, col = "red")
legend("topright", c("k-means", "Ward"), col = c("black", "red"), pch = c(19, 17))
par(mfrow = c(1, 1))

# Final choice: k = 3 (elbow + highest silhouette for both methods)
K <- 3
set.seed(2026)
km <- kmeans(S, centers = K, nstart = 25)
hc_cl <- cutree(hc, K)
print(table(kmeans = km$cluster, hierarchical = hc_cl))

# Adjusted Rand index (1 = identical partitions), same function as in the course
ari <- function(a, b) {
  tab <- table(a, b)
  n   <- sum(tab)
  s   <- sum(choose(tab, 2))
  sa  <- sum(choose(rowSums(tab), 2))
  sb  <- sum(choose(colSums(tab), 2))
  e   <- sa * sb / choose(n, 2)
  (s - e) / ((sa + sb) / 2 - e)
}
cat("Adjusted Rand index k-means vs Ward:", round(ari(km$cluster, hc_cl), 3), "\n")


section("3d. Profile of the segments")

d1$cluster <- km$cluster
profile <- aggregate(cbind(spend_grocery, spend_electronics, annual_spend, age,
                           f_price, f_service, premium_upgrade) ~ cluster,
                     data = d1, FUN = mean)
profile$size <- as.vector(table(d1$cluster))
print(round(profile, 2))

# Names read from the cluster centres printed above:
#   cluster 1: high grocery, low electronics  -> Grocery Loyalists
#   cluster 2: low grocery, high electronics  -> Tech Enthusiasts
#   cluster 3: about 50 EUR on each           -> Balanced Shoppers
seg_name <- c("Grocery Loyalists", "Tech Enthusiasts", "Balanced Shoppers")
d1$segment <- factor(seg_name[d1$cluster],
                     levels = c("Balanced Shoppers", "Grocery Loyalists", "Tech Enthusiasts"))
print(table(d1$segment, d1$preferred_channel))

cols <- c("grey50", "darkgreen", "steelblue")
plot(d1$spend_grocery, d1$spend_electronics, col = cols[d1$segment], pch = 19,
     xlab = "Monthly grocery spend (EUR)", ylab = "Monthly electronics spend (EUR)",
     main = "Customer segments (k-means, k = 3)")
legend("topright", levels(d1$segment), col = cols, pch = 19, bty = "n")


# Data sets for Tasks 4 and 5: for each of the 20 imputed data sets we add
# the factor scores (recomputed with the imputed sat_value) and the segment
# (from data set 1). Models are then fitted on each data set and pooled.
data_list <- list()
for (i in 1:20) {
  di <- complete(imp, i)
  sc <- make_scores(di)
  di$f_price   <- sc$f_price
  di$f_service <- sc$f_service
  di$segment   <- d1$segment
  di$education <- factor(di$education_level, ordered = FALSE,
                         levels = c("high_school", "bachelor", "graduate"))
  data_list[[i]] <- di
}


###############################################################################
## ---- 4. TASK 4 - What drives annual spend? -------------------------------
###############################################################################

section("4a. Linear effect of age")

# age and annual_spend are complete: ordinary lm on the original data
m_lin <- lm(annual_spend ~ age, data = d)
print(summary(m_lin))


section("4b. Non-linear effect of age")

m_quad <- lm(annual_spend ~ age + I(age^2), data = d)
m_ns   <- lm(annual_spend ~ ns(age, df = 3), data = d)

par(mfrow = c(1, 2))
plot(d$age, residuals(m_lin), pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Residual (EUR)", main = "Residuals of linear model")
abline(h = 0, col = "red")
lines(lowess(d$age, residuals(m_lin)), col = "blue", lwd = 2)
plot(d$age, d$annual_spend, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Annual spend (EUR)", main = "Fitted age effects")
age_grid <- data.frame(age = 18:75)
lines(age_grid$age, predict(m_lin, age_grid), col = "red", lwd = 2)
lines(age_grid$age, predict(m_quad, age_grid), col = "darkgreen", lwd = 2)
lines(age_grid$age, predict(m_ns, age_grid), col = "blue", lwd = 2, lty = 2)
legend("bottomleft", c("linear", "quadratic", "spline df=3"),
       col = c("red", "darkgreen", "blue"), lty = c(1, 1, 2), lwd = 2, bty = "n")
par(mfrow = c(1, 1))

# partial F-tests (nested models) and AIC
print(anova(m_lin, m_quad))
print(anova(m_lin, m_ns))
print(AIC(m_lin, m_quad, m_ns))
print(summary(m_quad)$coefficients)

b <- coef(m_quad)
cat("Peak age (quadratic, -b1 / (2 b2)):", round(-b[2] / (2 * b[3]), 1), "\n")
p_ns <- predict(m_ns, age_grid)
cat("Peak age (spline):", age_grid$age[which.max(p_ns)],
    " predicted spend:", round(max(p_ns), 1), "\n")
cat("Quadratic predictions at age 20 / 45 / 70:",
    round(predict(m_quad, data.frame(age = c(20, 45, 70))), 1), "\n")


section("4c. Full model, pooled over the 20 imputed data sets")

fits_full  <- list()
fits_noseg <- list()
fits_nosat <- list()
fits_age   <- list()
for (i in 1:20) {
  di <- data_list[[i]]
  fits_full[[i]]  <- lm(annual_spend ~ age + I(age^2) + f_price + f_service + segment, data = di)
  fits_noseg[[i]] <- lm(annual_spend ~ age + I(age^2) + f_price + f_service, data = di)
  fits_nosat[[i]] <- lm(annual_spend ~ age + I(age^2) + segment, data = di)
  fits_age[[i]]   <- lm(annual_spend ~ age + I(age^2), data = di)
}
# as.mira() turns our list of models into an object that pool() accepts
fits_full  <- as.mira(fits_full)
fits_noseg <- as.mira(fits_noseg)
fits_nosat <- as.mira(fits_nosat)
fits_age   <- as.mira(fits_age)

pool_full <- pool(fits_full)
print(summary(pool_full, conf.int = TRUE))
print(pool.r.squared(fits_full))
cat("R2 with age and age^2 only:\n")
print(pool.r.squared(fits_age))

# pooled Wald tests: do the segments matter? do the satisfaction scores matter?
print(D1(fits_full, fits_noseg))
print(D1(fits_full, fits_nosat))

est <- pool_full$pooled$estimate
cat("Peak age in the full model:", round(-est[2] / (2 * est[3]), 1), "\n")


###############################################################################
## ---- 5. TASK 5 - Who upgrades to premium? --------------------------------
###############################################################################

section("5a. Logit and probit, pooled over the 20 imputed data sets")

cat("Share upgraded:", round(mean(d$premium_upgrade), 3), "\n")

fits_logit  <- list()
fits_probit <- list()
for (i in 1:20) {
  di <- data_list[[i]]
  fits_logit[[i]]  <- glm(premium_upgrade ~ education + spend_grocery + f_price + f_service,
                          data = di, family = binomial(link = "logit"))
  fits_probit[[i]] <- glm(premium_upgrade ~ education + spend_grocery + f_price + f_service,
                          data = di, family = binomial(link = "probit"))
}
pool_logit  <- pool(as.mira(fits_logit))
pool_probit <- pool(as.mira(fits_probit))

cat("\nLogit - odds ratios:\n")
print(summary(pool_logit, conf.int = TRUE, exponentiate = TRUE))
b_logit <- pool_logit$pooled$estimate
cat("Odds ratio for +10 EUR grocery spend:", round(exp(10 * b_logit[4]), 3), "\n")

cat("\nProbit - coefficients:\n")
print(summary(pool_probit, conf.int = TRUE))

# Logit coefficients are about 1.6-1.8 times the probit ones
b_probit <- pool_probit$pooled$estimate
compare <- data.frame(term   = pool_logit$pooled$term,
                      logit  = round(b_logit, 3),
                      probit = round(b_probit, 3),
                      ratio  = round(b_logit / b_probit, 2))
print(compare)


section("5b. Predicted probabilities")

# average of the fitted probabilities over the 20 data sets
p_logit  <- rep(0, nrow(d))
p_probit <- rep(0, nrow(d))
for (i in 1:20) {
  p_logit  <- p_logit  + fitted(fits_logit[[i]])  / 20
  p_probit <- p_probit + fitted(fits_probit[[i]]) / 20
}
cat("Correlation of predicted probabilities:", round(cor(p_logit, p_probit), 4), "\n")
cat("Largest absolute difference:", round(max(abs(p_logit - p_probit)), 4), "\n")

plot(p_logit, p_probit, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "P(upgrade) - logit", ylab = "P(upgrade) - probit",
     main = "Logit vs probit predictions")
abline(0, 1, col = "red")


section("5c. ROC curve, AUC and threshold")

roc_logit  <- roc(d$premium_upgrade, p_logit, quiet = TRUE)
roc_probit <- roc(d$premium_upgrade, p_probit, quiet = TRUE)
cat("AUC logit:", round(auc(roc_logit), 3),
    " 95% CI:", round(ci.auc(roc_logit)[c(1, 3)], 3), "\n")
cat("AUC probit:", round(auc(roc_probit), 3), "\n")

# Youden threshold: maximises sensitivity + specificity
best <- coords(roc_logit, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
print(best)
thr <- best$threshold[1]

plot(roc_logit, main = paste("ROC curve - AUC =", round(auc(roc_logit), 3)),
     legacy.axes = TRUE)
points(best$specificity[1], best$sensitivity[1], pch = 19, col = "red")

# confusion matrices at 0.5 and at the Youden threshold
y <- d$premium_upgrade
for (t in c(0.5, thr)) {
  pred_class <- as.numeric(p_logit >= t)
  TP <- sum(pred_class == 1 & y == 1)
  TN <- sum(pred_class == 0 & y == 0)
  FP <- sum(pred_class == 1 & y == 0)
  FN <- sum(pred_class == 0 & y == 1)
  cat("\nConfusion matrix at threshold", round(t, 3), "\n")
  print(table(predicted = pred_class, actual = y))
  cat("Accuracy:",    round((TP + TN) / length(y), 3),
      " Sensitivity:", round(TP / (TP + FN), 3),
      " Specificity:", round(TN / (TN + FP), 3),
      " Precision:",   round(TP / (TP + FP), 3), "\n")
}


###############################################################################
## ---- 6. TASK 6 - Did the retention program work? ------------------------
###############################################################################

# all the variables used here are complete: no imputation needed

section("6a. Raw numbers")

print(table(loyalty = d$loyalty_program, churn = d$churn_event))
print(round(prop.table(table(d$loyalty_program, d$churn_event), 1), 3))
print(aggregate(months_observed ~ loyalty_program, data = d, FUN = mean))
print(table(loyalty = d$loyalty_program,
            censored_before_6 = d$churn_event == 0 & d$months_observed < 6))


section("6b. Kaplan-Meier curves and log-rank test")

km_fit <- survfit(Surv(months_observed, churn_event) ~ loyalty_program, data = d)
print(km_fit)
print(summary(km_fit, times = c(3, 6, 9, 11.8)))

plot(km_fit, col = c("red", "blue"), lwd = 2, conf.int = TRUE,
     xlab = "Months since start of observation", ylab = "Share still a customer",
     main = "Kaplan-Meier retention curves")
legend("bottomleft", c("Not enrolled", "Enrolled in program"),
       col = c("red", "blue"), lwd = 2, bty = "n")

lr <- survdiff(Surv(months_observed, churn_event) ~ loyalty_program, data = d)
print(lr)


section("6c. Cox model and proportional hazards")

cox <- coxph(Surv(months_observed, churn_event) ~ loyalty_program + age, data = d)
print(summary(cox))
hr <- exp(coef(cox))
cat("Enrolled customers leave at a rate", round(100 * (1 - hr[1]), 1), "% lower\n")
cat("Each extra 10 years of age changes the churn rate by",
    round(100 * (exp(10 * coef(cox)[2]) - 1), 1), "%\n")

zph <- cox.zph(cox)
print(zph)                       # p > 0.05: no evidence against proportional hazards

par(mfrow = c(1, 2))
plot(zph[1], main = "Schoenfeld residuals: loyalty_program")
abline(h = coef(cox)[1], col = "red", lty = 2)
plot(zph[2], main = "Schoenfeld residuals: age")
abline(h = coef(cox)[2], col = "red", lty = 2)
par(mfrow = c(1, 1))


section("6d. Discrete time: did the customer leave within 6 months?")

# Customers censored before month 6 are removed: we do not know whether
# they would have left before month 6 or not.
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
cat("Ratio logit / probit coefficients:", round(coef(lg6) / coef(pb6), 2), "\n")

# predicted probability of leaving within 6 months at the mean age
nd <- data.frame(loyalty_program = c(0, 1), age = mean(d6$age))
cat("P(leave within 6 months) at mean age, not enrolled / enrolled\n")
cat("  logit :", round(predict(lg6, nd, type = "response"), 3), "\n")
cat("  probit:", round(predict(pb6, nd, type = "response"), 3), "\n")
cat("  Kaplan-Meier 1 - S(6):", round(1 - summary(km_fit, times = 6)$surv, 3), "\n")


section("Session info")
print(sessionInfo())
