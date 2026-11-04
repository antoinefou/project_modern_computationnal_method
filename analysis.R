###############################################################################
##  Customer Analytics at ShopCo -- individual project
##  Modern Computational Methods, Fall 2026
##  Antoine Foucart -- data file: project_data_Y0PV0T.csv
##
##  Run it section by section in RStudio: the plots appear in the Plots pane.
##  All numbers and figures of the report come from this script.
##  The comments give the values obtained with this data file.
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

data <- read.csv("data/project_data_Y0PV0T.csv", stringsAsFactors = TRUE)

# education is an ordered variable: high school < bachelor < graduate
data$education_level <- factor(data$education_level, ordered = TRUE,
                               levels = c("high_school", "bachelor", "graduate"))

str(data)
nrow(data)            # 560 customers

# the six satisfaction items (1-7 scale)
sat_items <- c("sat_price", "sat_value", "sat_promotions",
               "sat_staff", "sat_store", "sat_checkout")


###############################################################################
## ---- 1. TASK 1 - Missing values ------------------------------------------
###############################################################################

## ---- 1a. How much is missing? ----

n_miss   <- colSums(is.na(data))
pct_miss <- round(100 * n_miss / nrow(data), 1)
miss_table <- data.frame(n_missing = n_miss, pct_missing = pct_miss)
miss_table[miss_table$n_missing > 0, ]
# education_level 66 (11.8 %), sat_value 57 (10.2 %), spend_electronics 70 (12.5 %)
# all the other variables are complete

sum(complete.cases(data))     # 384 complete rows out of 560 (69 %)

# missing data pattern (1 = observed, 0 = missing)
md.pattern(data[, c("education_level", "sat_value", "spend_electronics")],
           plot = FALSE)


## ---- 1b. Is missingness related to the observed variables? ----

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
  is_miss <- is.na(data[[v]])

  for (x in num_vars) {
    m_miss <- mean(data[[x]][is_miss])
    m_obs  <- mean(data[[x]][!is_miss])
    p <- t.test(data[[x]] ~ is_miss)$p.value
    new_row <- data.frame(incomplete_var = v, compared_with = x,
                          mean_if_missing = round(m_miss, 1),
                          mean_if_observed = round(m_obs, 1),
                          test = "t-test", p_value = round(p, 4))
    miss_tests <- rbind(miss_tests, new_row)
  }

  for (x in cat_vars) {
    p <- chisq.test(table(data[[x]], is_miss))$p.value
    new_row <- data.frame(incomplete_var = v, compared_with = x,
                          mean_if_missing = NA, mean_if_observed = NA,
                          test = "chi-square", p_value = round(p, 4))
    miss_tests <- rbind(miss_tests, new_row)
  }
}
miss_tests
# spend_electronics: missing much more often for OLDER customers
#                    (mean age 53.6 vs 44.4, p < 0.001)
# education_level  : missing much more often in the WEST (p < 0.001)
# sat_value        : no test significant -> MCAR plausible

# Same question with one model per variable: logistic regression of the
# missingness indicator (1 = missing) on all the observed variables
for (v in incomplete) {
  data$miss_ind <- as.numeric(is.na(data[[v]]))
  m <- glm(miss_ind ~ age + gender + region + spend_grocery + annual_spend +
             sat_price + sat_staff + loyalty_program,
           data = data, family = binomial)
  print(v)
  print(round(summary(m)$coefficients, 3))
}
data$miss_ind <- NULL
# education_level  : regionWest significant (p < 0.001), age borderline (p = 0.044)
# sat_value        : nothing significant
# spend_electronics: age significant (p < 0.001)

# share of missing spend_electronics by age group
round(tapply(is.na(data$spend_electronics), cut(data$age, c(17, 35, 50, 65, 75)), mean), 3)
# 0.057 (18-35)  0.111  0.170  0.222 (65-75)

# share of missing education by region
round(tapply(is.na(data$education_level), data$region, mean), 3)
# East 0.079  North 0.063  South 0.087  West 0.281

# Figure: age of customers with / without a value
par(mfrow = c(1, 3))
boxplot(data$age ~ is.na(data$spend_electronics), names = c("observed", "missing"),
        main = "spend_electronics", ylab = "Age", xlab = "", col = "grey85")
boxplot(data$age ~ is.na(data$sat_value), names = c("observed", "missing"),
        main = "sat_value", ylab = "Age", xlab = "", col = "grey85")
boxplot(data$age ~ is.na(data$education_level), names = c("observed", "missing"),
        main = "education_level", ylab = "Age", xlab = "", col = "grey85")
par(mfrow = c(1, 1))


## ---- 1c. Multiple imputation with mice ----

# One method per incomplete variable, chosen from its type:
#   education_level  : ordered factor   -> "polr" (proportional odds model)
#   sat_value        : rating 1-7       -> "pmm"  (only imputes values that exist)
#   spend_electronics: skewed, >= 0     -> "pmm"
meth <- make.method(data)
meth["education_level"]   <- "polr"
meth["sat_value"]         <- "pmm"
meth["spend_electronics"] <- "pmm"
meth[meth != ""]

# all variables are predictors, except the identifier
pred <- make.predictorMatrix(data)
pred[, "customer_id"] <- 0

imp <- mice(data, m = 20, maxit = 20, method = meth, predictorMatrix = pred,
            seed = 2026, printFlag = FALSE)
imp$loggedEvents          # NULL = no problem

# Convergence check: the 20 chains should mix, with no trend
print(plot(imp, c("sat_value", "spend_electronics"), layout = c(2, 2)))

# Observed (blue) vs imputed (red) densities
print(densityplot(imp, ~ spend_electronics + sat_value))

# imp$imp$education_level has one row per missing customer and one column
# per imputed data set: we compare the shares of the imputed values with
# the observed shares
edu_imputed <- unlist(imp$imp$education_level)
round(prop.table(table(data$education_level)), 3)   # observed: 0.298 0.431 0.271
round(prop.table(table(edu_imputed)), 3)            # imputed : 0.289 0.442 0.269


## ---- 1d. Mean of spend_electronics: complete cases vs multiple imputation ----

# (a) complete cases
x_obs   <- data$spend_electronics[!is.na(data$spend_electronics)]
cc_mean <- mean(x_obs)
cc_se   <- sd(x_obs) / sqrt(length(x_obs))
cc_mean                                             # 55.34
c(cc_mean - 1.96 * cc_se, cc_mean + 1.96 * cc_se)   # 95% CI: 53.18 - 57.51

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
Q_bar                             # 55.90
sqrt(T_var)                       # se = 1.037

# ... then with mice: same estimate and standard error
mi_fit  <- with(imp, lm(spend_electronics ~ 1))
summary(pool(mi_fit), conf.int = TRUE)
# estimate 55.90, 95% CI 53.87 - 57.94

# mean of the imputed values only
mean(as.matrix(imp$imp$spend_electronics))          # 59.81

# why are the imputed values higher? electronics spend is not related to age
# but it is the mirror image of grocery spend, and the customers with a
# missing value spend a bit less on groceries
cor(data$age, data$spend_electronics, use = "complete.obs")            # 0.003
cor(data$spend_grocery, data$spend_electronics, use = "complete.obs")  # -0.893


###############################################################################
## ---- 2. TASK 2 - What does "satisfaction" mean? -------------------------
###############################################################################

# Exploratory methods -> we use the FIRST completed data set
d1 <- complete(imp, 1)
X  <- d1[, sat_items]

## ---- 2a. PCA on the six satisfaction items ----

R <- cor(X)
round(R, 2)
# two blocks: price / value / promotions (0.52-0.59)
#             staff / store / checkout   (0.54-0.56)
# correlations between the blocks are below 0.16

# PCA by hand: eigen decomposition of the correlation matrix ...
eig <- eigen(R, symmetric = TRUE)
lambda <- eig$values

# ... and with prcomp (scaled data): the eigenvalues are the same
pca <- prcomp(X, scale. = TRUE)
max(abs(lambda - pca$sdev^2))     # about 1e-15

var_tab <- data.frame(component  = 1:6,
                      eigenvalue = round(lambda, 3),
                      pct_var    = round(100 * lambda / sum(lambda), 1),
                      cum_pct    = round(100 * cumsum(lambda) / sum(lambda), 1))
var_tab
# eigenvalues 2.33 and 1.92, then all about 0.4-0.5
# the first two components explain 70.8 % of the variance

sum(lambda > 1)                   # Kaiser rule: keep 2 components
round(pca$rotation[, 1:2], 3)
# PC1: all loadings of the same sign (general satisfaction level)
# PC2: price items (+) against service items (-)

plot(lambda, type = "b", pch = 19, xlab = "Component", ylab = "Eigenvalue",
     main = "Scree plot - six satisfaction items")
abline(h = 1, lty = 2, col = "red")


## ---- 2b. Factor analysis: 2 factors, varimax rotation ----

fa <- factanal(X, factors = 2, rotation = "varimax", scores = "regression")
print(fa, cutoff = 0.3)
# test "2 factors are sufficient": chi2 = 2.14, df = 4, p = 0.71 -> not rejected
# Factor 1: price 0.80, value 0.71, promotions 0.74  -> "Price & Value"
# Factor 2: staff 0.76, store 0.73, checkout 0.74    -> "Store Experience"

# Which column is the "price" factor? (the order of factors is arbitrary)
L <- fa$loadings[, 1:2]
if (abs(L["sat_price", 1]) > abs(L["sat_price", 2])) {
  price_col <- 1
} else {
  price_col <- 2
}
serv_col <- 3 - price_col
price_col        # 1 = Price & Value
serv_col         # 2 = Store Experience


## ---- 2c. Factor scores ----

# Regression (Thomson) scores: F = Z R^{-1} L, with Z the standardised items.
# We write the formula ourselves so we can apply the SAME rule (loadings and
# correlation matrix of data set 1) to each of the 20 imputed data sets
# in Tasks 4 and 5.
mean_X <- colMeans(X)
sd_X   <- apply(X, 2, sd)

make_scores <- function(dat) {
  Z <- scale(dat[, sat_items], center = mean_X, scale = sd_X)
  S <- Z %*% solve(R) %*% L
  data.frame(f_price = S[, price_col], f_service = S[, serv_col])
}

sc1 <- make_scores(d1)
# check: our formula gives exactly the factanal scores
cor(sc1$f_price, fa$scores[, price_col])     # 1
d1$f_price   <- sc1$f_price
d1$f_service <- sc1$f_service


## ---- 2d. Canonical correlation: satisfaction items vs spending ----

Y  <- d1[, c("spend_grocery", "spend_electronics", "annual_spend")]
cc <- cancor(scale(X), scale(Y))
round(cc$cor, 3)        # 0.503 0.166 0.039

# Bartlett's sequential test: row k tests "correlations k, k+1, ... are all 0"
n <- nrow(X); p <- ncol(X); q <- ncol(Y)
r <- length(cc$cor)
bartlett <- data.frame()
for (k in 1:r) {
  wilks <- prod(1 - cc$cor[k:r]^2)
  chi2  <- -(n - 1 - (p + q + 1) / 2) * log(wilks)
  df    <- (p - k + 1) * (q - k + 1)
  bartlett <- rbind(bartlett, data.frame(from = k, chi2 = round(chi2, 2), df = df,
                                         p_value = signif(1 - pchisq(chi2, df), 3)))
}
bartlett
# all 3 = 0 : chi2 = 177.95, df = 18, p < 0.001 -> rejected
# 2 and 3   : chi2 = 16.33,  df = 10, p = 0.09  -> not rejected
# -> only the first canonical correlation (0.50) is significant

# canonical variates and structure correlations (read them like loadings)
U <- scale(X) %*% cc$xcoef[, 1]
V <- scale(Y) %*% cc$ycoef[, 1]
round(cor(X, U), 3)     # price 0.91, value 0.74, promotions 0.82, service items ~ 0
round(cor(Y, V), 3)     # annual_spend 0.98, grocery 0.14, electronics -0.14
# -> satisfaction with prices goes with a higher ANNUAL spend

plot(U, V, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Satisfaction variate U1", ylab = "Spending variate V1",
     main = paste("First canonical pair, r =", round(cc$cor[1], 2)))


###############################################################################
## ---- 3. TASK 3 - Channels and customer segments --------------------------
###############################################################################

## ---- 3a. Region x preferred channel ----

# region and channel have no missing values: we use the original data
tab <- table(data$region, data$preferred_channel)
tab
round(100 * prop.table(tab, 1), 1)       # row percentages
# West: 36.0 % Catalog (about 20 % elsewhere)
# East: 32.4 % Mobile App, North: 30.4 % Online

chi <- chisq.test(tab)
chi                                      # chi2 = 18.35, df = 9, p = 0.031

cramer_v <- sqrt(chi$statistic / (sum(tab) * (min(dim(tab)) - 1)))
cramer_v                                 # 0.105 -> weak association

# standardised residuals (|value| > 2 = cell that departs from independence)
round(chi$stdres, 2)
# West-Catalog +3.26, East-Mobile App +2.19, North-Online +2.04


## ---- 3b. Correspondence analysis ----

# CA by hand: SVD of the matrix of standardised residuals
P_mat <- tab / sum(tab)                  # relative frequencies
r_mass <- rowSums(P_mat)                 # row masses
c_mass <- colSums(P_mat)                 # column masses
S_res <- diag(1 / sqrt(r_mass)) %*% (P_mat - r_mass %*% t(c_mass)) %*% diag(1 / sqrt(c_mass))
inertias <- svd(S_res)$d^2
round(inertias, 6)                       # 0.021036 0.011733 (third ~ 0)
sum(inertias)                            # total inertia 0.0328 ...
chi$statistic / sum(tab)                 # ... = chi2 / n

# ... and with the ca package (same inertias, plus contributions)
ca_res <- ca(tab)
summary(ca_res)
# dimension 1 (64.2 %): West and Catalog against North and Online
# dimension 2 (35.8 %): East and Mobile App

plot(ca_res, main = "Correspondence analysis: region x channel")


## ---- 3c. Clustering on grocery and electronics spend ----

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
data.frame(k = 1:8, within_SS = round(wss, 1), avg_silhouette = round(sil, 3))
# within SS: 1118, 355, 120, 102, ... -> elbow at k = 3
# silhouette is highest at k = 3 (0.654)

sil_hc <- numeric(8)
for (k in 2:8) {
  sil_hc[k] <- mean(silhouette(cutree(hc, k), dist(S))[, 3])
}
round(sil_hc[2:8], 3)     # Ward: also highest at k = 3 (0.652)

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
table(kmeans = km$cluster, hierarchical = hc_cl)   # only 2 customers differ

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
ari(km$cluster, hc_cl)    # 0.989 -> the two methods agree


## ---- 3d. Profile of the segments ----

d1$cluster <- km$cluster
profile <- aggregate(cbind(spend_grocery, spend_electronics, annual_spend, age,
                           f_price, f_service, premium_upgrade) ~ cluster,
                     data = d1, FUN = mean)
profile$size <- as.vector(table(d1$cluster))
round(profile, 2)

# Names read from the cluster centres printed above:
#   cluster 1: grocery 80, electronics 20 (115 customers) -> Grocery Loyalists
#   cluster 2: grocery 19, electronics 80 (219 customers) -> Tech Enthusiasts
#   cluster 3: about 50 EUR on each       (226 customers) -> Balanced Shoppers
seg_name <- c("Grocery Loyalists", "Tech Enthusiasts", "Balanced Shoppers")
d1$segment <- factor(seg_name[d1$cluster],
                     levels = c("Balanced Shoppers", "Grocery Loyalists", "Tech Enthusiasts"))
table(d1$segment, d1$preferred_channel)

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

## ---- 4a. Linear effect of age ----

# age and annual_spend are complete: ordinary lm on the original data
m_lin <- lm(annual_spend ~ age, data = data)
summary(m_lin)
# annual_spend = 579 + 2.95 * age, R2 = 0.13


## ---- 4b. Non-linear effect of age ----

m_quad <- lm(annual_spend ~ age + I(age^2), data = data)
m_ns   <- lm(annual_spend ~ ns(age, df = 3), data = data)

par(mfrow = c(1, 2))
plot(data$age, residuals(m_lin), pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Residual (EUR)", main = "Residuals of linear model")
abline(h = 0, col = "red")
lines(lowess(data$age, residuals(m_lin)), col = "blue", lwd = 2)
plot(data$age, data$annual_spend, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Age", ylab = "Annual spend (EUR)", main = "Fitted age effects")
age_grid <- data.frame(age = 18:75)
lines(age_grid$age, predict(m_lin, age_grid), col = "red", lwd = 2)
lines(age_grid$age, predict(m_quad, age_grid), col = "darkgreen", lwd = 2)
lines(age_grid$age, predict(m_ns, age_grid), col = "blue", lwd = 2, lty = 2)
legend("bottomleft", c("linear", "quadratic", "spline df=3"),
       col = c("red", "darkgreen", "blue"), lty = c(1, 1, 2), lwd = 2, bty = "n")
par(mfrow = c(1, 1))
# the residuals have an inverted-U shape: the linear model is not adequate

# partial F-tests (nested models) and AIC
anova(m_lin, m_quad)        # F = 51.8, p < 0.001
anova(m_lin, m_ns)          # F = 25.4, p < 0.001
AIC(m_lin, m_quad, m_ns)    # 7013 / 6966 / 6968 -> keep the quadratic model
summary(m_quad)$coefficients

b <- coef(m_quad)
-b[2] / (2 * b[3])          # peak age 54.7 (quadratic)

p_ns <- predict(m_ns, age_grid)
age_grid$age[which.max(p_ns)]    # peak age 55 (spline) ...
max(p_ns)                        # ... with a predicted spend of 774 EUR

predict(m_quad, data.frame(age = c(20, 45, 70)))   # 575.8  754.8  732.6


## ---- 4c. Full model, pooled over the 20 imputed data sets ----

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
summary(pool_full, conf.int = TRUE)
# f_price                 +68.1 EUR per SD  (p < 0.001)
# f_service                -1.4 EUR         (p = 0.75, no effect)
# Grocery Loyalists       +56.8 EUR vs Balanced Shoppers
# Tech Enthusiasts        -56.5 EUR vs Balanced Shoppers

pool.r.squared(fits_full)     # R2 = 0.52
pool.r.squared(fits_age)      # R2 = 0.20 with age and age^2 only

# pooled Wald tests: do the segments matter? do the satisfaction scores matter?
D1(fits_full, fits_noseg)     # segments    : F = 56,  p < 0.001
D1(fits_full, fits_nosat)     # satisfaction: F = 113, p < 0.001

est <- pool_full$pooled$estimate
-est[2] / (2 * est[3])        # peak age still 54.7


###############################################################################
## ---- 5. TASK 5 - Who upgrades to premium? --------------------------------
###############################################################################

## ---- 5a. Logit and probit, pooled over the 20 imputed data sets ----

mean(data$premium_upgrade)    # 35.7 % upgraded

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

# logit - odds ratios
summary(pool_logit, conf.int = TRUE, exponentiate = TRUE)
# bachelor 1.86, graduate 2.81 (vs high school)
# grocery spend 1.024 per EUR, f_price 0.88 (n.s.), f_service 2.99

b_logit <- pool_logit$pooled$estimate
exp(10 * b_logit[4])          # odds ratio 1.272 for +10 EUR grocery spend

# probit - coefficients (same signs and significance)
summary(pool_probit, conf.int = TRUE)

# logit coefficients are about 1.6-1.8 times the probit ones
b_probit <- pool_probit$pooled$estimate
compare <- data.frame(term   = pool_logit$pooled$term,
                      logit  = round(b_logit, 3),
                      probit = round(b_probit, 3),
                      ratio  = round(b_logit / b_probit, 2))
compare                       # ratios between 1.63 and 1.76


## ---- 5b. Predicted probabilities ----

# average of the fitted probabilities over the 20 data sets
p_logit  <- rep(0, nrow(data))
p_probit <- rep(0, nrow(data))
for (i in 1:20) {
  p_logit  <- p_logit  + fitted(fits_logit[[i]])  / 20
  p_probit <- p_probit + fitted(fits_probit[[i]]) / 20
}
cor(p_logit, p_probit)            # 0.9996
max(abs(p_logit - p_probit))      # 0.015 -> the link does not matter

plot(p_logit, p_probit, pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "P(upgrade) - logit", ylab = "P(upgrade) - probit",
     main = "Logit vs probit predictions")
abline(0, 1, col = "red")


## ---- 5c. ROC curve, AUC and threshold ----

roc_logit  <- roc(data$premium_upgrade, p_logit, quiet = TRUE)
roc_probit <- roc(data$premium_upgrade, p_probit, quiet = TRUE)
auc(roc_logit)                    # 0.769
ci.auc(roc_logit)                 # 95% CI 0.728 - 0.810
auc(roc_probit)                   # 0.769 (same)

# Youden threshold: maximises sensitivity + specificity
best <- coords(roc_logit, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
best                              # threshold 0.311, sens. 0.775, spec. 0.625
thr <- best$threshold[1]

plot(roc_logit, main = paste("ROC curve - AUC =", round(auc(roc_logit), 3)),
     legacy.axes = TRUE)
points(best$specificity[1], best$sensitivity[1], pch = 19, col = "red")

# confusion matrices at 0.5 and at the Youden threshold
y <- data$premium_upgrade
for (t in c(0.5, thr)) {
  pred_class <- as.numeric(p_logit >= t)
  TP <- sum(pred_class == 1 & y == 1)
  TN <- sum(pred_class == 0 & y == 0)
  FP <- sum(pred_class == 1 & y == 0)
  FN <- sum(pred_class == 0 & y == 1)
  print(table(predicted = pred_class, actual = y))
  print(round(c(threshold   = t,
                accuracy    = (TP + TN) / length(y),
                sensitivity = TP / (TP + FN),
                specificity = TN / (TN + FP),
                precision   = TP / (TP + FP)), 3))
}
# threshold 0.5  : accuracy 0.736, sensitivity 0.475, specificity 0.881
# threshold 0.311: accuracy 0.679, sensitivity 0.775, specificity 0.625


###############################################################################
## ---- 6. TASK 6 - Did the retention program work? ------------------------
###############################################################################

# all the variables used here are complete: no imputation needed

## ---- 6a. Raw numbers ----

table(loyalty = data$loyalty_program, churn = data$churn_event)
round(prop.table(table(data$loyalty_program, data$churn_event), 1), 3)
# churn_event = 1 for 53.9 % (not enrolled) and 35.5 % (enrolled)
# but this ignores censoring: the follow-up time is not the same for everybody

aggregate(months_observed ~ loyalty_program, data = data, FUN = mean)   # 6.6 / 7.6 months
table(loyalty = data$loyalty_program,
      censored_before_6 = data$churn_event == 0 & data$months_observed < 6)


## ---- 6b. Kaplan-Meier curves and log-rank test ----

km_fit <- survfit(Surv(months_observed, churn_event) ~ loyalty_program, data = data)
km_fit                                   # median 8.5 months (not enrolled), not reached (enrolled)
summary(km_fit, times = c(3, 6, 9, 11.8))
# still a customer at 6 months   : 64.1 % vs 77.7 %
# still a customer at 11.8 months: 33.5 % vs 56.4 %

plot(km_fit, col = c("red", "blue"), lwd = 2, conf.int = TRUE,
     xlab = "Months since start of observation", ylab = "Share still a customer",
     main = "Kaplan-Meier retention curves")
legend("bottomleft", c("Not enrolled", "Enrolled in program"),
       col = c("red", "blue"), lwd = 2, bty = "n")

survdiff(Surv(months_observed, churn_event) ~ loyalty_program, data = data)
# log-rank: chi2 = 21.3, df = 1, p < 0.001


## ---- 6c. Cox model and proportional hazards ----

cox <- coxph(Surv(months_observed, churn_event) ~ loyalty_program + age, data = data)
summary(cox)
# loyalty_program: HR = 0.545 (95% CI 0.42 - 0.70), p < 0.001
# age            : HR = 1.009 per year, p = 0.019

hr <- exp(coef(cox))
100 * (1 - hr[1])                        # enrolled customers leave at a rate 45.5 % lower
100 * (exp(10 * coef(cox)[2]) - 1)       # +9.7 % churn rate per 10 years of age

zph <- cox.zph(cox)
zph               # p = 0.15 (program), 0.52 (age), 0.31 (global)
                  # -> no evidence against proportional hazards

par(mfrow = c(1, 2))
plot(zph[1], main = "Schoenfeld residuals: loyalty_program")
abline(h = coef(cox)[1], col = "red", lty = 2)
plot(zph[2], main = "Schoenfeld residuals: age")
abline(h = coef(cox)[2], col = "red", lty = 2)
par(mfrow = c(1, 1))


## ---- 6d. Discrete time: did the customer leave within 6 months? ----

# Customers censored before month 6 are removed: we do not know whether
# they would have left before month 6 or not.
d6 <- data[!(data$churn_event == 0 & data$months_observed < 6), ]
d6$left6 <- as.numeric(d6$churn_event == 1 & d6$months_observed <= 6)
nrow(d6)                  # 463 customers kept (97 removed)
round(prop.table(table(d6$loyalty_program, d6$left6), 1), 3)

lg6 <- glm(left6 ~ loyalty_program + age, data = d6, family = binomial(link = "logit"))
pb6 <- glm(left6 ~ loyalty_program + age, data = d6, family = binomial(link = "probit"))
summary(lg6)$coefficients     # loyalty -0.677 (p = 0.001), age 0.0133 (p = 0.034)
summary(pb6)$coefficients     # loyalty -0.413 (p < 0.001), age 0.0081 (p = 0.032)

# odds ratios (logit)
round(exp(cbind(OR = coef(lg6), confint.default(lg6))), 3)   # program: OR = 0.51 (0.34 - 0.76)
round(coef(lg6) / coef(pb6), 2)                               # ratio logit / probit: about 1.64

# predicted probability of leaving within 6 months at the mean age
# (not enrolled / enrolled)
nd <- data.frame(loyalty_program = c(0, 1), age = mean(d6$age))
predict(lg6, nd, type = "response")       # logit : 0.381  0.238
predict(pb6, nd, type = "response")       # probit: 0.382  0.238
1 - summary(km_fit, times = 6)$surv       # Kaplan-Meier: 0.359  0.223


## ---- Session info ----
sessionInfo()
