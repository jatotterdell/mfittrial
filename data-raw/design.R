library(bayestestR)

# Full Model
int <- 4
pref <- 4
com <- int * pref
int_lab <- factor(rep(1:4, each = 4))
pref_lab <- factor(rep(1:4, times = 4))

Xint  <- kronecker(diag(1, int), rep(1, pref))
Xpref <- kronecker(rep(1, int), diag(1, pref))
Xcom  <- cbind(Xpref[, 1] * Xint, Xpref[, 2] * Xint, Xpref[, 3] * Xint, Xpref[, 4] * Xint)
X <- cbind(1, Xint, Xpref, Xcom)
colnames(X) <- c("intercept", paste0("int", 1:int), paste0("pref", 1:pref),
                 paste0(rep(paste0("int", 1:int), each = pref),
                        rep(paste0("pref", 1:pref), times = int)))
rownames(X) <- paste0(rep(paste0("int", 1:int), each = pref),
                      rep(paste0("pref", 1:pref), times = int))
X_con <- model.matrix( ~ int_lab*pref_lab, contrasts = list(int_lab = "contr.bayes", pref_lab = "contr.bayes"))[,]
X_con_t <- t(X_con)

Qint <- contr.bayes(4)
Qpref <- contr.bayes(4)
Qcom <- kronecker(Qint, Qpref)
Q <- as.matrix(Matrix::bdiag(1, Qint, Qpref, Qcom))
rownames(Q) <- colnames(X)
Q_t <- t(Q)

# Reduce model with no interaction
X_con_red <- model.matrix( ~ int_lab + pref_lab, contrasts = list(int_lab = "contr.bayes", pref_lab = "contr.bayes"))[,]
X_con_red_t <- t(X_con_red)
Q_red <- as.matrix(Matrix::bdiag(1, Qint, Qpref))
rownames(Q_red) <- colnames(X)[1:9]
Q_red_t <- t(Q_red)

# Equivalently
# D <- expand.grid(pref = factor(1:4), int = factor(1:4))
# X_con <- model.matrix( ~ pref * int, data = D, contrasts = list(pref = contr.bayes, int = contr.bayes))[,]



# Reduced form where we only care if they did or did not get their preference
#=========================================================================
# D <- expand.grid(pref = factor(1:4), int = factor(1:4))
# D$have_pref <- factor(as.numeric(D$pref != 1))
# D$got_pref <- factor(as.numeric(D$pref == D$int & D$int != 1))
# Xgot_pref <- cbind(got_pref = model.matrix( ~ 0 + got_pref, D))
# Xgot_pref[1:4, 1] <- 0
# X_red <- cbind(
#   intercept = 1,
#   model.matrix( ~ 0 + int, D),
#   Xgot_pref,
#   Xint[, 1] * Xgot_pref, Xint[, 2] * Xgot_pref, Xint[, 3] * Xgot_pref, Xint[, 4] * Xgot_pref)
# colnames(X_red)[8:15] <- paste0(rep(paste0("int", 1:int), each = 2),
#                                 rep(paste0("got_pref", 1:2), times = int))
# X_red <- unique(X_red)
# rownames(X_red) <- 1:nrow(X_red)
#
# Xgot_pref_con <- cbind(got_pref = model.matrix( ~ got_pref, D, contrasts = list(got_pref = contr.bayes))[, 2])
# X_red_con <- cbind(
#   intercept = 1,
#   Xint_con,
#   Xgot_pref_con,
#   Xint_con[, 1] * Xgot_pref_con, Xint_con[, 2] * Xgot_pref_con, Xint_con[, 3] * Xgot_pref_con
# )
# colnames(X_red_con)[-1] <- c("int1", "int2", "int3", "got_pref", paste0("int", 1:3, "got_pref"))
# X_red_con <- unique(X_red_con)
# rownames(X_red_con) <- 1:nrow(X_red_con)
#
# X_red_con_t <- t(X_red_con)
#
# Qgot_pref <- contr.bayes(2)
# Qredcom <- kronecker(Qint, Qgot_pref)
# Q_red <- as.matrix(Matrix::bdiag(1, Qint, Qgot_pref, Qredcom))
# Q_red_t <- t(Q_red)

usethis::use_data(X_con, X_con_t, Q, Q_t, X_con_red, X_con_red_t, Q_red, Q_red_t, overwrite = T)
