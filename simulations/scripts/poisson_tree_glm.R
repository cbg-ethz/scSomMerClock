library('addreg')
library('glm2')

LAMBDA_MIN <- 1e-6
df_H0 <- read.csv('/home/hw/Desktop/molClock_project/MolClockAnalysis/simulations/glm_poisson_data.csv')
y_no <- dim(df_H0)[1]
par_no <- dim(df_H0)[2]
df_H0$Y <- df_H0$Y + 1

init <- rep(mean(df_H0$Y), par_no - 1)

X_H0 <- data.matrix(df_H0[,2:par_no])
df_H1 <- data.frame(diag(y_no))
df_H1$Y <- df_H0$Y

H0_formula <- paste('Y ~ -1 ', paste(names(df_H0)[2:par_no], collapse = ' + '), sep = ' + ')
H0_formula_intercept <- paste('Y ~ ', paste(names(df_H0)[2:par_no], collapse = ' + '), sep = ' ')
H1_formula_intercept <- paste('Y ~ ', paste(names(df_H1)[1:y_no], collapse = ' + '), sep = ' ')

# Standard glm with poisson and identity link - problem of negative mu
model <- glm(H0_formula, family=poisson(link = "identity"), data=df, start = init)
estimates = X_H0 %*%pmax(model$coefficients, LAMBDA_MIN)
#summary(model)

ll_H0 <- sum(dpois(df_H0$Y, estimates, log = TRUE))
ll_H1 <- sum(dpois(df_H0$Y, pmax(df_H0$Y, LAMBDA_MIN), log = TRUE))

LR <- -2 * (ll_H0 - ll_H1)
dof <- y_no - model$rank
p_val <- pchisq(LR, dof, lower.tail=FALSE)
print(paste('P-value glm: ', p_val))

#  Improved glm2 with poisson and identity link
model2 <- glm2(H0_formula, family=poisson(link = "identity"), data=df_H0, start=init, maxit = 1000)
estimates2 = X_H0 %*% pmax(model2$coefficients, LAMBDA_MIN)
#summary(model2)

ll_H0_2 <- sum(dpois(df_H0$Y, estimates2, log = TRUE))
LR_2 <- -2 * (ll_H0_2 - ll_H1)
dof_2 <- y_no - model2$rank
p_val_2 <- pchisq(LR_2, dof_2, lower.tail=FALSE)
print(paste('P-value glm2: ', p_val_2))

# Non negative poisson regression using an EM algorithm
model_nn_H0 <- addreg(H0_formula_intercept, family=poisson(link = "identity"), data=df_H0,
                   method='em', accelerate = "squarem", maxit = 1000)
estimates_nn_H0 = X_H0 %*% pmax(model_nn$coefficients[2:par_no] + model_nn$coefficients[1], LAMBDA_MIN)

model_nn_H1 <- addreg(H1_formula_intercept, family=poisson(link = "identity"), data=df_H1,
                   method='em', accelerate = "squarem", maxit = 1000)
estimates_nn_H1 = pmax(model_nn_H1$coefficients[2:(y_no + 1)] + model_nn_H1$coefficients[1], LAMBDA_MIN)

# summary(model_nn)
ll_H0_nn <- sum(dpois(df_H0$Y, as.vector(estimates_nn_H0), log = TRUE))
ll_H1_nn <- sum(dpois(df_H0$Y, as.vector(estimates_nn_H1), log = TRUE))

LR_nn <- -2 * (ll_H0_nn - ll_H1)
dof_nn <- model_nn_H1$rank - model_nn_H0$rank
p_val_nn <- pchisq(LR_nn, dof_nn, lower.tail=FALSE)
print(paste('P-value addreg: ', p_val_nn))

# negative binomial model
phi <- 0.5

model_nbinom_nn_H0 <- addreg(H0_formula_intercept, family=negbin1(link = "identity", phi = phi), data=df_H0,
                             method='em', accelerate = "squarem", maxit = 1000)
estimates_nbinom_H0 = X_H0 %*% pmax(model_nbinom_nn_H0$coefficients[2:par_no] + model_nn$coefficients[1], LAMBDA_MIN)
prob_H0 <- 1 / (model_nbinom_nn_H0$scale + 1)

model_nbinom_nn_H1 <- addreg(H1_formula_intercept, family=negbin1(link = "identity", phi = phi), data=df_H1,
                             method='em', accelerate = "squarem", maxit = 1000)
estimates_nbinom_H1 = pmax(model_nbinom_nn_H1$coefficients[2:(y_no + 1)] + model_nn_H1$coefficients[1], LAMBDA_MIN)

ll_H0_nbinom <- sum(dnbinom(df_H0$Y, prob=prob_H0, estimates_nbinom_H0, log = TRUE))
ll_H1_nbinom <- sum(dnbinom(df_H0$Y, prob=prob_H1, estimates_nbinom_H1, log = TRUE))

LR_nbinom <- -2 * (ll_H0_nbinom - ll_H1_nbinom)
dof_nbinom <- model_nbinom_nn_H1$rank - model_nbinom_nn_H0$rank
p_val_nbinom <- pchisq(LR_nbinom, dof_nbinom, lower.tail=FALSE)
print(paste('P-value addreg: ', p_val_nbinom))
