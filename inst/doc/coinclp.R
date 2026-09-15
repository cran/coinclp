## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(coinclp)

## ----version------------------------------------------------------------------
clp_version()$version

## ----solve--------------------------------------------------------------------
A <- rbind(material = c(120, 210),
           labour   = c(110,  30),
           capacity = c(  1,   1))
b <- c(15000, 4000, 75)

fit <- clp_solve(c(143, 60), A, "<=", b, max = TRUE,
                 col_names = c("x", "y"), row_names = rownames(A))
fit

## ----solution-----------------------------------------------------------------
fit$solution
fit$duals
fit$reduced_costs

## ----activity-----------------------------------------------------------------
data.frame(row = rownames(A), activity = fit$row_activity, limit = b,
           dual = fit$duals)

## ----ranged-------------------------------------------------------------------
clp_solve(c(1, 1), rbind(c(1, 1)),
          row_lower = 2, row_upper = 5,
          lower = -Inf, upper = Inf)$objval

## ----sparse, eval = requireNamespace("Matrix", quietly = TRUE)----------------
set.seed(1)
big <- Matrix::rsparsematrix(200, 500, density = 0.01)
res <- clp_solve(rep(1, 500), big, ">=", rep(-1, 200), lower = 0, upper = 10)
res$status_message

## ----model--------------------------------------------------------------------
model <- clp_model()
clp_set_log_level(model, 0)
clp_load_problem(model, ncols = 2, nrows = 3,
                 start = c(0L, 3L, 6L),
                 index = c(0L, 1L, 2L, 0L, 1L, 2L),
                 value = c(120, 110, 1, 210, 30, 1),
                 obj   = c(-143, -60),
                 rowub = b)
clp_initial_solve(model)
c(objective = clp_objective_value(model), iterations = clp_iterations(model))

## ----warmstart----------------------------------------------------------------
basis <- clp_status_array(model)

clp_set_row_upper(model, c(15000, 4000, 70))
clp_copyin_status(model, basis)
clp_dual_simplex(model)

c(objective = clp_objective_value(model), iterations = clp_iterations(model))

## ----algorithms---------------------------------------------------------------
for (alg in c("auto", "primal", "dual", "barrier")) {
    fit <- clp_solve(c(143, 60), A, "<=", b, max = TRUE,
                     control = clp_control(algorithm = alg))
    cat(sprintf("%-8s %.4f\n", alg, fit$objval))
}

## ----mps----------------------------------------------------------------------
path <- system.file("extdata", "productmix.mps", package = "coinclp")
from_file <- clp_model()
clp_set_log_level(from_file, 0)
clp_read_mps(from_file, path)
clp_col_names(from_file)
clp_initial_solve(from_file)
clp_objective_value(from_file)

## ----features-----------------------------------------------------------------
clp_features()

## ----cleanup, include = FALSE-------------------------------------------------
clp_free(model)
clp_free(from_file)

## ----compat-------------------------------------------------------------------
lp <- initProbCLP()
setLogLevelCLP(lp, 0)
loadProblemCLP(lp, 2, 3, c(0, 3, 6), c(0, 1, 2, 0, 1, 2),
               c(120, 110, 1, 210, 30, 1),
               lb = c(0, 0), ub = c(1e30, 1e30), obj_coef = c(143, 60),
               rlb = rep(-1e30, 3), rub = b)
setObjDirCLP(lp, -1)
solveInitialCLP(lp)
status_codeCLP(getSolStatusCLP(lp))
getObjValCLP(lp)
delProbCLP(lp)

