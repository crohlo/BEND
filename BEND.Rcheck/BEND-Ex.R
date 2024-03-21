pkgname <- "BEND"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BEND-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BEND')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Bayes_BPREM")
### * Bayes_BPREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Bayes_BPREM
### Title: Bayesian Bivariate Piecewise Random Effects Model (BPREM)
### Aliases: Bayes_BPREM

### ** Examples

## No test: 
# load simulated data
data(SimData_BPREM)
# plot observed data
plot_BEND(data = SimData_BPREM,
          id_var = "id",
          time_var = "time",
          y_var = "y1",
          y2_var = "y2")
# fit Bayes_BPREM()
results_bprem <- Bayes_BPREM(data = SimData_BPREM,
                             id_var = "id",
                             time_var = "time",
                             y1_var = "y1",
                             y2_var = "y2")
# result summary
summary(results_bprem)
# plot fitted results
plot_BEND(data = SimData_BPREM,
          id_var = "id",
          time_var = "time",
          y_var = "y1",
          y2_var = "y2",
          results = results_bprem)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Bayes_BPREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Bayes_CREM")
### * Bayes_CREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Bayes_CREM
### Title: Bayesian Crossed Random Effects Model (CREM)
### Aliases: Bayes_CREM

### ** Examples

## No test: 
# load simulated data
data(SimData_PCREM)
# plot observed data
plot_BEND(data = SimData_PCREM,
          id_var = "id",
          time_var = "time",
          y_var = "y")
# fit Bayes_CREM()
results_pcrem <- Bayes_CREM(data = SimData_PCREM,
                            ind_id_var = "id",
                            cross_id_var = "teacherid",
                            time_var = "time",
                            y_var = "y",
                            form="piecewise")
# result summary
summary(results_pcrem)
# plot fitted results
plot_BEND(data = SimData_PCREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_pcrem)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Bayes_CREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Bayes_PREM")
### * Bayes_PREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Bayes_PREM
### Title: Bayesian Piecewise Random Effects Model (PREM) + Extensions
### Aliases: Bayes_PREM

### ** Examples

## No test: 
# load simulated data
data(SimData_PREM)
# plot observed data
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y")

# PREM ---------------------------------------------------------------------------------
# fit Bayes_PREM()
results_prem <- Bayes_PREM(data = SimData_PREM,
                           id_var = "id",
                           time_var = "time",
                           y_var = "y")
# result summary
summary(results_prem)
# plot fitted results
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_prem)

# CI-PREM ---------------------------------------------------------------------------------
# fit Bayes_PREM()
results_ciprem <- Bayes_PREM(data = SimData_PREM,
                             id_var = "id",
                             time_var = "time",
                             y_var = "y",
                             outcome_predictive_vars = "outcome_pred_1")
# result summary
summary(results_ciprem)
# plot fitted results
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_ciprem)

# PREMM ---------------------------------------------------------------------------------
# fit Bayes_PREM()
results_premm <- Bayes_PREM(data = SimData_PREM,
                            id_var = "id",
                            time_var = "time",
                            y_var = "y",
                            n_class = 2)
# result summary
summary(results_premm)
# plot fitted results
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_premm)


# CI-PREMM ---------------------------------------------------------------------------------
# fit Bayes_PREM()
results_cipremm <- Bayes_PREM(data = SimData_PREM,
                              id_var = "id",
                              time_var = "time",
                              y_var = "y",
                              n_class = 2,
                              class_predictive_vars = c("class_pred_1", "class_pred_2"),
                              outcome_predictive_vars = "outcome_pred_1")
# result summary
summary(results_cipremm)
# plot fitted results
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_cipremm)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Bayes_PREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_BEND")
### * plot_BEND

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_BEND
### Title: Plot a BEND Model (PREM, CREM, BPREM)
### Aliases: plot_BEND

### ** Examples

# load simulated data
data(SimData_PREM)
# plot observed data
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y")
# load fitted model results
data(results_prem)
# plot fitted results
plot_BEND(data = SimData_PREM,
          id_var = "id",
          time_var = "time",
          y_var = "y",
          results = results_prem)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_BEND", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.BPREM")
### * summary.BPREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.BPREM
### Title: Summarize the results of a bivariate piecewise random effects
###   model (BPREM)
### Aliases: summary.BPREM

### ** Examples

# load fitted model results
data(results_bprem)
# result summary
summary(results_bprem)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.BPREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.CREM")
### * summary.CREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.CREM
### Title: Summarize the results of a crossed random effects model (CREM)
### Aliases: summary.CREM

### ** Examples

# load fitted model results
data(results_pcrem)
# result summary
summary(results_pcrem)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.CREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.PREM")
### * summary.PREM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.PREM
### Title: Summarize the results of a piecewise random effects model (PREM)
### Aliases: summary.PREM

### ** Examples

# load fitted model results
data(results_prem)
# result summary
summary(results_prem)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.PREM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
