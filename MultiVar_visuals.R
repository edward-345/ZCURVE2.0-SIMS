source('helper_functions.R')
source('MultiVar_sims.R')

#-------------------------------------------------------------------------------
# Behavior as n increases
#-------------------------------------------------------------------------------
# Sequence of n values
n_vals <- seq(10, 100, by = 10)

# Under Null Hypothesis
# Simulated EDRs
null_simEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0, 5),
      exp_mu = rep(0, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
null_trueEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = rep(0, 5),
      wgt = c(1,0,0,0,0)
    )
  ),
  type = "Theoretical"
)

null_EDR.hm.n <- bind_rows(null_simEDR.hm.n, null_trueEDR.hm.n)

nplot_null.EDR.hm.n <- ggplot(null_EDR.hm.n,
                              aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Weak Effect Size, homogenous
# Simulated EDRs
weak_simEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = c(0,0,0,0,0),
      exp_mu = c(0.2,0.2,0.2,0.2,0.2),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
weak_trueEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = rep(0.2, 5),
      wgt = rep(0.2, 5)
    )
  ),
  type = "Theoretical"
)

weak_EDR.hm.n <- bind_rows(weak_simEDR, weak_trueEDR)

nplot_null.hm.n <- ggplot(weak_EDR.hm.n, aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR with Homogenous ES = 0.2",
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Medium Effect Size, homogenous
# Simulated EDRs
med_simEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0,5),
      exp_mu = rep(0.5, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
med_trueEDR.hm.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = rep(0.5, 5),
      wgt = rep(0.2, 5)
    )
  ),
  type = "Theoretical"
)

med_EDR.hm.n <- bind_rows(med_simEDR, med_trueEDR)

nplot_med.hm.n <- ggplot(med_EDR.hm.n, aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR with Homogenous ES = 0.5",
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Strong Effect Size, homogenous
# Simulated EDRs
str_simEDR.hm <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0,5),
      exp_mu = rep(0.8, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
str_trueEDR.hm <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = c(0.8,0,0,0,0),
      wgt = c(1,0,0,0,0)
    )
  ),
  type = "Theoretical"
)

str_EDR.hm <- bind_rows(str_simEDR.hm, str_trueEDR.hm)

nplot_str.hm <- ggplot(str_EDR.hm, aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR with Homogenous ES = 0.8",
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Under Null Hypothesis
# Simulated ERRs
null_simERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    err_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0, 5),
      exp_mu = rep(0, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
null_trueERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    true_err(
      n = n,
      alpha = .05,
      pop.es = rep(0, 5),
      wgt = c(1,0,0,0,0)
    )
  ),
  type = "Theoretical"
)

null_ERR.hm.n <- bind_rows(null_simERR.hm.n, null_trueERR.hm.n)

nplot_null.ERR.hm.n <- ggplot(null_ERR.hm.n,
                              aes(x = n, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR with Homogenous Null",
    x = "Sample size (n)",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Weak Effect size, homogenous
# Simulated ERRs
weak_simERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    err_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0, 5),
      exp_mu = rep(0.2, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
weak_trueERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    true_err(
      n = n,
      alpha = .05,
      pop.es = c(0.2,0,0,0,0),
      wgt = c(1,0,0,0,0)
    )
  ),
  type = "Theoretical"
)

weak_ERR.hm.n <- bind_rows(weak_simERR.hm.n, weak_trueERR.hm.n)

nplot_weak.ERR.hm.n <- ggplot(weak_ERR.hm.n,
                              aes(x = n, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR with Homogenous ES = 0.2",
    x = "Sample size (n)",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Medium Effect size, homogenous
# Simulated ERRs
med_simERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    err_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = rep(0, 5),
      exp_mu = rep(0.5, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
med_trueERR.hm.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    true_err(
      n = n,
      alpha = .05,
      pop.es = c(0.5,0,0,0,0),
      wgt = c(1,0,0,0,0)
    )
  ),
  type = "Theoretical"
)

med_ERR.hm.n <- bind_rows(med_simERR.hm.n, med_trueERR.hm.n)

nplot_med.ERR.hm.n <- ggplot(med_ERR.hm.n,
                              aes(x = n, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR with Homogenous ES = 0.5",
    x = "Sample size (n)",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))
