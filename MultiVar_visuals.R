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

#-------------------------------------------------------------------------------
# Behavior as n increases, heterogeneous
#-------------------------------------------------------------------------------
# Weak Effect Size, heterogeneous
# Simulated EDRs
incr_simEDR.ht.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = c(0,0,0,0,0),
      exp_mu = c(0.2,0.4,0.6,0.8,1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
incr_trueEDR.ht.n <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = c(0.2,0.4,0.6,0.8,1),
      wgt = rep(0.2, 5)
    )
  ),
  type = "Theoretical"
)

incr_EDR.ht.n <- bind_rows(incr_simEDR.ht.n, incr_trueEDR.ht.n)

nplot_incr.ht.n <- ggplot(incr_EDR.ht.n, aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR with Heterogenous ES = [0.2, 1]",
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Weak Effect Size, heterogeneous
# Simulated ERRs
incr_simERR.ht.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    err_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = c(0,0,0,0,0),
      exp_mu = c(0.2,0.4,0.6,0.8,1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
incr_trueERR.ht.n <- tibble(
  n = n_vals,
  err = sapply(n_vals, function(n)
    true_err(
      n = n,
      alpha = .05,
      pop.es = c(0.2,0.4,0.6,0.8,1),
      wgt = rep(0.2, 5)
    )
  ),
  type = "Theoretical"
)

incr_ERR.ht.n <- bind_rows(incr_simERR.ht.n, incr_trueERR.ht.n)

nplot_incr.ht.n <- ggplot(incr_ERR.ht.n, aes(x = n, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR with Heterogenous ES = [0.2, 1]",
    x = "Sample size (n)",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Behavior as correlation r increases
#-------------------------------------------------------------------------------
# Sequence of r values
r_vals <- seq(0, 1, by = 0.1)

# Under Null Hypothesis
# Simulated EDRs
null_simEDR.hm.r <- tibble(
  r = r_vals,
  edr = sapply(r_vals, function(r)
    edr_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = rep(0, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
null_trueEDR.hm.r <- tibble(
  r = r_vals,
  edr = rep(true_edr(
    n = 20,
    alpha = .05,
    pop.es = rep(0, 5),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

null_EDR.hm.r <- bind_rows(null_simEDR.hm.r, null_trueEDR.hm.r)

nplot_null.EDR.hm.r <- ggplot(null_EDR.hm.r,
                              aes(x = r, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR as Correlation r increases (homogenous)",
    x = "Correlation between DVs r = [0, 1]",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Under Null Hypothesis
# Simulated EDRs
incr_simEDR.hm.r <- tibble(
  r = r_vals,
  edr = sapply(r_vals, function(r)
    edr_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = c(0.2, 0.4, 0.6, 0.8, 1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
incr_trueEDR.hm.r <- tibble(
  r = r_vals,
  edr = rep(true_edr(
    n = 20,
    alpha = .05,
    pop.es = c(0.2, 0.4, 0.6, 0.8, 1),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

incr_EDR.hm.r <- bind_rows(incr_simEDR.hm.r, incr_trueEDR.hm.r)

plot_incr.EDR.hm.r <- ggplot(incr_EDR.hm.r,
                              aes(x = r, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR vs Correlation r ES = [0.2, 1]",
    x = "Correlation between DVs r = [0, 1]",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

plot_incr.EDR.hm.r

#------------------------                       -------------------------
# Med effect size
# Simulated EDRs
med_simEDR.hm.r <- tibble(
  r = r_vals,
  edr = sapply(r_vals, function(r)
    edr_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = c(0.2, 0.4, 0.6, 0.8, 1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
med_trueEDR.hm.r <- tibble(
  r = r_vals,
  edr = rep(true_edr(
    n = 20,
    alpha = .05,
    pop.es = c(0.2, 0.4, 0.6, 0.8, 1),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

med_EDR.hm.r <- bind_rows(med_simEDR.hm.r, med_trueEDR.hm.r)

plot_med.EDR.hm.r <- ggplot(med_EDR.hm.r,
                             aes(x = r, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "EDR vs Correlation r ES = 0.5 (homogenous)",
    x = "Correlation between DVs r = [0, 1]",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

plot_med.EDR.hm.r

#-------------------------------------------------------------------------------
# Under Null Hypothesis
# Simulated ERRs
null_simERR.hm.r <- tibble(
  r = r_vals,
  err = sapply(r_vals, function(r)
    err_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = rep(0, 5),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
null_trueERR.hm.r <- tibble(
  r = r_vals,
  err = rep(true_err(
    n = 20,
    alpha = .05,
    pop.es = rep(0, 5),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

null_ERR.hm.r <- bind_rows(null_simERR.hm.r, null_trueERR.hm.r)

plot_null.ERR.hm.r <- ggplot(null_ERR.hm.r,
                              aes(x = r, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR as Correlation r increases (homogenous)",
    x = "Correlation between DVs r = [0, 1]",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

plot_null.ERR.hm.r

#-------------------------------------------------------------------------------
# Under Null Hypothesis
# Simulated ERRs
incr_simERR.hm.r <- tibble(
  r = r_vals,
  err = sapply(r_vals, function(r)
    err_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = c(0.2, 0.4, 0.6, 0.8, 1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
incr_trueERR.hm.r <- tibble(
  r = r_vals,
  err = rep(true_err(
    n = 20,
    alpha = .05,
    pop.es = c(0.2, 0.4, 0.6, 0.8, 1),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

incr_ERR.hm.r <- bind_rows(incr_simERR.hm.r, incr_trueERR.hm.r)

plot_incr.ERR.hm.r <- ggplot(incr_ERR.hm.r,
                             aes(x = r, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR vs Correlation r ES = [0.2, 1]",
    x = "Correlation between DVs r = [0, 1]",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------
# Med effect size
# Simulated ERRs
med_simERR.hm.r <- tibble(
  r = r_vals,
  err = sapply(r_vals, function(r)
    err_m.sim(
      k_sims = 15000,
      n = 20,
      r = r,
      control_mu = rep(0, 5),
      exp_mu = c(0.2, 0.4, 0.6, 0.8, 1),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True ERRs
med_trueERR.hm.r <- tibble(
  r = r_vals,
  err = rep(true_err(
    n = 20,
    alpha = .05,
    pop.es = c(0.2, 0.4, 0.6, 0.8, 1),
    wgt = c(1,0,0,0,0)
  ), length(r_vals)),
  type = "Theoretical"
)

med_ERR.hm.r <- bind_rows(med_simERR.hm.r, med_trueERR.hm.r)

plot_med.ERR.hm.r <- ggplot(med_ERR.hm.r,
                            aes(x = r, y = err, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "ERR vs Correlation r ES = 0.5 (homogenous)",
    x = "Correlation between DVs r = [0, 1]",
    y = "ERR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))
