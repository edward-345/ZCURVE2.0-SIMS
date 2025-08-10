source('helper_functions.R')
source('MultiVar_sims.R')

#-------------------------------------------------------------------------------
# Behavior as n increases
#-------------------------------------------------------------------------------
# Sequence of n values
n_vals <- seq(10, 100, by = 10)

# Under Null Hypothesis
# Simulated EDRs
sim_edr <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    edr_m.sim(
      k_sims = 15000,
      n = n,
      r = 0.5,
      control_mu = c(0,0,0,0,0),
      exp_mu = c(0,0,0,0,0),
      sd = 1
    )
  ),
  type = "Simulated"
)

# True EDRs
true_edr <- tibble(
  n = n_vals,
  edr = sapply(n_vals, function(n)
    true_edr(
      n = n,
      alpha = .05,
      pop.es = numeric(5),
      wgt = c(0.2,0.2,0.2,0.2,0.2)
    )
  ),
  type = "Theoretical"
)

df <- bind_rows(df_sim, df_theory)

ggplot(df, aes(x = n, y = edr, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Sample size (n)",
    y = "EDR",
    color = "Curve type"
  ) +
  scale_color_manual(values = c(Simulated = "blue", Theoretical = "red"))

#-------------------------------------------------------------------------------