% Type I
scenario_null <- list(
    n_clusters = 100,
    cluster_size = 5,
    beta = c(0.5, 0.2),
    Gamma = matrix(0, 2, 2),
    theta = 2
)

sim_null <- run_lrt_simulation(
    scenario = scenario_null,
    R = 100,
    alpha = 0.05,
    weights = c("oracle", "plugin"),
    seed = 123
)

sim_null$summary
sim_null$diagnostics


% Power
scenario_alt <- list(
    n_clusters = 100,
    cluster_size = 5,
    beta = c(0.5, 0.2),
    Gamma = matrix(
        c(0.5, 0.1,
          0.1, 0.3),
        nrow = 2
    ),
    theta = 2
)

sim_alt <- run_lrt_simulation(
    scenario = scenario_alt,
    R = 100,
    alpha = 0.05,
    weights = c("oracle", "plugin"),
    seed = 456
)

sim_alt$summary
