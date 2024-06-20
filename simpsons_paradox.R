library(tibble)
library(tidybayes)
library(tidyr)
library(dplyr)
library(brms)
library(psych)
library(faintr)

#### helper functions

summarize_posterior <- function(post_do0, post_do1, 
                                label_do0 = "refuse drug",
                                label_do1 = "take drug") {
  
  # Ensure the inputs are tibbles (this is ugly but necessary to ensure the 
  # function can handle tibbles and normal lists as input
  post_do0 <- tibble(value = pull(as_tibble(post_do0)))
  post_do1 <- tibble(value = pull(as_tibble(post_do1)))
  
  #calculate difference between conditions
  diff <- post_do1 - post_do0
  
  # combine data into single tibble
  combined_data <- bind_rows(
    post_do1 |> mutate(condition = label_do1),
    post_do0 |> mutate(condition = label_do0),
    diff |> mutate(condition = "causal effect")
  )
  
  # calculate means and confidence intervals
  post_sum <- combined_data |>
    group_by(condition) |>
    summarize(
      CI_lower = quantile(value, 0.025),
      mean = mean(value),
      CI_upper = quantile(value, 0.975)
    )
  
  # reorder columns
  order <- c(label_do1, label_do0, "causal effect")
  post_sum <- post_sum |>
    mutate(condition = factor(condition, levels = order)) |>
    arrange(condition)
  
  return(post_sum)
  
}




# generate dataset with Simpsons paradox
data_simpsons_paradox <- tibble(
  gender = c("Male", "Male", "Female", "Female"),
  bloodP = c("Low", "Low", "High", "High"),
  drug = c("Take", "Refuse", "Take", "Refuse"),
  k = c(81, 234, 192, 55),
  N = c(87, 270, 263, 80),
  proportion = k/N
)

# cast into long format
data_SP_long <- rbind(
  data_simpsons_paradox |> 
    uncount(k) |>
    mutate(recover = TRUE) |>
    select(-N, -proportion), 
  data_simpsons_paradox |>
    uncount(N-k) |>
    mutate(recover=FALSE) |>
    select(-k, -N, -proportion)
  )


#### Causal inference with gender as confound
n_iter = 2000

# fist, estimate the distribution of gender (by intercept-only regression)
# we need this later when we marginalize over gender
fit_gender <- brm(
  formula = gender ~ 1,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

# then estimate the distribution of recovery rates as predicted
# by gender and treatment (this is before the do-calculus step)
fit_recovery_gender <- brm(
  formula = recover ~ gender * drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)



# sample from estimated gender dist
# why? because for the do(treatment) operation, we make gender 
# independent from treatment (remove causal connection)
# in the formula: we no longer compute P(G=g|D=d), but only 
# P(G=g)
posterior_gender_sample <- tidybayes::predicted_draws(
  # predicted_draws draws from posterior predictive dist
  # P(y_new | x_new, y_obs), because we want data points
  object = fit_gender,
  newdata = tibble(Intercept = 1),
  value = "gender",
  ndraws = n_iter*2
) |> 
  ungroup() |>
  mutate(gender = ifelse(gender, "Male", "Female")) |>
  select(gender)



# posterior predictive samples for D=1 (do(D=1))
posterior_DrugTaken_g <- tidybayes::epred_draws(
  # epred_draws draws from the expectation of the post pred dist 
  # E(y_new | x_new, y_obs), because we want the mean
  object = fit_recovery_gender,
  newdata = posterior_gender_sample |> mutate(drug="Take"),
  value = "taken", #can we just leave this argument out and call them all "value"?
  ndraws = n_iter * 2
) |> ungroup() |>
  select(taken)

# posterior predictive samples for D=0
posterior_DrugRefused_g <- tidybayes::epred_draws(
  object = fit_recovery_gender,
  newdata = posterior_gender_sample |> mutate(drug="Refuse"),
  value = "refused",
  ndraws = n_iter * 2
) |> ungroup() |>
  select(refused)

# summarize results 
summarize_posterior(posterior_DrugRefused_g, posterior_DrugTaken_g)



#### Causal inference with blood pressure as mediator
# in this formula, we ignore bp alltogether. We assume that bp is either fully 
# determined by drug intake, therefore all changes in bp depend on the drug
# and we can ignore the effect of bp on recovery. Or, if bp is not only dependent
# on drug intake, and has other causes, we assume these causes to be independent
# from drug intake, and therefore equally distributed between the drug and placebo
# groups. This means that if we find a significant difference of recovery rates 
# between the drug and placebo groups, this difference can only be explained by 
# drug intake
fit_recovery_bp <- brm(
  formula = recover ~ drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

posterior_DrugTaken_bp <-
  faintr::filter_cell_draws(fit_recovery_bp, drug == "Take") |>
  pull(draws) |>
  logistic()


posterior_DrugRefused_bp <-
  faintr::filter_cell_draws(fit_recovery_bp, drug == "Refuse") |>
  pull(draws) |>
  logistic()














