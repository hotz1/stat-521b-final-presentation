# load packages
library(knitr)
library(kableExtra)
library(tidyverse)
library(here)

set.seed(1234)

# Define an "unknown" function to show off the GP regression model
f_tilde <- function(x){
  # Compute function
  fun_vals <-  0.125 * sin(pi*x)^2 * sqrt(x + 4) - sin(pi*x/2) * cos(pi*x^2/7)
  
  # Replace f(x) for any values outside of the range [0, 10] with NA
  fun_vals[which(!between(x, 0, 10))] <- NA
  return(fun_vals)
}

# Plot the true function 
f_tilde_plot <- ggplot() + 
  stat_function(fun = f_tilde, n = 501, color = "#ff9900",
                linewidth = 0.75, xlim = c(0,10)) + 
  labs(x = "x", y = latex2exp::TeX("\\tilde{f}(x)"),
       title = "Plot of the true \"unknown\" function") +
  scale_x_continuous(breaks = seq(0, 10, 2)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "#fafafa"))

# Find max of the function using optimize
true_optim <- optimize(f_tilde, interval = c(7,8), maximum = TRUE)

# Add a single point at the global optimum value to the plot
(f_tilde_plot <- f_tilde_plot +
  geom_point(aes(x = true_optim$maximum, y = true_optim$objective),
             fill = "blue", size = rel(2), shape = 23))
# ggsave(here("images/true_example_fn.jpg"), f_tilde_plot)

# squared exponential kernel function
sq_exp_kernel <- function(X1, X2, len = 1, prior_sd = 1){
  # Squared exponential kernel distance between X1 and X2
  # len is length parameter (scales the distances)
  # prior_se is the standard deviation of the functional values from GP prior
  exps <- outer(X1, X2, \(x1, x2) exp(-0.5 / len^2 * (x1-x2)^2))
  return(prior_sd^2 * exps)
}

# Mean and covariance matrix of posterior predictive dist. from GP regression 
posterior_pred <- function(new_x, datapoints, noise_sd = 1e-6, len = 1, prior_sd = 1){
  # new_x contains value(s) of x which we want to predict
  # datapoints is a collection of observed (x,y) pairs
  # noise_sd is the standard deviation of iid noise for observing y = f(x) + e
  # len and prior_sd are passed to the kernel function
  
  # Get observed data from the matrix
  obs_x <- datapoints[,1]
  obs_y <- datapoints[,2]
  
  # kernels
  obs_to_obs <- sq_exp_kernel(obs_x, obs_x, len, prior_sd) + noise_sd^2 * diag(nrow(datapoints))
  obs_to_new <- sq_exp_kernel(obs_x, new_x, len, prior_sd)
  new_to_new <- sq_exp_kernel(new_x, new_x, len, prior_sd) + noise_sd^2 * diag(length(new_x))
  
  # Posterior predictive mean and variance
  post_mean <- t(obs_to_new) %*% solve(obs_to_obs) %*% obs_y
  post_covmat <- new_to_new - t(obs_to_new) %*% solve(obs_to_obs) %*% obs_to_new
  return(list(mean = post_mean, cov = post_covmat))
}

# Function to plot GP estimates and queried points
plot_GP_regression <- function(GP_estimates, GP_observed){
  ggplot() +
    geom_line(aes(x = x, y = mean), data = GP_estimates, 
              colour = "black", linewidth = 1) +
    geom_ribbon(aes(x = x, ymin = mean - qnorm(0.975) * se, ymax = mean + qnorm(0.975) * se),
                data = GP_estimates, fill = "#92929240", color = "black", lty = "dashed") +
    geom_point(aes(x = x_obs, y = y_obs), data = GP_observed, 
               colour = "red", size = 2) +
    labs(x = "x", y = latex2exp::TeX("\\tilde{f}(x)"),
         title = latex2exp::TeX("95% Credible Intervals for \\tilde{f}(x)"),
         subtitle = bquote("Based on"~.(nrow(GP_observed))~"samples and a Gaussian Process with prior"~mu~"= 0,"~sigma^2~"= 1")) +
    scale_x_continuous(breaks = seq(0, 10, 2)) +   
    theme_minimal() +
    ylim(c(-2.5, 2.5)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "#fafafa"))
}

# Function to update estimated GP mean and standard error based on new observations
update_GP_estimates <- function(GP_ests, datapoints, noise_sd = 1e-6, len = 1, prior_sd = 1){
  posterior <- posterior_pred(GP_ests$x, datapoints, noise_sd, len, prior_sd)
  return(data.frame(x = GP_ests$x, mean = posterior$mean, se = sqrt(diag(posterior$cov))))
}

# Create matrix of x-values from 0 to 10 and estimates of the function
GP_estimates <- data.frame(x = seq(0, 10, 0.02), mean = 0, se = 1)

# Matrix of "observed" data
GP_observed <- data.frame(x_obs = runif(20, 0, 10), zscore = rnorm(20, 0, 1))
GP_observed <- GP_observed %>%
  mutate(y_obs_low = f_tilde(x_obs) + 0.025 * zscore, 
         y_obs_high = f_tilde(x_obs) + 0.25 * zscore)

# Plot the functions
for(n in c(1, 2, 5, 10, 20)){
  # Low noise
  # temp_low <- GP_observed %>% 
  #   filter(row_number() <= n) %>%
  #   select(x_obs, y_obs = y_obs_low)
  # 
  # plot_GP_regression(update_GP_estimates(GP_estimates, temp_low, 0.025), temp_low) %>%
  #   ggsave(filename = here(paste0("images/GP", n, "_noise_low.jpg")),
  #          width = 6, height = 4, units = "in")
  
  # High noise
  temp_high <- GP_observed %>% 
    filter(row_number() <= n) %>%
    select(x_obs, y_obs = y_obs_high)
  
  plot_GP_regression(update_GP_estimates(GP_estimates, temp_high, 0.25), temp_high) %>%
    ggsave(filename = here(paste0("images/GP", n, "_noise_high.jpg")),
           width = 6, height = 4, units = "in")
}

# Compute the three acquisition functions for the data
GP_updated_high20 <- update_GP_estimates(GP_estimates, temp_high, 0.25)
y_best <- max(temp_high$y_obs)
GP_updated_high20 <- GP_updated_high20 %>%
  mutate(gamma = (mean - y_best)/se, 
         prob_imp = pnorm(gamma),
         expec_imp = se * (gamma * pnorm(gamma) + dnorm(gamma)),
         GPUCB = mean + se * (2*log(1^(1/2 + 2) * pi^2 * 3/0.05))^(1/2))

(acquisition_plots <- GP_updated_high20 %>% 
  pivot_longer(cols = c(prob_imp, expec_imp, GPUCB),
               names_to = "Func", values_to = "Acquisition") %>%
  mutate(Func = factor(Func, levels = c("prob_imp", "expec_imp", "GPUCB"),
                       labels = c("Probability of Improvement", 
                                  "Expected Improvement", "GP-UCB"))) %>%
  ggplot() +
  geom_line(aes(x = x, y = Acquisition, color = Func), linewidth = 0.75) +
  facet_wrap(~Func, nrow = 3, scales = "free_y") +
  labs(y = "Acquisition Function", title = "Values of the three primary GP acquisition functions") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 10, 2)) +   
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "#fafafa"),
        strip.background = element_rect(fill = "#FFDF0040"),
        legend.position = "none"))
 
# ggsave(here("images/acquisition_plots.jpg"), acquisition_plots,
#        width = 6, height = 4, units = "in")
        
