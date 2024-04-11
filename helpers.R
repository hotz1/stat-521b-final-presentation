# load packages
library(knitr)
library(kableExtra)
library(tidyverse)
library(here)

# Define an "unknown" function to show off the GP regression model
f_tilde <- function(x){
  # Compute function
  fun_vals <-  0.125 * sin(pi*x)^2 * sqrt(x + 4) - sin(pi*x/2) * cos(pi*x^2/7)
  
  # Replace f(x) for any values outside of the range [0, 10] with NA
  fun_vals[which(!between(x, 0, 10))] <- NA
  return(fun_vals)
}

# Plot the function 
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
ggsave(here("images/example_function.jpg"))

# squared exponential kernel function
sq_exp_kernel <- function(X1, X2, len = 1){
  # Squared exponential kernel distance between X1 and X2
  # len is length parameter (scales the distances)
  # noise is the standard deviation of iid noise (scales variances)
  exps <- outer(X1, X2, \(x1, x2) exp(-(x1-x2)^2 / (2 * len^2)))
  return(exps)
}

# Mean and covariance matrix of posterior predictive dist. from GP regression 
posterior_pred <- function(new_x, datapoints, noise_sd = 1, len = 1){
  # new_x contains value(s) of x which we want to predict
  # datapoints is a collection of observed (x,y) pairs
  # kernel is a kernel function
  # noise_sd is the standard deviation of iid noise for observing y = f(x) + e
  # len is length parameter to pass to the kernel
  
  # Get observed data from the matrix
  obs_x <- datapoints[,1]
  obs_y <- datapoints[,2]
  
  # kernels
  obs_to_obs <- sq_exp_kernel(obs_x, obs_x, len)
  obs_to_new <- sq_exp_kernel(obs_x, new_x, len)
  new_to_new <- sq_exp_kernel(new_x, new_x, len)
  
  # Posterior predictive mean and variance
  post_mean <- t(obs_to_new) %*% solve(obs_to_obs + noise_sd^2 * diag(nrow(obsmat))) %*% obs_y
  post_covmat <- new_to_new - t(obs_to_new) %*% solve(obs_to_obs) %*% obs_to_new
  return(list(mean = post_mean, cov = post_covmat))
}

# s