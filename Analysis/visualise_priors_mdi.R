

library(ggplot2)
library(MASS)

# Generate gamma rvs

x <- rgamma(1000000, shape = 1, rate = 0.2)

den <- density(x)

dat <- data.frame(x = den$x, y = den$y)

# Plot density as points

phi_prior_plt <- ggplot(data = dat, aes(x = x, y = y)) + 
  geom_line() +
  theme_bw() +
  labs(
    title = bquote("Prior distribution of" ~ phi ~ "parameter for MDI"),
    x = bquote(phi[ij]),
    y = "Density"
  )

ggsave("~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Data_inspection/phi_prior.png", phi_prior_plt)
  
# geom_point(size = 3) +
  # theme_classic()

x <- rgamma(1000000, shape = 2, rate = 4)

den <- density(x)

dat <- data.frame(x = den$x, y = den$y)


alpha_prior_plt <- ggplot(data = dat, aes(x = x, y = y)) + 
  geom_line() +
  theme_bw() +
  labs(
    title = "Prior distribution of mass parameter for MDI",
    x = bquote(alpha[l]),
    y = "Density"
  )

ggsave("~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Data_inspection/alpha_prior.png", alpha_prior_plt)
