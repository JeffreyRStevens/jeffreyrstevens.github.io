library(tidyverse)
library(here)
library(patchwork)

set.seed(123)
small_large <- data.frame(x = 1:5, y = 1:5) |>
  complete(x, y) |>
  mutate(gr = "(a)")
small_small <-data.frame(x = 1:5, y = 1:5) |>
  complete(x, y) |>
  mutate(gr = "(b)") |>
  slice_sample(n = 24)
small_samples <- bind_rows(small_large, small_small)

large_large <- data.frame(x = 1:8, y = 1:8) |>
  complete(x, y) |>
  mutate(gr = "(a)") |>
  slice_sample(n = 50)
large_small <-data.frame(x = 1:8, y = 1:8) |>
  complete(x, y) |>
  mutate(gr = "(b)") |>
  slice_sample(n = 48)
large_samples <- bind_rows(large_large, large_small)

small <- small_samples |>
  ggplot(aes(x = x, y = y)) +
  geom_jitter() +
  facet_wrap(vars(gr)) +
  coord_cartesian(xlim = c(0,6), ylim = c(0,6)) +
  theme_void() +
  theme(aspect.ratio = 1,
        plot.background = element_rect(fill = "white",
                                       color = NA))
ggsave(here("posts/dog_number/images/weber_small.png"), height = 3, width = 6)

large <- large_samples |>
  ggplot(aes(x = x, y = y)) +
  geom_jitter() +
  facet_wrap(vars(gr)) +
  coord_cartesian(xlim = c(0,9), ylim = c(0,9)) +
  theme_void() +
  theme(aspect.ratio = 1,
        plot.background = element_rect(fill = "white",
                                       color = NA))
ggsave(here("posts/dog_number/images/weber_large.png"), height = 3, width = 6)

small / large
ggsave(here("posts/dog_number/images/weber.png"), height = 6, width = 6)
