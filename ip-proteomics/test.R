library(samr)
library(ROTS)
library(limma)
library(DEP)
library(tidyverse)

View(diamonds)

d <- diamonds[sample(nrow(diamonds), 2000), ]


ggplot(data = d, aes(
        x = carat,
        y = price,
        col = carat,
        size = carat
)) +
        geom_point()


ggplot(data = d, aes(
        x = carat,
        y = price,
        col = cut,
        size = carat
)) +
        geom_point()

ggplot(data = d, aes(
        x = carat,
        y = price,
        col = cut
)) +
        geom_point()




