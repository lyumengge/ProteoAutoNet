##################################################
## Project: Thyroid cell line
## Script purpose: Main figure - SEC trace
## Date: 2024-12-09
## Author: Mengge LYU
## Version: 1.0.0
##################################################

# Library loading ---------------------------------------------------------
library(rio)
library(CCprofiler)
library(data.table)
library(dplyr)
library(ggplot2)

# Database loading --------------------------------------------------------

## TPC ----------------------------------------------------------------------
setwd("E:\\MenggeLYU\\SECtrace_plot")
TPC1 <- rio::import("Thyroid_cell_line_TPC1.txt")
TPC2 <- rio::import("Thyroid_cell_line_TPC2.txt")
TPC3 <- rio::import("Thyroid_cell_line_TPC3.txt")
TPC1 <- TPC1[,c(1,3)]
TPC2 <- TPC2[,c(1,3)]
TPC3 <- TPC3[,c(1,3)]

TPC1$Replicate <- "Replicate 1"
TPC2$Replicate <- "Replicate 2"
TPC3$Replicate <- "Replicate 3"

combined_data <- bind_rows(
  TPC1 %>% rename(X = 1, Y = 2),
  TPC2 %>% rename(X = 1, Y = 2),
  TPC3 %>% rename(X = 1, Y = 2)
)

ggplot(combined_data, aes(x = X, y = Y)) +
  geom_line(linewidth = 1, color = "black") +  # Set all points to black
  facet_wrap(~Replicate, ncol = 1) +  # Arrange facets in one column
  labs(
    title = "TPC-1",
    x = "Retention time (min)",
    y = "Absorbance (280 nm)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_blank(),  # Hide the replicate labels in the facet
    legend.position = "none"  # Hide the legend
  )


## Nthy --------------------------------------------------------------------


Nthy1 <- rio::import("Thyroid_cell_line_Nthy1.txt")
Nthy2 <- rio::import("Thyroid_cell_line_Nthy2.txt")
Nthy3 <- rio::import("Thyroid_cell_line_Nthy3.txt")
Nthy1 <- Nthy1[,c(1,3)]
Nthy2 <- Nthy2[,c(1,3)]
Nthy3 <- Nthy3[,c(1,3)]

Nthy1$Replicate <- "Replicate 1"
Nthy2$Replicate <- "Replicate 2"
Nthy3$Replicate <- "Replicate 3"

combined_data <- bind_rows(
  Nthy1 %>% rename(X = 1, Y = 2),
  Nthy2 %>% rename(X = 1, Y = 2),
  Nthy3 %>% rename(X = 1, Y = 2)
)

ggplot(combined_data, aes(x = X, y = Y)) +
  geom_line(linewidth = 1, color = "black") +  # Set all points to black
  facet_wrap(~Replicate, ncol = 1) +  # Arrange facets in one column
  labs(
    title = "Nthy",
    x = "Retention time (min)",
    y = "Absorbance (280 nm)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_blank(),  # Hide the replicate labels in the facet
    legend.position = "none"  # Hide the legend
  )

## FTC238 ------------------------------------------------------------------
FTC2381 <- rio::import("Thyroid_cell_line_FTC2381.txt")
FTC2382 <- rio::import("Thyroid_cell_line_FTC2382.txt")
FTC2383 <- rio::import("Thyroid_cell_line_FTC2381.txt")
FTC2381 <- FTC2381[,c(1,3)]
FTC2382 <- FTC2382[,c(1,3)]
FTC2383 <- FTC2383[,c(1,3)]

FTC2381$Replicate <- "Replicate 1"
FTC2382$Replicate <- "Replicate 2"
FTC2383$Replicate <- "Replicate 3"

combined_data <- bind_rows(
  FTC2381 %>% rename(X = 1, Y = 2),
  FTC2382 %>% rename(X = 1, Y = 2),
  FTC2383 %>% rename(X = 1, Y = 2)
)

# Plot using ggplot2
ggplot(combined_data, aes(x = X, y = Y)) +
  geom_line(linewidth = 1, color = "black") +  # Set all points to black
  facet_wrap(~Replicate, ncol = 1) +  # Arrange facets in one column
  labs(
    title = "FTC238",
    x = "Retention time (min)",
    y = "Absorbance (280 nm)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_blank(),  # Hide the replicate labels in the facet
    legend.position = "none"  # Hide the legend
  )


## Standard ----------------------------------------------------------------
FTC2381 <- rio::import("Thyroid_cell_line_SP1.txt")
FTC2382 <- rio::import("Thyroid_cell_line_SP2.txt")
FTC2383 <- rio::import("Thyroid_cell_line_SP3.txt")
FTC2381 <- FTC2381[,c(1,3)]
FTC2382 <- FTC2382[,c(1,3)]
FTC2383 <- FTC2383[,c(1,3)]

FTC2381$Replicate <- "Replicate 1"
FTC2382$Replicate <- "Replicate 2"
FTC2383$Replicate <- "Replicate 3"

combined_data <- bind_rows(
  FTC2381 %>% rename(X = 1, Y = 2),
  FTC2382 %>% rename(X = 1, Y = 2),
  FTC2383 %>% rename(X = 1, Y = 2)
)

ggplot(combined_data, aes(x = X, y = Y)) +
  geom_line(linewidth = 1, color = "black") +  # Set all points to black
  facet_wrap(~Replicate, ncol = 1) +  # Arrange facets in one column
  labs(
    title = "Standard Mix",
    x = "Retention time (min)",
    y = "Absorbance (280 nm)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_blank(),  # Hide the replicate labels in the facet
    legend.position = "none"  # Hide the legend
  )



