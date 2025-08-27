library(bregr)
library(ggplot2)

lung <- survival::lung |>
  dplyr::filter(ph.ecog != 3)
lung$ph.ecog <- factor(lung$ph.ecog)

mds <- br_pipeline(
  lung,
  y = c("time", "status"),
  x = colnames(lung)[6:10],
  x2 = c("age", "sex"),
  method = "coxph"
)

p1 = br_show_forest(mds, tab_headers = c("Focal", "Variable", "Level", "N", "", "Hazard Ratio (95% CI)", "P"))
p1
ggsave("fig1_comb_forest_plot.pdf", plot = p1, width = 7, height = 6)

mds2 <- br_pipeline(
  survival::lung,
  y = c("time", "status"),
  x = colnames(survival::lung)[6:10],
  x2 = c("age", "sex"),
  method = "coxph"
)

p2 = br_show_risk_network(mds2)
ggsave("fig1_risk_network_plot.pdf", plot = p2, width = 7, height = 6)

br_show_table(mds)
br_show_table(mds, export = TRUE, args_table_export = list(format = "html"))

# subgroup
data <- survival::lung
data <- data |>
  dplyr::mutate(
    ph.ecog = factor(ph.ecog),
    sex = ifelse(sex == 1, "Male", "Female")
  )
mds <- br_pipeline(
  data,
  y = c("time", "status"),
  x = "ph.ecog",
  group_by = "sex",
  method = "coxph"
)

p3 = br_show_forest(
  mds,
  drop = 2,
  subset = !(Group_variable == "Female" & variable == "ph.ecog" & label == 3),
  tab_headers = c("Focal", "Variable", "Level", "N", "", "Hazard Ratio (95% CI)", "P"),
  #log_first = TRUE
  x_trans = "log"
)
p3
ggsave("fig1_subgroup_forest_plot.pdf", plot = p3, width = 7, height = 6)

# linear
m <- br_pipeline(mtcars,
                 y = "mpg",
                 x = colnames(mtcars)[2:4],
                 x2 = "vs",
                 method = "gaussian"
)
p4 = br_show_forest_ggstats(m)
p4

ggsave("fig1_mod_comparison.pdf", plot = p4, width = 7, height = 4)


pdf("fig1_fitted_line.pdf", width = 7, height = 4)
br_show_fitted_line(m, xvar = "cyl")
dev.off()

pdf("fig1_fitted_line_2d.pdf", width = 7, height = 4)
br_show_fitted_line_2d(m, xvar = "cyl", yvar = "mpg")
dev.off()
