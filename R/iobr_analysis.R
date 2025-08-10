# iobr signature 到luad中分析，然后挑多个画连接图，和单一多变量图，再看泛瘤种某个sigture的影响差异，用一个附图放一些可拓展性的图表

#remotes::install_github("wanglabcsu/bregr")
library(bregr)
library(UCSCXenaShiny)
library(data.table)
library(dplyr)

data_path = "/Users/wsx/Nutstore Files/舒晨阳/回归分析可视化项目/reg_project/data"

dir(data_path)

load(file.path(data_path, "IOBR.all.signature_integration_tpm.TCGA_LUAD.Rdata"))
#load(file.path(data_path, "IOBR.expr.tme_combine.TCGA_LUAD.Rdata"))
#info = fread(file.path(data_path, "TCGA_info.tsv"))

luad_cli = tcga_clinical_fine |>
  filter(Cancer == "LUAD") |>
  inner_join(tcga_surv, by = c("Sample"="sample")) |>
  filter(Code == "TP") |>
  select(Sample, Age, Gender, Stage_ajcc,
         OS, OS.time)

# luad_sig = expr.tme_combine |>
#   mutate(ID = substr(ID, 1, 15),
#          ID = gsub("\\.", "-", ID)) |>
#   distinct(ID, .keep_all = TRUE)

luad_sig = expr.all.signature_integration |>
  mutate(ID = substr(ID, 1, 15),
         ID = gsub("\\.", "-", ID)) |>
  select(-Index) |>
  select(ID, ends_with("ssGSEA")) |>
  distinct(ID, .keep_all = TRUE) |>
  # mutate_at(vars(ends_with("ssGSEA")),
  #           ~ifelse(.>median(., na.rm = TRUE), 1L, 0L)) |>
  rename_at(vars(ends_with("ssGSEA")), ~gsub("_ssGSEA", "", .))

data_luad = inner_join(
  luad_cli, luad_sig,
  by = c("Sample" = "ID")
)

result = br_pipeline(
  data_luad,
  y = c("OS.time", "OS"),
  x = setdiff(colnames(luad_sig), "ID"),
  x2 = c("Age", "Gender", "Stage_ajcc"),
  method = "coxph"
)

#View(result@results_tidy)
br_get_results(result, tidy = TRUE, p.value < 0.05)
br_get_results(result, tidy = TRUE, p.value < 0.05 & Focal_variable == term) |> View()

results_all = br_get_results(result, tidy = TRUE)
write.csv(results_all, file = "LUAD_IOBR_signature_scan_results.csv", quote = FALSE)

result_sig = br_get_results(result, tidy = TRUE, p.value < 0.05 & Focal_variable == term)
top_risk = result_sig |>
  #filter(estimate > 1 & estimate < 200) |>
  filter(estimate > 1) |>
  slice_min(p.value, n = 5) |>
  pull(term)

top_prot = result_sig |>
  #filter(estimate < 1 & estimate > 0.05) |>
  filter(estimate < 1) |>
  slice_min(p.value, n = 5) |>
  pull(term)


br_show_forest(
  result, rm_controls = TRUE, subset = p.value < 0.05
)

#devtools::load_all("../bregr")
#debug(br_show_forest)

p1 = br_show_forest(
  result, rm_controls = TRUE,
  subset = term %in% c(top_risk, top_prot),
  drop = c(1, 3),
  log_first = TRUE
  #ticks_at = c(0.001, 1, 10, 1e8),
)
p1

ggplot2::ggsave(filename = "p1.pdf", plot = p1,
                #dpi = 300,
                width = 9, height = 4, units = "in")

p2 = br_show_forest_ggstats(
  result, idx = "Winter_hypoxia_signature"
)
p2
ggplot2::ggsave(filename = "p2.pdf", plot = p2,
                width = 7, height = 3.5, units = "in")


p3 = br_show_forest_ggstats(
  result, idx = "Nature_metabolism_Hypoxia"
)
p3
ggplot2::ggsave(filename = "p3.pdf", plot = p3,
                width = 7, height = 3.5, units = "in")

p4 = br_show_forest_ggstats(
  result, idx = "Hu_hypoxia_signature"
)
p4
ggplot2::ggsave(filename = "p4.pdf", plot = p4,
                width = 7, height = 3.5, units = "in")

br_show_table(result, term %in% c(top_risk, top_prot),
              export = TRUE, args_table_export = list(format = "html"))


br_show_table_gt(result, idx = c(top_risk, top_prot))

tb1 = br_show_table_gt(result, idx = "Winter_hypoxia_signature") |>
  gtsummary::add_n()

tb1 |>
  gtsummary::as_gt() |>
  gt::gtsave(filename = "p2_table.png")

br_get_results(result)

devtools::load_all("../bregr/")
#undebug(br_show_risk_network)
c(top_risk, top_prot)
p = br_show_risk_network(result,
                         Focal_variable %in% c(top_risk, top_prot)) +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 8),
                 axis.text.x = element_text(angle = 5))
p
ggplot2::ggsave(filename = "p_risk_network.pdf", plot = p,
                width = 9, height = 9, units = "in")
