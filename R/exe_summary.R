execute_res <- arrow::read_csv_arrow(paste0(
  here::here(),
  "/data/execute_res.csv"
))

execute_res |>
  dplyr::group_by(company, model) |>
  dplyr::summarise(execution_ratio = sum(ok) / length(ok)) |>
  dplyr::arrange(dplyr::desc(execution_ratio)) |>
  tinytable::tt()
