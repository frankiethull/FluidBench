# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   ~~~~ FluidBench ~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# each question was asked once,
# in series, 01-05,
# one-shot correct or failed,
# benchmark will be three-parter:
# 1) does the code execute?
#  --- an easy first-pass at least executable code
# 2) 'solved problem'
#  --- is the output as expected ?
# 3) 'correctness'
#  --- eval code + maths of simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

llm_answers <- list.files(
  path = "bench",
  full.names = TRUE,
  recursive = TRUE,
  pattern = ".R"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### round one: execution

executed_df <- tibble::tibble(path = llm_answers) |>
  tidyr::separate(
    path,
    into = c("bench", "company", "model", "file"),
    sep = "/",
    remove = FALSE
  ) |>
  dplyr::mutate(
    script = file,
    result = purrr::map(path, purrr::safely(~ source(.x, local = TRUE)))
  )

execute_res <- executed_df |>
  dplyr::mutate(
    ok = map_lgl(result, ~ is.null(.x$error)),
    error = map_chr(
      result,
      ~ ifelse(is.null(.x$error), NA_character_, as.character(.x$error))
    )
  ) |>
  dplyr::select(company, model, script, path, ok, error)

arrow::write_csv_arrow(
  execute_res,
  paste0(here::here(), "/data/execute_res.csv")
)

#### round two: solution

#### round three: accuracy

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
