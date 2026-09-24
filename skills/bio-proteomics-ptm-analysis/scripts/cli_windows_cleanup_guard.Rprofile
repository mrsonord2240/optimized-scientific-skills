# Opt-in audit/batch workaround for the Windows cli 3.6.6 teardown crash.
# Use only through R_PROFILE_USER for one disposable process; do not install globally.
mark_cli_unload_guarded <- function(...) {
  if (.Platform$OS.type != "windows" || !("cli" %in% loadedNamespaces())) return(invisible(NULL))
  cli_namespace <- asNamespace("cli")
  if (!exists("clienv", envir = cli_namespace, inherits = FALSE)) return(invisible(NULL))
  cli_state <- get("clienv", envir = cli_namespace, inherits = FALSE)
  if (is.environment(cli_state) && identical(cli_state$unloaded, FALSE)) cli_state$unloaded <- TRUE
  receipt <- Sys.getenv("CLI_WINDOWS_CLEANUP_GUARD_RECEIPT", unset = "")
  if (nzchar(receipt) && is.environment(cli_state) && identical(cli_state$unloaded, TRUE)) {
    dir.create(dirname(receipt), recursive = TRUE, showWarnings = FALSE)
    connection <- file(receipt, open = "wx")
    on.exit(close(connection), add = TRUE)
    writeLines(c("cli_unloaded=TRUE", paste0("pid=", Sys.getpid())), connection, useBytes = TRUE)
  }
  invisible(NULL)
}
setHook(packageEvent("cli", "onLoad"), mark_cli_unload_guarded, action = "append")
