# Cold-start launcher: shiny::runApp() needs shiny pre-installed, but app.R's
# auto-installer (which handles every other dependency) only runs once the app
# is sourced. Bootstrap shiny here, then launch. LAUNCH_BROWSER=false for
# headless / agent / CI use.
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny", repos = "https://cloud.r-project.org")
}
open_browser <- !identical(tolower(Sys.getenv("LAUNCH_BROWSER")), "false")
shiny::runApp(getwd(), launch.browser = open_browser)
