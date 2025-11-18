#' Suppress terminal output
#'
#'
#' @examples
#'
#' # use with tidyverse, for instance
#' group_by(...) %>% silence
#' # to suppress terminal messages
#'

silence <- function(...){
  capture.output({
    suppressPackageStartupMessages({
      suppressMessages({
        out <- suppressWarnings(...)
      })
    })
  }, type = "output")

  return(out)
}


#' Clear global environment of variables, functions and plots as well as console
#'
#' @param ...        list of variables/functions to keep
#'
#' @examples
#'
#' f <- function(x,y) x+y
#' clean_up(f) # f is still available
#'

clean_up <- function(...){
  # Aux function to clear work environment

  ## Clear console and plots, and remove variables from workspace unless given as arguments
  for(i in dev.list()) dev.off(which=i)
  keep = as.list(match.call())
  rm(list=setdiff(ls(envir = .GlobalEnv), keep), envir = .GlobalEnv)
  cat("\014")  # Clear console
}


#' Add transparency in plot colors
#'
#' @param col          specify a color, or vector of colors
#' @param alpha        range in $[0,1]$, 0 is complete transparency, 1 is no transparency
#'
#' @return color(s) with added transparency in hex code
#'
#' @examples
#'
#' add_alpha("red",.5)

add_alpha <- function(col, alpha=1){
  ## Add an alpha value to a color
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

