#' @importFrom utils menu
#' @importFrom tcltk .TkUp tk_select.list

select.list.edit<-function (choices, preselect = NULL, multiple = TRUE, title = NULL, graphics = getOption("menu.graphics")) {
  #FLAG 
  C_selectlist<-NULL
  if (!interactive()) 
    stop("select.list() cannot be used non-interactively")
  if (!is.null(title) && (!is.character(title) || length(title) != 
                          1)) 
    stop("'title' must be NULL or a length-1 character vector")
  if (isTRUE(graphics)) {
    if (.Platform$OS.type == "windows" || .Platform$GUI == 
        "AQUA") 
      return(.External2(C_selectlist, choices, preselect, 
                        multiple, title))
    else if (graphics && capabilities("tcltk") && capabilities("X11") && 
             suppressWarnings(tcltk::.TkUp)) 
      return(tcltk::tk_select.list(choices, preselect, 
                                   multiple, title))
  }
  if (!multiple) {
    res <- utils::menu(choices, FALSE, title)
    if (res < 1L || res > length(choices)) 
      return("")
    else return(choices[res])
  }
  else {
    nc <- length(choices)
    if (length(title) && nzchar(title[1L])) 
      cat(title, "\n", sep = "")
    def <- if (is.null(preselect)) 
      rep(FALSE, nc)
    else choices %in% preselect
    op <- paste0(format(seq_len(nc)), ": ", ifelse(def, "+", 
                                                   " "), " ", choices)
    if (nc > 10L) {
      fop <- format(op)
      nw <- nchar(fop[1L], "w") + 2L
      ncol <- getOption("width")%/%nw
      if (ncol > 1L) 
        op <- paste0(fop, c(rep("  ", ncol - 1L), "\n"), 
                     collapse = "")
      cat("", op, sep = "\n")
    }
    else cat("", op, "", sep = "\n")
    cat(gettext("Enter one or more numbers separated by spaces, or an empty line to cancel\n"))
    repeat {
      res <- tryCatch(scan("", what = 0, quiet = TRUE, 
                           nlines = 1), error = identity)
      if (!inherits(res, "error")) 
        break
      cat(gettext("Invalid input, please try again\n"))
    }
    if (!length(res) || (length(res) == 1L && !res[1L])) 
      return(character())
    #res <- sort(res[1 <= res && res <= nc])
    return(choices[res])
  }
}
