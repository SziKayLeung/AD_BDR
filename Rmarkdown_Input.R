pilltabs <- function(tab, count = TRUE, rows = TRUE, cols = TRUE, chisq = TRUE, resid = TRUE, row.names = TRUE) {
  
  if (!requireNamespace("questionr", quietly = TRUE))
    stop("the questionr package is needed for the pilltabs() function to work. Please install it.",
         call. = FALSE)  
  
  res <- list()
  
  if (count) res[["count"]] <- kable(tab, output = FALSE, row.names = row.names)
  if (rows)  res[["rows"]] <- kable(round(questionr::rprop(tab, n = TRUE),1), output = FALSE, row.names = row.names)
  if (cols)  res[["cols"]] <- kable(round(questionr::cprop(tab, n = TRUE),1), output = FALSE, row.names = row.names)
  if (resid) res[["resid"]] <- kable(round(questionr::chisq.residuals(tab),2), output = FALSE, row.names = row.names)
  if (chisq) {
    test <- stats::chisq.test(tab)
    res[["chisq"]] <- paste0('X-squared = ', round(test$statistic, 4), 
                             ', df = ', test$parameter,
                             ', p = ', format.pval(test$p.value, digits = 4)) 
  }
  class(res) <- "pilltabs"
  res
}

# plot theme
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.80, 0.80),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="CM Roman"),
                 axis.text.x= element_text(size=12,family="CM Roman"),
                 axis.text.y= element_text(size=12,family="CM Roman"))