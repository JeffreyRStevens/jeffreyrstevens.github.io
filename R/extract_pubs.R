
cvfile <- readLines("cv.qmd")
startline <- which(cvfile == "## Publications")
endline <- which(cvfile == "## Presentations")
body <- cvfile[(startline+1):(endline-2)]

prelim <- c("---", 'title: "Publications"', 'format:', '  html:', "    css: flipmargin.css", "    toc: true", "    toc-location: left", "---")

newfile <- c(prelim, body)
writeLines(newfile, "publications.qmd")
