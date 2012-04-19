LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	*.map\
	*.sh\
	sequences.tex\
        *.out

vignette:
	$(R_HOME)/bin/R CMD Sweave ddgraph.Rnw
	pdflatex ddgraph.tex
	bibtex ddgraph
	pdflatex ddgraph.tex
	pdflatex ddgraph.tex
	rm -f Rplots.pdf
clean:
	rm -f $(LATEXFILES)
	rm -f *~

