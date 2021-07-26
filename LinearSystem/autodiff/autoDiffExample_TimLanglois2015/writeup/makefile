all: mls.pdf

# builds the pdf
mls.pdf: mls.tex
	latexmk -pdf mls.tex

# removes all intermediate files
clean:
	latexmk -c

# purge will also remove the output file (.pdf)
purge:
	latexmk -C
