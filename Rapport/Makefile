FILE=RapportTER

.PHONY: all clean $(FILE).pdf

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	rubber --pdf $(FILE);

clean:
	@rm -rf *.aux *.blg *.out *.bbl *.log *.toc
