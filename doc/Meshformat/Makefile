# "make" builds all

SRC=Meshformat.tex

all: meshformat
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

meshformat:
	pdflatex ${SRC}
	bibtex ${SRC}||true
	pdflatex ${SRC}
	pdflatex ${SRC}

.PHONY: userguide

