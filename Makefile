all:
	@pdflatex main
	@pdflatex main
	@bibtex main
	@pdflatex main
	@pdflatex main
	@mkdir -p aux
	@mv `ls main.* | grep -v tex | grep -v pdf` aux
clean:
	@rm -f `ls main.* | grep -v tex | grep -v pdf`
	@rm -rf aux
