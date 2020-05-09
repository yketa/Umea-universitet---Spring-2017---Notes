all: main.pdf

%.pdf: dir
	@pdflatex $*
	@pdflatex $*
	-@biber $* || bibtex $*
	@pdflatex $*
	@pdflatex $*
	@mv `ls $*.* | grep -v '\.tex\|\.pdf'` aux

dir:
	@mkdir -p aux

clean:
	@rm -f `ls main.* | grep -v '\.tex\|\.pdf'`
	@rm -rf aux

mrproper: clean
	@rm -f *converted-to.pdf
	@rm -f main.pdf

