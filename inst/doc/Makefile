all: 

github:
	./knit github
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape_github.md

md:
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape.md

docx:
	./knit docx
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape.docx


pdf:
	./knit pdf
	pandoc -s -S --biblio wrightscape.bib wrightscape.md --template=elsarticle.latex --latex-engine=xelatex -o wrightscape.tex
	xelatex wrightscape.tex
	xelatex wrightscape.tex
	xelatex wrightscape.tex

epub:
	./knit epub
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape.epub

bib:
	bibtool -s -d -x citations.aux ~/Documents/Mendeley/bib/*.bib > wrightscape.bib


