all: 

github:
	./knit github
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape_github.md

docx:
	./knit docx
	pandoc -s -S --biblio wrightscape.bib wrightscape.md --reference-docx=template.docx -o wrightscape.docx


pdf:
	./knit pdf
	pandoc -s -S --biblio wrightscape.bib wrightscape.md --template=elsarticle.latex -o wrightscape.tex
	xelatex wrightscape.tex
	xelatex wrightscape.tex
	xelatex wrightscape.tex

epub:
	./knit epub
	pandoc -s -S --biblio wrightscape.bib wrightscape.md -o wrightscape.epub

bib:
	bibtool -s -d --preserve.key.case=on -x citations.aux ~/Documents/Mendeley/bib/*.bib > wrightscape.bib

clean:
	rm *.tex *.spl *.out *.log *.aux *.docx 
