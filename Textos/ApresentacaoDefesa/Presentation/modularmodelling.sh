#modularmodelling.sh
xelatex --src -interaction=nonstopmode modularmodelling
#bibtex modularmodelling
#latex --src -interaction=nonstopmode modularmodelling
#latex --src -interaction=nonstopmode modularmodelling
xelatex --src -interaction=nonstopmode -synctex=1 modularmodelling.tex
# rm -f *.aux *.log *.blg *.dvi *.ps *.toc *.lot *.lof *.idx *.nlo *.nls *.ilg *.out
