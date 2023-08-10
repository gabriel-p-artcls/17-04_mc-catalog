#!/bin/bash
#
# Script to generate a PDF with all differences marked
# via the latexdiff package, between two .tex files.
# See: http://www.ctan.org/pkg/latexdiff
#
# All output is stored in a diff/ sub-folder that is
# created if it doesn't exist.
#
# Run with:
#
# ./latexdiff.sh old.tex new.tex [-s]
#
# The -s flag is optional and can be added at the end to
# turn on silent mode (no output to screen).
#
# Pausing:
# read -p "Press [Enter] key to resume..."

clear

# Check if silent flag is on.    
if [[ "$3" == "-s" ]]; then
	echo "Silent mode ON."
    outvar="/dev/null"
else
    outvar="/dev/stdout"
fi

# Check if diffs/ dir exists.
DIRECTORY="diffs/"
if [ ! -d "$DIRECTORY" ]; then
	mkdir "$DIRECTORY"
  	echo "diff/ dir created."
fi

# Copy necessary .bib, .bst, .cls and .sty files to /diff
# folder.
cp *.bib diffs/ 2> /dev/null
cp *.bst diffs/ 2> /dev/null
cp *.cls diffs/ 2> /dev/null
cp *.sty diffs/ 2> /dev/null

# launch latexdiff script.
echo "Running latexdiff"
latexdiff $1 $2 > diffs/diff.tex
echo "latexdiff done"

# Enter /diffs folder.
cd diffs

echo "Compiling"

file_t="diff.tex"
file_a="diff.aux"
# Run LaTeX commands sequentially.
pdflatex -draftmode "$file_t" > $outvar
echo "pdflatex"
bibtex "$file_a" > $outvar
echo "bibtex"
pdflatex -draftmode "$file_t" > $outvar
echo "pdflatex"
pdflatex "$file_t" > $outvar

# Open PDF with Okular.
echo "Opening PDF."
setsid okular diff.pdf 2> /dev/null &

echo "Finished."