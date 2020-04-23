#!/bin/tcsh

#cleanFile cover-65.dat2 > bE378K21-raw.dat.cleaned

function plotcmd
{
        printf "set logscale y\n"
	printf "set terminal svg\n"
	printf "set style line 1 lt 1 lw 2 pt 2 linecolor rgb \"red\"\n"
	printf "set style line 2 lt 1 lw 2 pt 2 linecolor rgb \"blue\"\n"
	printf "set xlabel \"COVID-19 genome coordinates\"\n"
	printf "set ylabel \"GAP Pileup Frequency / 100 \"\n"
	printf "plot [ 0 to 30000 ] [ 1 to 100 ] \"public-pileupGAP.freq\" title \"COVID-19 Public: 368 genomes\" with lines ls 1, \"sample-pileupGAP.freq\" title \"COVID-19 Sample: 1439 genomes\" with lines ls 2"
}

plotcmd | gnuplot > data.svg
inkscape -z --export-text-to-path --export-pdf data.pdf data.svg
gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=data.png data.pdf

rm -f bE378K21-raw.dat.cleaned bE378K21-screen.dat.cleaned data.svg
