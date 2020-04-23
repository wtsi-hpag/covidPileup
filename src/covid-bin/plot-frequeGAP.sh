#!/bin/tcsh

#cleanFile cover-65.dat2 > bE378K21-raw.dat.cleaned

function plotcmd
{
	printf "set terminal svg\n"
        printf "set logscale x\n"
        printf "set logscale y\n"
	printf "set style line 1 lt 1 lw 2 pt 2 linecolor rgb \"red\"\n"
	printf "set style line 2 lt 1 lw 2 pt 2 linecolor rgb \"blue\"\n"
	printf "set style line 3 lt 1 lw 2 pt 2 linecolor rgb \"green\"\n"
	printf "set xlabel \"Number of GAPs \"\n"
	printf "set ylabel \" SNP Frequency / 100 \"\n"
	printf "plot [ 1 to 10000 ] [ 1 to 100 ] \"public-frequeGAP.freq\" title \"COVID-19 Public: 369 genomes\" with lines ls 1, \"sample-frequeGAP.freq\" title \"COVID-19 Sample: 1439 genomes\" with lines ls 2"
}

plotcmd | gnuplot > data.svg
inkscape -z --export-text-to-path --export-pdf data.pdf data.svg
gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=data.png data.pdf

rm -f bE378K21-raw.dat.cleaned bE378K21-screen.dat.cleaned data.svg
