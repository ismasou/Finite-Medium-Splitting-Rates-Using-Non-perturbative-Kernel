#!/bin/bash
cd GToGG
gnuplot "Plot.gp"
cd .. 
\cp -f GToGG/Rate.pdf RateGToGG.pdf


cd QToGQ
gnuplot "Plot.gp"
cd .. 
\cp -f QToGQ/Rate.pdf RateQToGQ.pdf


cd GToQQ
gnuplot "Plot.gp"
cd .. 
\cp -f GToQQ/Rate.pdf RateGToQQ.pdf