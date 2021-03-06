#!/bin/bash
# Compile and run C++ code
g++ halfsin.cpp
./a.out
# Open gnuplot and print graphs
gnuplot<<EOF
set terminal jpeg small font arial size 270,205
unset label
set yrange [-0.25:1.25]
set xlabel "Node"
set ylabel "Field"
set output "fft0.jpeg"
plot "A" with lines
set output "fft1.jpeg"
set title " FT: Velocity Field at time t=1"
plot "B" with lines
set output "fft2.jpeg" 
set title " FT: Velocity Field at time t=2"
plot "C" with lines
set output "fft3.jpeg"
set title " FT: Velocity Field at time t=3"
plot "D" with lines	
set output "fft4.jpeg"
set title " FT: Velocity Field at time t=4"
plot "E" with lines	
set output "fft5.jpeg"
set title " FT: Velocity Field at time t=5"
plot "F" with lines
set output "fft6.jpeg"
set title " FT: Velocity Field at time t=6"
plot "G" with lines
set output "fft7.jpeg"
set title " FT: Velocity Field at time t=7"
plot "H" with lines
set output "fft8.jpeg"
set title " FT: Velocity Field at time t=8"
plot "I" with lines
set output "fft9.jpeg"
set title " FT: Velocity Field at time t=9"
plot "J" with lines
set output "fft10.jpeg"
set title "After One Time Period"
plot "K" with lines
set output "fft11.jpeg"
set title " FT: Velocity Field at time t=11"
plot "L" with lines
set output "fft12.jpeg"
set title " FT: Velocity Field at time t=11"
plot "M" with lines
set output "fft13.jpeg" 
set title " FT: Velocity Field at time t=12"
plot "N" with lines
set output "fft14.jpeg"
set title " FT: Velocity Field at time t=13"
plot "O" with lines
EOF
# Create the report pdf and open it
pdflatex report
evince report.pdf
