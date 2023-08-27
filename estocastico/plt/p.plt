set terminal png size 1000, 1000
ext="png"
unset xtics
unset ytics
unset colorbox

unset key
set xrange [0:1]
set yrange [0:1]
set size ratio -1

set palette defined ( 1 "#ff0000",\
                      2 "#0000ff",\
                      3 "#ffff00" )

i=0
while (i <= 10 ){
	set output sprintf("p-%d.%s", i, ext)
	plot sprintf("../dat/p-%d.dat", i) u 1:2:3 pt 7 lc palette,\
	i= i+1
}
unset output
