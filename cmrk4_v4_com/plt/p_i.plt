Nx=  ` sed -n "/#define Nx / p" ../rps.h | cut -f3 -d" " `
Ny=  ` sed -n "/#define Ny / p" ../rps.h | cut -f3 -d" " `
NF=  ` sed -n "/#define NF / p" ../rps.h | cut -f3 -d" " `

set terminal png size 500+55, 500+55
set xtics offset 0.0, 0.2
set ytics offset 0.0, 0.2 
set xrange [1-0.5:Nx+0.5]
set yrange [1-0.5:Ny+0.5]
set size ratio -1
unset key
set palette defined ( 0 0 0 0, 1 1 1 1 )

i=0
while (i <= 99){
	set output sprintf("p_i-%d.png", i)
	plot sprintf("../dat/p_i-%d.dat", i) u ($1+1):($2+1):($3) matrix w image
	unset output
	i=i+1
}
