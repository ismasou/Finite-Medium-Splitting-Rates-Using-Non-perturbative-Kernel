

TofmOverC = 0.394
TofmOverC2 = 0.716
TimeConversion(P,z) = 2*P*z*(1.-z)

set term epslatex color font ",14in" size 4.8,4 standalone lw 3 ps 2 header '\newcommand{\hl}[1]{\setlength{\fboxsep}{0.75pt}\colorbox{white}{#1}}'
set style line 1 lc rgb '#110141';   set style line 2 lc rgb '#710162';   set style line 11 lc rgb '#a12a5e';
set style line 4 lc rgb '#ed0345';   set style line 5 lc rgb '#ef6a32';   set style line 12 lc rgb '#fbbf45';
set style line 7 lc rgb '#aad962';   set style line 6 lc rgb '#03c383';   set style line 9 lc rgb '#017351';
set style line 10 lc rgb '#01545a';   set style line 3 lc rgb '#26294a';   set style line 8 lc rgb '#1a1334';
set style line 13 lc rgb '#58587a';
set style line 14 lc rgb '#9B4D91';
set style line 15 lc rgb '#BF5428';
set style line 16 lc rgb '#03c383';

lsLO  = 5
lsNLO = 16
lsEQCD = 14

lsLO2 = 15
lsNLO2 = 9
lsEQCD2= 2

OutputFold= "Temp"
File = OutputFold."/Rate"
system("mkdir ".OutputFold)


gs2 = 2.763515753;
gSqr250 = 3.725027366;
CA = 3.; CF= 4./3.;
Pgg(z) = 2.*CA*(1.-z*(1.-z))**2/(z*(1.-z));
Pqq(z) = 0.5*(z**2+(1.-z)**2);
Pqg(z) = CF*(1+(1.-z)**2)/z;

FONTSIZE = '\large '
set log
set key at graph 0.99, graph 0.45 width 1 height 3 font ',30'
set xlabel FONTSIZE.' Evolution time: $t [fm/c]$ '

set xrange[0.12:10]

Theta(x,min) = ( x > min ? 1:1/0)
Theta1(x,max) = 1 #( TimeConversion(300.,0.1)*x < max ? 1:1/0)


P = 300.0
z = 0.25
set output File.".tex"

set ylabel FONTSIZE.' Splitting rate: $\frac{1}{g^{4}T}\frac{d\Gamma (P,z)}{dz}$' offset 3,0
set format y FONTSIZE.' $10^{%T}$'
set yrange[2e-5:1]
set log y
set ytics autofreq
set label 1 FONTSIZE.' $q\to gq$' at graph 0.05, graph 0.95
set label 2 gprintf(FONTSIZE.' $z=%g$',z) at graph 0.05, graph 0.85
EQCDDat   = system(sprintf(" awk '( ($2/%g-1)*($2/%g-1) < 0.002 && ($1/%g-1)*($1/%g-1) < 0.001 ){ print $4}' ../InfiniteLData/EQCDFullRate-T500.txt ",z,z,P,P))

EQCD(x) = EQCDDat/gs2**2*Theta(x,1)*Pqg(z)

pl "FullRate/Rate-P300-z".gprintf("%g",z).".txt" u (TofmOverC*TimeConversion(P,z)*$1):($2/gs2**2*Pqg(z)) w l ls lsEQCD lw 3 ti FONTSIZE.' $T=500$MeV',\
	"Opacity/Rate-P300-z".gprintf("%g",z).".txt"   u (TofmOverC*TimeConversion(P,z)*$1):($2/gs2**2*Theta1($1,100.5)*Pqg(z)) w l ls lsLO lw 3 dt 3 ti FONTSIZE.' Opacity N=1',\
	"OpacityImproved/Rate-P300-z".gprintf("%g",z).".txt"   u (TofmOverC*TimeConversion(P,z)*$1):($2/gs2**2*Theta1($1,100.5)*Pqg(z)) w l ls lsNLO lw 3 dt 3 ti FONTSIZE.' Opacity N=x',\
	"HO/Rate-P300-z".gprintf("%g",z).".txt"   u ( TofmOverC*TimeConversion(P,z)*$1):(($2)/gs2**2*Theta1($1,100.5)*Pqg(z)) smo cs ls 4 lw 3 dt 4 ti FONTSIZE.' NLO-HO',\
	EQCD(x) w l lw 3 lc rgb "#aaaaaa" dt "-" ti FONTSIZE.' AMY'

set output

sys("latexmk -f -ps -jobname=".File." ".File.".tex " )
sys("ps2pdf ".File.".ps ".File.".pdf")
sys("mv ".File.".pdf .")

sys("rm -rf ".OutputFold)
sys("clear")
