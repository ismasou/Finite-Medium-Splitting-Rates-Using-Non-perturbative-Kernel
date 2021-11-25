
Process="GToGG"
MakeFile="PlotMaking/MakefileGToGG"

P="300" # 300xTemperature
z="0.25"

mkdir PlotMaking/$Process

make -f $MakeFile OpacityRate
./OpacityRate.out -P $P -z $z
mkdir PlotMaking/$Process/Opacity
mv -f Opacity/Rate-P300-z0.25.txt PlotMaking/$Process/Opacity/

make -f $MakeFile ImprovedOpacity
./ImprovedOpacity.out -P $P -z $z
mkdir PlotMaking/$Process/OpacityImproved
mv -f OpacityImproved/Rate-P300-z0.25.txt PlotMaking/$Process/OpacityImproved/

make -f $MakeFile Harmonic
./HO.out -P $P -z $z
mkdir PlotMaking/$Process/HO
mv -f HO/Rate-P300-z0.25.txt PlotMaking/$Process/HO/

make -f $MakeFile FullRate
./FullRate.out -P $P -z $z
mkdir PlotMaking/$Process/FullRate
mv -f OUTPUT/Rate-P300-z0.25.txt PlotMaking/$Process/FullRate/
 