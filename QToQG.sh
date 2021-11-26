
Process="QToGQ"
MakeFile="PlotMaking/MakefileQToGQ"

P="300" # 300xTemperature
z="0.25"

mkdir -p PlotMaking/$Process

make -f $MakeFile OpacityRate
./OpacityRate.exe -P $P -z $z
mkdir -p PlotMaking/$Process/Opacity
mv -f Opacity/Rate-P300-z0.25.txt PlotMaking/$Process/Opacity/

make -f $MakeFile ImprovedOpacity
./ImprovedOpacity.exe -P $P -z $z
mkdir -p PlotMaking/$Process/OpacityImproved
mv -f OpacityImproved/Rate-P300-z0.25.txt PlotMaking/$Process/OpacityImproved/

make -f $MakeFile Harmonic
./HO.exe -P $P -z $z
mkdir -p PlotMaking/$Process/HO
mv -f HO/Rate-P300-z0.25.txt PlotMaking/$Process/HO/

make -f $MakeFile FullRate
./FullRate.exe -P $P -z $z
mkdir -p PlotMaking/$Process/FullRate
mv -f OUTPUT/Rate-P300-z0.25.txt PlotMaking/$Process/FullRate/