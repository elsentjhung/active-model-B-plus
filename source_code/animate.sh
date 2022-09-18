#!/bin/bash

nt=0
nT=1000000

one="$"1
two="$"2
three="$"3
four="$"4

while [ "$nt" -le "$nT" ];
do
  
B=99999999
if [ "$nt" -le "$B" ];
    then Output=image0"$t".png;
fi

B=9999999
if [ "$nt" -le "$B" ];
    then Output=image00"$t".png;
fi

B=999999
if [ "$nt" -le "$B" ];
    then Output=image000"$t".png;
fi

B=99999
if [ "$nt" -le "$B" ];
    then Output=image0000"$t".png;
fi

B=9999
if [ "$nt" -le "$B" ];
    then Output=image00000"$t".png;
fi

B=999
if [ "$nt" -le "$B" ];
    then Output=image000000"$t".png;
fi

B=99
if [ "$nt" -le "$B" ];
    then Output=image0000000"$t".png;
fi

B=9
if [ "$nt" -le "$B" ];
    then Output=image00000000"$t".png;
fi
  
echo $nt
t=$(($nt/100));
  
gnuplot << EOF

# initialisation du terminal
reset
set term pngcairo
set pm3d map
set size sq
set cbrange[-1.5:1.5]
set xr[0:128]
set yr[0:128]

set title "t = $t"
set output "$Output"

splot "./data.txt.$nt" u 1:2:5

EOF

nt=$(($nt+1000))
done

#fps = 30
ffmpeg -r 30 -f image2 -pattern_type glob -i "*?png" -vcodec libx264 -crf 20 -pix_fmt yuv420p output.mp4
rm *png

# fps : frame per sec
# vbitrate sets the quality


