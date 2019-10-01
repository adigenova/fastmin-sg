egrep -v "^@" $1 |awk '{fwd=$1 ;ctg=$3; start=$4; getline; if(ctg != $3){}else{ v=start-$4; if(v < 0){v=v*-1;};print v}}' | head -n 900000 > $1.insert.txt
