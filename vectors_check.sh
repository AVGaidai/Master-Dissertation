#!/bin/bash

if [[ $# -ge 2 ]] ;
then
    NAME1=$1
    NAME2=$2
else
    exit 1
fi

EPS=0.0001

if [[ $3 ]] ;
then
    EPS=$3
fi
    
VEC1=`cat $NAME1`
VEC2=`cat $NAME2`

let i=0
for word in $VEC1 ;
do
    ARR1[$i]=$word
    let i=$i+1
done

let i=0
for word in $VEC2 ;
do
    ARR2[$i]=$word
    let i=$i+1
done


if [[ ${ARR1[0]} -ne ${ARR2[0]} ]] ;
then
    echo "Vectors is not equal! (Accuracy: $EPS)"
    exit 2
fi

let i=1

while [[ $i -lt ${ARR1[0]} ]] ;
do
    dif="${ARR1[$i]} - ${ARR2[$i]}"
    dif=`echo $dif | bc -l`
    cmp=`echo "$dif < 0" | bc -l`
    
    if [[ $cmp -eq 1 ]] ;
    then
        dif=`echo "0 - $dif" | bc -l`
    fi
    
    cmp=`echo "$dif > $EPS" | bc -l`
    if [[ $cmp -eq 1 ]] ;
    then
        echo "Vectors is not equal! (Accuracy: $EPS)"
        exit 2
    fi

    let i=$i+1
done

echo "Vectors is equal! (Accuracy: $EPS)"

