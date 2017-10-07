#!/bin/bash

# If the number of arguments is greater than or equal to 2
if [[ $# -ge 2 ]] ;
then
    NAME1=$1   # The first filename
    NAME2=$2   # The second filename
else
    exit 1
fi

# Accuracy
EPS=0.0001

# If the thrid argument is exist
if [[ $3 ]] ;
then
    EPS=$3  # Setting new value of accuracy
fi

# Reading vectors
VEC1=`cat $NAME1`
VEC2=`cat $NAME2`

# Processing of first vector
let i=0
for word in $VEC1 ;
do
    ARR1[$i]=$word
    let i=$i+1
done

# Processing of second vector
let i=0
for word in $VEC2 ;
do
    ARR2[$i]=$word
    let i=$i+1
done

# If size of the first vector is not equal size of the second vector
if [[ ${ARR1[0]} -ne ${ARR2[0]} ]] ;
then
    echo "Vectors are not equal! (Accuracy: $EPS)"
    exit 2
fi

# Comparing vectors
let i=1
while [[ $i -lt ${ARR1[0]} ]] ;
do
    # Difference of values the first vector and the second vector
    dif="${ARR1[$i]} - ${ARR2[$i]}"
    dif=`echo $dif | bc -l`
    cmp=`echo "$dif < 0" | bc -l`

    # ABS(dif)
    if [[ $cmp -eq 1 ]] ;
    then
        dif=`echo "0 - $dif" | bc -l`
    fi

    # If difference is greater than EPS 
    cmp=`echo "$dif > $EPS" | bc -l`
    if [[ $cmp -eq 1 ]] ;
    then
        echo "Vectors are not equal! (Accuracy: $EPS)"
        exit 2
    fi

    let i=$i+1
done

echo "Vectors are equal! (Accuracy: $EPS)"

