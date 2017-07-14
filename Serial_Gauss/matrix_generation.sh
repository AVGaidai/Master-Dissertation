#!/bin/bash

if [[ $1 ]] ;
then
    NAME=$1
else
    NAME="data.txt"
fi

echo "" > $NAME

if [[ $2 ]] ;
then
    let ROWS=$2
else
    let ROWS=3
fi

if [[ $3 ]] ;
then
    let COLS=$3
else
    let COLS=4
fi


string=$ROWS" "$COLS
echo $string > $NAME

let i=0

while [[ $i -ne $ROWS ]] ;
do
    let j=0
    string=""
    while [[ $j -ne $COLS ]] ;
    do
	if [[ $RANDOM%10 < 5 ]] ;
	then
	    NUM="-"
	else
	    NUM=""
	fi
	let int=$RANDOM%20
	let frac=$RANDOM%100
	NUM=$NUM$int"."$frac
	string=$string$NUM" "
	let j=$j+1
    done

    echo $string >> $NAME

    let i=$i+1
done
