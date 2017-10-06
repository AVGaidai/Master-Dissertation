#!/bin/bash

if [[ $# -lt 2 ]] ;
then
    exit
fi

if [[ $1 ]] ;
then
    let ROWS=$1
fi

if [[ $2 ]] ;
then
    let COLS=$2
fi

if [[ $3 ]] ;
then
    NAME=$3
else
    NAME="data.txt"
fi

echo "" > $NAME

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
