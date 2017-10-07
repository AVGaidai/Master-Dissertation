#!/bin/bash

# If the number of arguments is less than 2
if [[ $# -lt 2 ]] ;
then
    exit
fi

# The first argument is the number of matrix rows
if [[ $1 ]] ;
then
    let ROWS=$1
fi

# The second argument is the number of matrix columns
if [[ $2 ]] ;
then
    let COLS=$2
fi

# The third argument is the output filename (default 'data.txt')
if [[ $3 ]] ;
then
    NAME=$3
else
    NAME="data.txt"
fi

# Creating output file
echo "" > $NAME

# Writing number of rows and columns into output file
string=$ROWS" "$COLS
echo $string > $NAME

# Generating matrix
let i=0
while [[ $i -ne $ROWS ]] ;
do
    let j=0
    string=""
    while [[ $j -ne $COLS ]] ;
    do
        # Generating sign of number
	if [[ $RANDOM%10 < 5 ]] ;
	then
	    NUM="-"
	else
	    NUM=""
	fi

        # Generating number
	let int=$RANDOM%20    # Integer part
	let frac=$RANDOM%100  # Float part
	NUM=$NUM$int"."$frac
        # Forming a row of the matrix
	string=$string$NUM" "
	let j=$j+1
    done

    # Writing row of matrix into output file
    echo $string >> $NAME

    let i=$i+1
done
