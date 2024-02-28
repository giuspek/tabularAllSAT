#!/bin/bash
for ((i=10;i<=50;i++));
do
	t=$(expr 1.5 * $i | bc)
	echo $t
done