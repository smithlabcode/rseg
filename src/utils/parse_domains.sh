#!/usr/bin/env bash

for f in $@
do
	classes=$(cut -f4 $f|sort|uniq)
	for c in $classes
	do
		t=${c/:/-}
		cat $f |grep $c > ${f/.bed/-$t.bed}
	done
done
