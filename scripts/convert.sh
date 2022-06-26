#!/bin/bash
cd ../figures/
for file in *.pdf;
	do inkscape "$file" --export-dpi=300 "${file%pdf}eps";
done

