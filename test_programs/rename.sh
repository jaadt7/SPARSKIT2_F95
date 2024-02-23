#!/bin/bash
for file in *.f90; do
    mv -- "$file" "${file%.f90}.f95"
done

