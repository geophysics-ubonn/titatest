#!/bin/bash
for i in $( ls *.mag); do
   echo nomer: $i
   ./R5_mag $i
done
for i in $( ls *.pha); do
   echo nomer: $i
   ./R5_pha $i
done
