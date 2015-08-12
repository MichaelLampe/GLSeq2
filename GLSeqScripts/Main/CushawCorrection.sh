#!/bin/bash

# File that needs to have certain sequences removed
UNCLEAN_FILE=$1
# File that will be produced as a result
CLEAN_FILE=$2
COUNT=$#
if [[ $COUNT != 2 ]]; then
  echo "Incorrect number of arguments"
  echo "Correct format: CushawCorrection.sh [UNCLEAN FILE] [CLEAN FILE TO BE MADE]"
  echo "Please provide these as absolute file paths."
  exit 1
fi
  
# https://www.biostars.org/p/108702/
awk '!/\t\*\t/' $UNCLEAN_FILE > $CLEAN_FILE

if [[ $? -ne 0 ]]; then
  # If awk didn't return 0
  echo "Cushaw_PE countable.sam file correction shell script has failed."
  exit 1
fi
# Good exit
exit 0