#!/bin/bash

# Fq file that needs to be split
FQ_FILE=$1
# The two files that it will be split into
FIRST_FILE=$2
SECOND_FILE=$3
COUNT=$#
if [[ $COUNT != 3 ]]; then
  echo "Incorrect number of arguments"
  echo "Correct format: DataPrepFileSplit.sh [FQ_FILE_TO_BE_SPLIT] [FIRST_SPLIT_FILE_NAME] [SECOND_SPLIT_FILE_NAME]"
  echo "Please provide these as absolute file paths."
  exit 2
fi
  
# Seqanswers.com
cat $FQ_FILE | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > $FIRST_FILE
if [[ $? -ne 0 ]]; then
  # If awk didn't return 0
  echo "Dataprep file splitting script has failed."
  exit 1
fi

cat $FQ_FILE | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > $SECOND_FILE
if [[ $? -ne 0 ]]; then
  # If awk didn't return 0
  echo "Dataprep file splitting script has failed."
  exit 1
fi

# Good exit
exit 0