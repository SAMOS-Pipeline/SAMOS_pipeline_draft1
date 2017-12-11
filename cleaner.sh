#!/bin/bash

echo "This script will delete the direcorites LMask1 and LMask2 along with their associated databases."
read -p  "Would you like to continue? (yes) " answer
answer=${answer:-yes}
DIR=$pwd
echo $DIR
if [ $answer =  "yes" ] || [ $answer = "y" ] 
 then
  echo "Files and directories will be cleaned... "
  if [ -d "LMask1" ]; then
   echo "Deleting LMask1"
   rm -r "LMask1"
  fi
  if [ -d "LMask2" ]; then
   echo "Deleting LMask2"
   rm -r "LMask2"
  fi
 echo "Deleting $DIR"/LMask*.db
 find . -type f -name "*.db" -exec rm -f {} \; 
 echo "Deleting $DIR"/LMask*.SMF
 find . -type f -name "*.SMF" -exec rm -f {} \;
 echo "Deleted associated LMask files"
  
else
 echo "Exiting without cleaning. "
 exit
fi
