#!/bin/bash
cd src

make
echo "Run tests? (y/n)"
read yntest
if [ "$yntest" == "y" ] # If y, compile and run both c++ codes with O3 optimization
then
  ./testcode.out
fi


echo "Generate new data before plotting? (y/n)"
read yn
if [ "$yn" == "y" ] # If y, compile and run both c++ codes with O3 optimization
then
  echo "Please provide matrix dimension"
  read dim
  echo "Please provide tolerance"
  read tol
  echo "Please provide RhoN"
  read rhoN
  ./mainprog.out $dim $tol $rhoN
fi
