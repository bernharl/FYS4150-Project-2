#!/bin/bash

cd src
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
  make
  ./mainprog.out $dim $tol $rhoN
fi
