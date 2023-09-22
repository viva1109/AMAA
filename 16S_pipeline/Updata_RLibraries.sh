#!/bin/bash
rm -r RFunctions_tmp
wget -r -l 99 -np -P RFunctions_tmp "http://viva1109.duckdns.org/RFunctions_docker/FunctionsTMAT"
mkdir /usr/local/RFunctions
cp ./RFunctions_tmp/viva1109.duckdns.org/RFunctions_docker/* /usr/local/RFunctions -f -r
rm -r RFunctions_tmp