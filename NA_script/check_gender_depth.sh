#"Y:2649521-59034050"
java -jar $GATKJAR \
   -T DepthOfCoverage \
   -R $REF \
   -o $input.depth \
   -I $inputlist
  -ct 4 -ct 6 -ct 10
  -L Y:2649521-59034050 

