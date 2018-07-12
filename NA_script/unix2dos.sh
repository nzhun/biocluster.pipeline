file=$1;
tr -s '\15\32' '\n' <  $file  >${file}.temp
mv ${file}.temp $file