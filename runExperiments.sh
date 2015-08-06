for i in *.txt; 
do  echo "$i"; ./TrianglePeel.exe  < "$i" > "$i-Triangle-results"; done
