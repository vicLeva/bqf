TARGET_STRING="tests failed"
SCRIPT="ctest"

for i in {1..100}
do 
    echo  "Run #$i"
    OUTPUT=$($SCRIPT)
    if echo "$OUTPUT" | grep -q "$TARGET_STRING"; then
        echo "Failure on run #$i"
    fi 
done