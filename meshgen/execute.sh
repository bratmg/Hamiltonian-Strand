clear
echo
echo "Clearing data from output..."
rm ./output/*.dat
rm ./output/*.tecdat
rm ./output/preplot/*.plt
rm -r ./output/domain*
rm *.dat
echo 
echo "Compiling..."
echo "====================================================================="
make
echo "====================================================================="
echo "Continue? (Enter to continue, Ctrl + C to exit)"
echo "====================================================================="
read dummy_variable
./bin/meshgen
echo

