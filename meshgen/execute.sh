clear
echo "====================================================================="
echo "Clearing data from output..."
echo "====================================================================="
echo
rm ./output/*.dat
rm ./output/*.tecdat
rm ./output/preplot/*.plt
rm *.dat
echo 
echo "====================================================================="
echo "Compiling..."
echo "====================================================================="
echo
make
echo
echo "====================================================================="
echo "Continue? (Enter to continue, Ctrl + C to exit)"
echo "====================================================================="
read dummy_variable
./bin/meshgen
echo
