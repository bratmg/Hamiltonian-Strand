clear
echo "#####################################################################"
echo "#                                                                   #"
echo "#                     PREPLOT CONVERTOR                             #" 
echo "#                                                                   #"
echo "#####################################################################"
echo
echo "Clearing data from preplot folder..."
rm ./output/preplot/*.plt
echo
cd output
#preplotRun="/Applications/Tecplot360EX/bin/preplot"
preplotRun="/home/bharath/Softwares/Tecplot/bin/preplot"
for f in *.dat
do
   echo "Processing $f"
   $preplotRun $f 
done
mv *.plt ./preplot/
echo "#####################################################################"
echo "#                                                                   #"
echo "#                             FIN.                                  #"
echo "#                                                                   #"
echo "#####################################################################"

