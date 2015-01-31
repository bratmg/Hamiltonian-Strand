#!/bin/bash
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

# preplot commands (UNIX and OSX path)
#preplotRun="/Applications/Tecplot360EX/bin/preplot"
preplotRun="/home/bharath/Softwares/Tecplot/bin/preplot"

echo "Enter total number of domains: "
read nDomain

cd output
# open folders in sequence
for (( i = 0; i < $nDomain; i++ ))
do
	echo
	echo "================================================================="
	echo " Domain$i"
	echo "================================================================="

	if [ $i -lt 10 ]
	then
		folder="domain00$i"
	elif [ $i -lt 100 ]
	then
		folder="domain0$i"
	else
		folder="domain$i"
	fi
	
	# enter the domain### folder
	cd $folder

	# process
	for f in *.tecdat
	do
   	echo "Processing $f"
   	$preplotRun $f 
	done
	mkdir preplot
	mv *.plt ./preplot/

	# move out
	cd ..
done

echo "#####################################################################"
echo "#                                                                   #"
echo "#                             FIN.                                  #"
echo "#                                                                   #"
echo "#####################################################################"

