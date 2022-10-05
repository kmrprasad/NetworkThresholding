#!/bin/bash

matlab_exec=matlab
mkdir -p 'Data'

#loop for the repeated tries
for try in 1 #{1..15..1} 
do
	#loop for different mut
	#for mut in 2 #3 4 5 
	#do
		mut=2
		muw=2
		echo CreateNetwork
		matlab -nodisplay -r 'CreatePowerLawNetwork('363','24','107','$mut','$muw','12','32');quit'
		echo Network created
		# loop for number of subject
		for nbsubj in 100
		do	
			# loop over the SNR
			for j in  35 
			do 
				#prepare for the results file name
				nameFolder="./Data/Simulation_"
				foo="Subj_"
				foo2="SNR"
				foo3="_Try"
				foo4="_Mut"
				echo $nameFolder $nbsubj $foo $foo2 $j $foo3 $try $foo4 $mut
				mkdir $nameFolder$nbsubj$foo$foo2$j$foo3$try$foo4$mut

				#loop for creating subj simulation
				echo Before loop Subject
				for (( isubj=1; isubj<=$nbsubj; isubj++ ))
				do
					echo $foo$isubj $nbsubj $j $try
					Rscript CreateTS.R $j $nameFolder$nbsubj$foo$foo2$j$foo3$try$foo4$mut $isubj
				done 
				echo After loop Subject
				foofinal=$nameFolder$nbsubj$foo$foo2$j$foo3$try
				echo $i $j $foofinal
			#echo Before Matlab compute Network
			#matlab -nodisplay -r 'ComputePacoNetworkSimulation('$nbsubj','$j','$try','$mut');quit'
			#echo After Matlab compute Network
			done
		done
	#done
done
