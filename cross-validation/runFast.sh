#### cross-validation scheme to narrow prior guess for hyper-parameter gamma

#!/bin/bash

model='add'    ### add or epi

start=1
ende=100
stepsize=0.001 # for gamma
   	    
    echo $(date +"%d-%m-%Y") $(date +%H:%M:%S)

    for ((sim=start;sim<=ende;sim++)) 
    do	   
    	gamma=`echo "$sim*$stepsize" | bc`
	out=$model'_'$sim
	echo $out $gamma

	cp inputWG_0.par temp0.par
	sed 's/DummyOut/'$out'/g' temp0.par >temp1.par
	sed 's/DummyGamma/'$gamma'/g' temp1.par >temp2.par

	case $model in
	    add) mv temp2.par inputWG.par ;;
	    epi) sed 's/solveEpi=.false./solveEpi=.true./g' temp2.par >inputWG.par ;;
    	esac
  
	rm temp*
	
	### estimation
	../fbayesWG > $out.out
	
	
	### prediction
	cp predict_0.par temp0.par
	sed 's/DummyOut/'$out'/g' temp0.par >temp1.par
	
	case $model in
	    add) mv temp1.par predict.par ;;
	    epi) sed 's/solveEpi=.false./solveEpi=.true./g' temp1.par >predict.par ;;
    	esac
    	
	rm temp*
	./predict_v2 >> $out.out
	    
    done
   
