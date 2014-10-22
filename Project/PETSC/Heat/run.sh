#!/bin/bash
# petsc/scr/sles/examples/tutorials/ex2.c  or heat_petsc.c
set appl='./heat_petsc'
set options='-ksp_monitor -sles_view -log_summary -options_table -options_left -m 10 -n 10'
set num='0'
for ((np= 1; np<=8; np *=2))  do
	echo "$np"
	for ksptype in "gmres" "bcgs" "tfqmr";  do 
		echo "$ksptype"
		set pctypes_parallel='bjacobi asm'
		set pctypes_serial='ilu'
		if [[ "$np" = 1]]; then 
			set pctype_list="$pctypes_serial $pctypes_parallel"
		else
			set pctype_list="$pctypes_parallel"
		fi

		for pctype in "$pctype_list"; do
			if [["$pctype" = ilu]] then
				for ((level=0; level<=2; level++))
				do
					echo '' 
					echo '**********Beginning new run***********'
					echo ''
					set cmd="mpirun -np $np $appl -ksp_type $ksptype -pc_type $pctype -pc_ilu_levels $level $options"
					set num='expr $num+1'; echo"$num: $cmd"
					eval $cmd
				done
			else
				for subpctype in "jacobi" "sor" "ilu"
				do
				if [["$subpctype"=ilu]] then
					for ((level=0; level<=2; level++))
					do
					echo ''
					echo '*********Begining new run***********'
					set cmd="mpirun -np $np $appl -ksp_type $ksptype -pc_type $pctype -sub_ksp_type preonly -sub_pc_type $subpctype -sub_pc_ilu_levels $level $options"
					set num='expr $num+1'; echo"$num:$cmd"
					eval $cmd
					done
				else
					echo ''
					echo '*********Begining new run*********'
					set cmd="mpirun -np $np $appl -ksp_type $ksptype -pc_type $pctype -sub_ksp_type preonly -sub_pc_type $subpctype $options"
					set num='expr $num+1'; echo "$num:$cmd"
					eval $cmd
				fi
				done
			fi

		done
	done
done

