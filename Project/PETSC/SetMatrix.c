#include"petscksp.h"
#include"petscda.h"

Mat FormLaplacianDA2d(DA grid, int n)
{
	Mat A;
	int r, i, j, is, js, ie, je, in, jn, nelm;
	MatStencil cols[5], row;
	double h, oneByh2, vals[5];

	h = 1.0/(n+1); oneByh2 = 1.0/(h*h);

	DAGetMatrix(grid, MATMPIAIJ, &A);
	DAGetCorners(grid, &is, &js, 0, &in, &jn, 0);
	ie = is + in - 1;
	je = js + jn - 1;
	for(i=is; i<=ie; i++){
		for(j=js; j<=je; j++){
			row.j = j;
			row.i = i;
			nelm = 0;
			if(j-1 > 0){
				vals[nelm] = oneByh2;
				cols[nelm].j = j-1;
				cols[nelm++].i = i;
			}
			if(i-1>0){
				vals[nelm] = oneByh2;
				cols[nelm].j = j;
				cols[nelm++].i = i-1;
			}
			vals[nelm] = -4*oneByh2;
			cols[nelm].j =j;
			cols[nelm++].i = i;
			if(i+1 < n-1){
				vals[nelm] = oneByh2;
				cols[nelm].j = j;
				cols[nelm++].i = i+1;
			}
			if(j+1 < n-1){
				vals[nelm] = oneByh2;
				cols[nelm].j = j+1;
				cols[nelm++].i = i;
			}
			MatSetValuesStencil(A,1,&row, nelm, cols, vals, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	return A;
}


