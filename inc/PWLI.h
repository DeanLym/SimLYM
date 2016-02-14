#ifndef __PWLI_H__
#define __PWLI_H__


double PiecewiseLinearInter(int n, double *pX1, double *pY1, double x){
	double ans = 0.0;
	double *pLx = new double[n];
	double xj,xj1,xj2; //xj1 means x(j-1) xj2 -- x(j+1)
		
	// for j = 0
	int j = 0;
	xj = pX1[j];
	xj2 = pX1[j+1];
	if (x >= xj && x <= xj2){
		pLx[j] = (x - xj2) / (xj - xj2);
		}
	else{
		pLx[j] = 0;
	}
	ans += pLx[j]*pY1[j];
			// for j=1:n-2
	for (j = 1; j < n-1; ++j){ // The j-th PiecewiseLinear base function of interpolation Lxj
		xj = pX1[j];
		xj1 = pX1[j-1];
		xj2 = pX1[j+1];
		if (x >= xj1 && x <= xj){
			pLx[j] = (x - xj1) / (xj - xj1);
		}
		else if (x > xj && x <= xj2){
			pLx[j] = (x - xj2) / (xj - xj2);
		}
		else{
			pLx[j] = 0;
		}
		ans += pLx[j]*pY1[j];
	}

	//for j=n-1;
	j = n-1 ;
	xj = pX1[j];
	xj1 = pX1[j-1];
	if (x >= xj1 && x <= xj){
		pLx[j] = (x - xj1) / (xj - xj1);
	}
	else{
		pLx[j] = 0;
	}
	ans += pLx[j]*pY1[j];

	delete[] pLx;
	return ans;
}


#endif