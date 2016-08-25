#include "mex.h"
#include <math.h>
#include "matrix.h"
#include<stdio.h>
#include <algorithm>
#include <iostream>
#include <omp.h>

//finds the weight and index for linear interpolation---------------------------------------------------------------------------------------------
void find_between(double inarr1[], int N, double val, int &ind_low, int &ind_high, double &weight_low, double &weight_high){
	ind_low = 0;
	ind_high = 1;
	weight_low = 1.0;
	weight_high = 0;

	for (int i = 0; i < N; i++){
		if (inarr1[i] > val){
			ind_high = i;
			ind_low = i - 1;
			weight_high = (val - inarr1[ind_low])/(inarr1[ind_high] - inarr1[ind_low]);
			weight_low = 1 - weight_high;
			break;
		}
		if (i == N-1){
			ind_high = i;
			ind_low = i - 1;
			weight_low = 0.0;
			weight_high = 1.0;
		}
	}
}

//finds numerical index---------------------------------------------------------------------------------------------
int find_in(double inarr1[], int N, double val){
	for (int i = 0; i < N; i++){
		if (inarr1[i] == val){
            return i;
		}
    }
}

//find capital points where mass is greater than 0
void cap_greater0_def(double cap1[], double cap2[], double cap3[],  double cap4[], double k1[], double k2[], double k3[], double k4[], 
	double m1[], double m2[], double m3[], double m4[], double finer_grid[], 
	int startadd, int N,  int &count1, int &count2, int &count3, int &count4){
	count1 = 0;
	count2 = 0;
	count3 = 0;
	count4 = 0;
	for (int i = startadd; i < N + startadd; i++){
		if (cap1[i] > 0){
			k1[count1] = finer_grid[i - startadd];
			m1[count1] = cap1[i];
			count1++;
		}
		if (cap2[i] > 0){
			k2[count2] = finer_grid[i - startadd];
			m2[count2] = cap2[i];
			count2++;
		}
		if (cap3[i] > 0){
			k3[count3] = finer_grid[i - startadd];
			m3[count3] = cap3[i];
			count3++;
		}
		if (cap4[i] > 0){
			k4[count4] = finer_grid[i - startadd];
			m4[count4] = cap4[i];
			count4++;
		}
	}
}

//sum wealth across age 
void sumAge(double inarr1[], double inarr2[], double inarr3[], double inarr4[], double outarr1[], double outarr2[], 
	double outarr3[], double outarr4[], double outarr5[], int G, int capN){
	for (int i = 0; i < capN; i++){
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		double temp4 = 0.0;
		for (int j = 0; j < G; j++){
			temp1 += inarr1[j*capN + i];
			temp2 += inarr2[j*capN + i];
			temp3 += inarr3[j*capN + i];
			temp4 += inarr4[j*capN + i];
		}
		outarr1[i] = temp1;
		outarr2[i] = temp2;
		outarr3[i] = temp3;
		outarr4[i] = temp4;
		outarr5[i] = temp1 + temp2 + temp3 + temp4;
	}
}

// dot product of arr1 and arr2---------------------------------------------------------------------------------------------
double ddot(double inarr1[], double inarr2[], int N){
	double temp = 0;
	for (int i = 0; i < N; i++)
		temp += inarr1[i]*inarr2[i];
	return temp;
}

//interpolation bilinear with valueslice between k (xdim) and b (ydim)---------------------------------------------------------------------------------------------
void bilinear_kopt(double kslice[], double k, double b, double aggregateK[], double b_grid[], double kopt[], \
        int prod_type, int aggShock_ind, int NK, int Nb, int Ny, int Nk, int G){      
    int ind_low_k, ind_high_k, ind_low_b, ind_high_b;
    double weight_low_k, weight_high_k, weight_low_b, weight_high_b;
    
    find_between(aggregateK, NK, k, ind_low_k, ind_high_k, weight_low_k, weight_high_k);
    find_between(b_grid, Nb, b, ind_low_b, ind_high_b, weight_low_b, weight_high_b);
    
    double denom1 = (aggregateK[ind_high_k] - aggregateK[ind_low_k])*(b_grid[ind_high_b] - b_grid[ind_low_b]);
    double xdev2 = aggregateK[ind_high_k] - k;
    double xdev1 = k - aggregateK[ind_low_k];
    double ydev2 = b_grid[ind_high_b] - b;
    double ydev1 = b - b_grid[ind_low_b];
    double tempy[2] = {ydev2, ydev1};
    
#pragma omp parallel for 
    for (int ii = 0; ii < G; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int ind11 = (ind_low_b*NK*Ny*4*G*Nk) + (ind_low_k*Ny*4*G*Nk) + (aggShock_ind*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            long long int ind21 = (ind_low_b*NK*Ny*4*G*Nk) + (ind_high_k*Ny*4*G*Nk) + (aggShock_ind*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;           
            long long int ind12 = (ind_high_b*NK*Ny*4*G*Nk) + (ind_low_k*Ny*4*G*Nk) + (aggShock_ind*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            long long int ind22 = (ind_high_b*NK*Ny*4*G*Nk) + (ind_high_k*Ny*4*G*Nk) + (aggShock_ind*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            
            double temp1[2] = {xdev2*kopt[ind11] + xdev1*kopt[ind21], xdev2*kopt[ind12] + xdev1*kopt[ind22]};            
            double num = ddot(temp1, tempy, 2);
            kslice[ii*Nk + jj] = num/denom1;
        }
    }
}

//lookup kdec slice---------------------------------------------------------------------------------------------
void get_kdec_slice_def(double kslice[], int iK, double kopt[], int prod_type, int iy, int NK, int Nb, int Ny, int Nk, int G){      
    
#pragma omp parallel for 
    for (int ii = 0; ii < G; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int address = (iK*Ny*G*4*Nk) + (iy*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            kslice[ii*Nk+jj] = kopt[address];
        }//jj
    }//ii
}

//get default individual from vnd and vd iterpolation---------------------------------------------------------------------------------------------
void get_value_slice(double vnd_slice[], double vd_slice[], int iK, double aggregateK[], int ib, double bgrid[], double vnd[], double vd[], \
        int prod_type, int iy, int NK, int Nb, int Ny, int Nk, int G){      
    
#pragma omp parallel for 
    for (int ii = 0; ii < G; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int add_nd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            long long int add_d = (iK*Ny*G*4*Nk) + (iy*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            
            vnd_slice[ii*Nk+jj] = vnd[add_nd];
            vd_slice[ii*Nk+jj] = vd[add_d];
        }//jj
    }//ii
}

//copies array ---------------------------------------------------------------------------------------------
void copyArray(double inarr[], double outarr[], int Nrow, int Ncol){
	for (int i = 0; i < Nrow*Ncol; i++)
		outarr[i] = inarr[i];
}

//MAIN-------------------------------------------------------------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){        
    double *k_grid, *mass, *surv_prob, *finer_grid, *kopt, *aggregateK, *grid_slice;
    
    k_grid = mxGetPr(prhs[0]);
    mass = mxGetPr(prhs[1]);
    surv_prob = mxGetPr(prhs[2]);
    finer_grid = mxGetPr(prhs[3]);
    kopt = mxGetPr(prhs[4]);
    grid_slice = mxGetPr(prhs[5]);
    
    double y_slice = grid_slice[0];
    double k_slice = grid_slice[1];
    double b_slice = grid_slice[2];
    
    double step = finer_grid[1] - finer_grid[0];

    const int Nk = 200;
    const int G = 30;
    const int Ny = 15;
    const int NK = 52;
    const int Nb = 25;
    
    const int N_finek = 500;
    const double probz[4] = {.9608, .0392, .0392, .9608};

    double *bigK_realized, *countMat; 
    bigK_realized = (double *)mxMalloc(NK*Ny*sizeof(double));
    countMat = (double *)mxMalloc(NK*Ny*sizeof(double));

    double *capital_slice0, *capital_slice1, *capital_slice2, *capital_slice3;
    
    capital_slice0 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice1 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice2 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice3 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    
    for (int iK = 0; iK < NK; iK++){
        for(int iy = 0; iy < Ny; iy++){    
            double *cap_dist0, *cap_dist1, *cap_dist2, *cap_dist3;
            cap_dist0 = (double *)mxCalloc(N_finek*G, sizeof(double));
            cap_dist1 = (double *)mxCalloc(N_finek*G, sizeof(double));
            cap_dist2 = (double *)mxCalloc(N_finek*G, sizeof(double));
            cap_dist3 = (double *)mxCalloc(N_finek*G, sizeof(double));

            cap_dist0[0] = mass[0]*.5;
            cap_dist1[0] = mass[0]*.5;
            cap_dist2[0] = mass[0]*.5;
            cap_dist3[0] = mass[0]*.5;

            double *kdec0, *kdec1, *kdec2, *kdec3;
            kdec0 = (double *)mxMalloc(Nk*G*sizeof(double));
            kdec1 = (double *)mxMalloc(Nk*G*sizeof(double));
            kdec2 = (double *)mxMalloc(Nk*G*sizeof(double));
            kdec3 = (double *)mxMalloc(Nk*G*sizeof(double));

            double *sum0, *sum1, *sum2, *sum3, *popSum;
            sum0 = (double *) mxMalloc(N_finek*sizeof(double));
            sum1 = (double *) mxMalloc(N_finek*sizeof(double));
            sum2 = (double *) mxMalloc(N_finek*sizeof(double));
            sum3 = (double *) mxMalloc(N_finek*sizeof(double));
            popSum = (double *) mxMalloc(N_finek*sizeof(double));

            //lookup kdec ---------------------------------------------------------------------------------------------
            get_kdec_slice_def(kdec0, iK, kopt, 0, iy, NK, Nb, Ny, Nk, G);      
            get_kdec_slice_def(kdec1, iK, kopt, 1, iy, NK, Nb, Ny, Nk, G);      
            get_kdec_slice_def(kdec2, iK, kopt, 2, iy, NK, Nb, Ny, Nk, G);      
            get_kdec_slice_def(kdec3, iK, kopt, 3, iy, NK, Nb, Ny, Nk, G);      


            double *k00, *k01, *k02, *k03, *mass0, *mass1, *mass2, *mass3;
            k00 = (double *)mxMalloc(N_finek*sizeof(double));
            k01 = (double *)mxMalloc(N_finek*sizeof(double));
            k02 = (double *)mxMalloc(N_finek*sizeof(double));
            k03 = (double *)mxMalloc(N_finek*sizeof(double));

            mass0 = (double *)mxMalloc(N_finek*sizeof(double));
            mass1 = (double *)mxMalloc(N_finek*sizeof(double));
            mass2 = (double *)mxMalloc(N_finek*sizeof(double));
            mass3 = (double *)mxMalloc(N_finek*sizeof(double));

            double tempcount = 0;

            for (int ia = 0; ia < G-1; ia++){                
                int count0 = 0, count1 = 0, count2 = 0, count3 = 0;
                int startadd = N_finek*ia;
                int startadd2 = N_finek*(ia + 1);

                cap_greater0_def(cap_dist0, cap_dist1, cap_dist2, cap_dist3, k00, k01, k02, k03, \
                    mass0, mass1, mass2, mass3, finer_grid, startadd, N_finek, count0, count1, count2, count3);

                if(count3 > (int)tempcount){
                    tempcount = (double)count3;
                    countMat[(iy*NK)+iK] = (double)count3;
                }

                //low low
                for (int ic = 0; ic < count1; ic++){
                    int ind_low1, ind_high1;
                    double weight_low, weight_high;
                    find_between(k_grid, Nk, k00[ic], ind_low1, ind_high1, weight_low, weight_high);
                    double newmass = mass0[ic]*surv_prob[ia];

                    double kp = kdec0[(ia*Nk)+ind_low1]*weight_low + kdec0[(ia*Nk)+ind_high1]*weight_high;
                    kp = std::min(kp, finer_grid[N_finek-1]);
                    kp = std::max(kp, finer_grid[0]);
                    int ind_low = (int)(kp/step);
                    int ind_high = std::min(ind_low + 1, N_finek - 1);
                    double w_high = (kp - finer_grid[ind_low])/step;
                    double w_low = 1.0 - w_high;
                    cap_dist0[startadd2+ind_high] = cap_dist0[startadd2+ind_high] + w_high*newmass*probz[0];
                    cap_dist0[startadd2+ind_low] = cap_dist0[startadd2+ind_low] + w_low*newmass*probz[0];
                    cap_dist1[startadd2+ind_high] = cap_dist1[startadd2+ind_high] + w_high*newmass*probz[1];
                    cap_dist1[startadd2+ind_low] = cap_dist1[startadd2+ind_low] + w_low*newmass*probz[1];
                }//lowlow

                //low high
                for (int ic = 0; ic < count1; ic++){
                    int ind_low1, ind_high1;
                    double weight_low, weight_high;
                    find_between(k_grid, Nk, k01[ic], ind_low1, ind_high1, weight_low, weight_high);
                    double newmass = mass1[ic]*surv_prob[ia];

                    double kp = kdec1[(ia*Nk)+ind_low1]*weight_low + kdec1[(ia*Nk)+ind_high1]*weight_high;
                    kp = std::min(kp, finer_grid[N_finek-1]);
                    kp = std::max(kp, finer_grid[0]);

                    int ind_low = (int)(kp/step);
                    int ind_high = std::min(ind_low + 1, N_finek - 1);
                    double w_high = (kp - finer_grid[ind_low])/step;
                    double w_low = 1.0 - w_high;

                    cap_dist1[startadd2+ind_high] = cap_dist1[startadd2+ind_high] + w_high*newmass*probz[3];
                    cap_dist1[startadd2+ind_low] = cap_dist1[startadd2+ind_low] + w_low*newmass*probz[3];
                    cap_dist0[startadd2+ind_high] = cap_dist0[startadd2+ind_high] + w_high*newmass*probz[2];
                    cap_dist0[startadd2+ind_low] = cap_dist0[startadd2+ind_low] + w_low*newmass*probz[2];
                }//lowhigh

                //highlow
                for (int ic = 0; ic < count2; ic++){
                    int ind_low1, ind_high1;
                    double weight_low, weight_high;
                    find_between(k_grid, Nk, k02[ic], ind_low1, ind_high1, weight_low, weight_high);
                    double newmass = mass2[ic]*surv_prob[ia];

                    double kp = kdec2[(ia*Nk)+ind_low1]*weight_low + kdec2[(ia*Nk)+ind_high1]*weight_high;
                    kp = std::min(kp, finer_grid[N_finek-1]);
                    kp = std::max(kp, finer_grid[0]);

                    int ind_low = (int)(kp/step);
                    int ind_high = std::min(ind_low + 1, N_finek - 1);
                    double w_high = (kp - finer_grid[ind_low])/step;
                    double w_low = 1.0 - w_high;

                    cap_dist2[startadd2+ind_high] = cap_dist2[startadd2+ind_high] + w_high*newmass*probz[0];
                    cap_dist2[startadd2+ind_low] = cap_dist2[startadd2+ind_low] + w_low*newmass*probz[0];
                    cap_dist3[startadd2+ind_high] = cap_dist3[startadd2+ind_high] + w_high*newmass*probz[1];
                    cap_dist3[startadd2+ind_low] = cap_dist3[startadd2+ind_low] + w_low*newmass*probz[1];
                }//highlow

                //highhigh
                for (int ic = 0; ic < count3; ic++){
                    int ind_low1, ind_high1;
                    double weight_low, weight_high;
                    find_between(k_grid, Nk, k03[ic], ind_low1, ind_high1, weight_low, weight_high);
                    double newmass = mass3[ic]*surv_prob[ia];

                    double kp = kdec3[(ia*Nk)+ind_low1]*weight_low + kdec3[(ia*Nk)+ind_high1]*weight_high;
                    kp = std::min(kp, finer_grid[N_finek-1]);
                    kp = std::max(kp, finer_grid[0]);

                    int ind_low = (int)(kp/step);
                    int ind_high = std::min(ind_low + 1, N_finek - 1);
                    double w_high = (kp - finer_grid[ind_low])/step;
                    double w_low = 1.0 - w_high;

                    cap_dist3[startadd2+ind_high] = cap_dist3[startadd2+ind_high] + w_high*newmass*probz[3];
                    cap_dist3[startadd2+ind_low] = cap_dist3[startadd2+ind_low] + w_low*newmass*probz[3];
                    cap_dist2[startadd2+ind_high] = cap_dist2[startadd2+ind_high] + w_high*newmass*probz[2];
                    cap_dist2[startadd2+ind_low] = cap_dist2[startadd2+ind_low] + w_low*newmass*probz[2];                    
                }//high
            }//closes ia

            sumAge(cap_dist0, cap_dist1, cap_dist2, cap_dist3, sum0, sum1, sum2, sum3, popSum, G, N_finek);

            bigK_realized[(iy*NK)+iK] = ddot(popSum, finer_grid, N_finek);

            if ((iy == y_slice) & (iK == k_slice)){
                copyArray(cap_dist0, capital_slice0, N_finek*G, 1);
                copyArray(cap_dist1, capital_slice1, N_finek*G, 1);
                copyArray(cap_dist2, capital_slice2, N_finek*G, 1);
                copyArray(cap_dist3, capital_slice3, N_finek*G, 1);                    
            }
            
            mxFree(k00);
            mxFree(k01);
            mxFree(k02);
            mxFree(k03);                
            mxFree(mass0);
            mxFree(mass1);
            mxFree(mass2);
            mxFree(mass3);                               
            mxFree(sum0);
            mxFree(sum1);
            mxFree(sum2);
            mxFree(sum3);
            mxFree(popSum);
            mxFree(cap_dist0);
            mxFree(cap_dist1);
            mxFree(cap_dist2);
            mxFree(cap_dist3);      
            mxFree(kdec0);
            mxFree(kdec1);
            mxFree(kdec2);
            mxFree(kdec3);                        
        } //closes iy
    } //closes iK

    std::cout << "distribution default done " << std::endl;

    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], bigK_realized);
    mxSetM(plhs[0], NK*Ny);
    mxSetN(plhs[0], 1);
    
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[1], countMat);
    mxSetM(plhs[1], NK*Ny);
    mxSetN(plhs[1], 1);
    
    plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[2], capital_slice0);
    mxSetM(plhs[2], G*N_finek);
    mxSetN(plhs[2], 1);
    
    plhs[3] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[3], capital_slice1);
    mxSetM(plhs[3], G*N_finek);
    mxSetN(plhs[3], 1);
    
    plhs[4] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[4], capital_slice2);
    mxSetM(plhs[4], G*N_finek);
    mxSetN(plhs[4], 1);
    
    plhs[5] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[5], capital_slice3);
    mxSetM(plhs[5], G*N_finek);
    mxSetN(plhs[5], 1);
}