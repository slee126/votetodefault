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
void cap_greater0(double cap1[], double cap2[], double cap3[],  double cap4[], double k1[], double k2[], double k3[], double k4[], 
	double m1[], double m2[], double m3[], double m4[], double finer_grid[], 
	int startadd, int N,  int &count1, int &count2, int &count3, int &count4, 
        int add1[], int add2[], int add3[], int add4[]){
	count1 = 0;
	count2 = 0;
	count3 = 0;
	count4 = 0;
	for (int i = startadd; i < N + startadd; i++){
		if (cap1[i] > 0){
			k1[count1] = finer_grid[i - startadd];
			m1[count1] = cap1[i];
            add1[count1] = i - startadd;
			count1++;
		}
		if (cap2[i] > 0){
			k2[count2] = finer_grid[i - startadd];
			m2[count2] = cap2[i];
            add2[count2] = i - startadd;
			count2++;
		}
		if (cap3[i] > 0){
			k3[count3] = finer_grid[i - startadd];
			m3[count3] = cap3[i];
            add3[count3] = i - startadd;
			count3++;
		}
		if (cap4[i] > 0){
			k4[count4] = finer_grid[i - startadd];
			m4[count4] = cap4[i];
            add4[count4] = i - startadd;
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
void get_kdec_slice(double kslice[], int iK, int ib, double kopt[], int prod_type, int iy, int NK, int Nb, int Ny, int Nk, int G){      
    
#pragma omp parallel for 
    for (int ii = 0; ii < G; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int address = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*4*G*Nk) + (ii*4*Nk) + (prod_type*Nk) + jj;
            kslice[ii*Nk+jj] = kopt[address];
        }//jj
    }//ii
}

//get default individual from vnd and vd iterpolation---------------------------------------------------------------------------------------------
void get_value_slice(double vnd_slice[], double vd_slice[], int iK, int ib,  double vnd[], double vd[], \
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
    double *k_grid, *mass, *surv_prob, *finer_grid, *kopt, *b_grid, *vnd, *vd, *grid_slice;
    
    k_grid = mxGetPr(prhs[0]);
    mass = mxGetPr(prhs[1]);
    surv_prob = mxGetPr(prhs[2]);
    finer_grid = mxGetPr(prhs[3]);
    kopt = mxGetPr(prhs[4]);
    vnd = mxGetPr(prhs[5]);
    vd = mxGetPr(prhs[6]);
    grid_slice = mxGetPr(prhs[7]);
    
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

    double *bigK_realized, *countMat, *default_mat, *vote_nd, *vote_d; 
    bigK_realized = (double *)mxMalloc(NK*Nb*Ny*sizeof(double));
    countMat = (double *)mxMalloc(NK*Nb*Ny*sizeof(double));
    default_mat = (double *)mxMalloc(NK*Nb*Ny*sizeof(double));
    vote_nd = (double *)mxMalloc(NK*Nb*Ny*sizeof(double));
    vote_d = (double *)mxMalloc(NK*Nb*Ny*sizeof(double));
    
    double *default_slice0, *capital_slice0, *nd_slice0, *default_slice1, *capital_slice1, *nd_slice1, 
            *default_slice2, *capital_slice2, *nd_slice2, *default_slice3, *capital_slice3, *nd_slice3;   
    
    default_slice0 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice0 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    nd_slice0 = (double *)mxCalloc(G*N_finek, sizeof(double));     
    default_slice1 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice1 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    nd_slice1 = (double *)mxCalloc(G*N_finek, sizeof(double));     
    default_slice2 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice2 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    nd_slice2 = (double *)mxCalloc(G*N_finek, sizeof(double));     
    default_slice3 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    capital_slice3 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    nd_slice3 = (double *)mxCalloc(G*N_finek, sizeof(double)); 
    
    for (int ib = 0; ib < Nb; ib++){
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

                double *vnd0, *vnd1, *vnd2, *vnd3;
                vnd0 = (double *)mxMalloc(Nk*G*sizeof(double));
                vnd1 = (double *)mxMalloc(Nk*G*sizeof(double));
                vnd2 = (double *)mxMalloc(Nk*G*sizeof(double));
                vnd3 = (double *)mxMalloc(Nk*G*sizeof(double));

                double *vd0, *vd1, *vd2, *vd3;
                vd0 = (double *)mxMalloc(Nk*G*sizeof(double));
                vd1 = (double *)mxMalloc(Nk*G*sizeof(double));
                vd2 = (double *)mxMalloc(Nk*G*sizeof(double));
                vd3 = (double *)mxMalloc(Nk*G*sizeof(double));

                double *sum0, *sum1, *sum2, *sum3, *popSum;
                sum0 = (double *) mxMalloc(N_finek*sizeof(double));
                sum1 = (double *) mxMalloc(N_finek*sizeof(double));
                sum2 = (double *) mxMalloc(N_finek*sizeof(double));
                sum3 = (double *) mxMalloc(N_finek*sizeof(double));
                popSum = (double *) mxMalloc(N_finek*sizeof(double));

                //lookup kdec ---------------------------------------------------------------------------------------------
                get_kdec_slice(kdec0, iK, ib, kopt, 0, iy, NK, Nb, Ny, Nk, G);      
                get_kdec_slice(kdec1, iK, ib, kopt, 1, iy, NK, Nb, Ny, Nk, G);      
                get_kdec_slice(kdec2, iK, ib, kopt, 2, iy, NK, Nb, Ny, Nk, G);      
                get_kdec_slice(kdec3, iK, ib, kopt, 3, iy, NK, Nb, Ny, Nk, G);      

                get_value_slice(vnd0, vd0, iK, ib, vnd, vd, 0, iy, NK, Nb, Ny, Nk, G);
                get_value_slice(vnd1, vd1, iK, ib, vnd, vd, 1, iy, NK, Nb, Ny, Nk, G);
                get_value_slice(vnd2, vd2, iK, ib, vnd, vd, 2, iy, NK, Nb, Ny, Nk, G);
                get_value_slice(vnd3, vd3, iK, ib, vnd, vd, 3, iy, NK, Nb, Ny, Nk, G);

                double *k00, *k01, *k02, *k03, *mass0, *mass1, *mass2, *mass3;
                int *add0, *add1, *add2, *add3;
                k00 = (double *)mxMalloc(N_finek*sizeof(double));
                k01 = (double *)mxMalloc(N_finek*sizeof(double));
                k02 = (double *)mxMalloc(N_finek*sizeof(double));
                k03 = (double *)mxMalloc(N_finek*sizeof(double));

                mass0 = (double *)mxMalloc(N_finek*sizeof(double));
                mass1 = (double *)mxMalloc(N_finek*sizeof(double));
                mass2 = (double *)mxMalloc(N_finek*sizeof(double));
                mass3 = (double *)mxMalloc(N_finek*sizeof(double));
              
                add0 = (int *)mxMalloc(N_finek*sizeof(int));
                add1 = (int *)mxMalloc(N_finek*sizeof(int));
                add2 = (int *)mxMalloc(N_finek*sizeof(int));
                add3 = (int *)mxMalloc(N_finek*sizeof(int));
                
                double tempcount = 0;
                double default_pop = 0.0;
                double nodefault_pop = 0.0;
                
                for (int ia = 0; ia < G; ia++){                
                    int count0 = 0, count1 = 0, count2 = 0, count3 = 0;
                    int startadd = N_finek*ia;
                    int startadd2 = N_finek*(ia + 1);
                    
                    cap_greater0(cap_dist0, cap_dist1, cap_dist2, cap_dist3, k00, k01, k02, k03, \
                        mass0, mass1, mass2, mass3, finer_grid, startadd, N_finek, count0, count1, count2, count3, add0, add1, add2, add3);

                    if(count3 > (int)tempcount){
                        tempcount = (double)count3;
                        countMat[(iy*NK*Nb)+(iK*Nb)+ib] = (double)count3;
                    }
                    
                    //low low
                    for (int ic = 0; ic < count1; ic++){
                        int ind_low1, ind_high1;
                        double weight_low, weight_high;
                        find_between(k_grid, Nk, k00[ic], ind_low1, ind_high1, weight_low, weight_high);
                        double newmass = mass0[ic]*surv_prob[ia];
                        
                        if (ia < G - 1){
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
                        }

                        double value_def_weight = vd0[(ia*Nk)+ind_low1]*weight_low + vd0[(ia*Nk)+ind_high1]*weight_high;
                        double value_nodef_weight = vnd0[(ia*Nk)+ind_low1]*weight_low + vnd0[(ia*Nk)+ind_high1]*weight_high;                            
                        if (value_def_weight > value_nodef_weight){
                            default_pop += mass0[ic];              
                            if ((iy == y_slice) & (iK == k_slice) & (ib == b_slice))
                                default_slice0[startadd+add0[ic]] = mass0[ic];
                        }
                        else{
                            nodefault_pop += mass0[ic];
                            if ((iy == y_slice) & (iK == k_slice) & (ib == b_slice))
                                nd_slice0[startadd+add0[ic]] = mass0[ic];
                        }
                    }//lowlow

                    //low high
                    for (int ic = 0; ic < count1; ic++){
                        int ind_low1, ind_high1;
                        double weight_low, weight_high;
                        find_between(k_grid, Nk, k01[ic], ind_low1, ind_high1, weight_low, weight_high);
                        double newmass = mass1[ic]*surv_prob[ia];

                        if (ia < G - 1){
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
                        }    
                        
                        double value_def_weight = vd1[(ia*Nk)+ind_low1]*weight_low + vd1[(ia*Nk)+ind_high1]*weight_high;
                        double value_nodef_weight = vnd1[(ia*Nk)+ind_low1]*weight_low + vnd1[(ia*Nk)+ind_high1]*weight_high;                           
                        if (value_def_weight > value_nodef_weight){
                            default_pop += mass1[ic];              
                            if (iy == y_slice & iK == k_slice & ib == b_slice)
                                default_slice1[startadd+add1[ic]] = mass1[ic];
                        }
                        else{
                            nodefault_pop += mass1[ic];
                            if (iy == y_slice & iK == k_slice & ib == b_slice)
                                nd_slice1[startadd+add1[ic]] = mass1[ic];
                        }
                    }//lowhigh

                    //highlow
                    for (int ic = 0; ic < count2; ic++){
                        int ind_low1, ind_high1;
                        double weight_low, weight_high;
                        find_between(k_grid, Nk, k02[ic], ind_low1, ind_high1, weight_low, weight_high);
                        double newmass = mass2[ic]*surv_prob[ia];
                        
                        if (ia < G - 1){
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
                        }                        
                                               
                        double value_def_weight = vd2[(ia*Nk)+ind_low1]*weight_low + vd2[(ia*Nk)+ind_high1]*weight_high;
                        double value_nodef_weight = vnd2[(ia*Nk)+ind_low1]*weight_low + vnd2[(ia*Nk)+ind_high1]*weight_high;                       
                        if (value_def_weight > value_nodef_weight){
                            default_pop += mass2[ic];              
                            if (iy == y_slice & iK == k_slice & ib == b_slice)
                                default_slice2[startadd+add2[ic]] = mass2[ic];
                        }
                        else{
                            nodefault_pop += mass2[ic];
                            if (iy == y_slice & iK == k_slice & ib == b_slice)
                                nd_slice2[startadd+add2[ic]] = mass2[ic];
                        }
                    }//highlow

                    //highhigh
                    for (int ic = 0; ic < count3; ic++){
                        int ind_low1, ind_high1;
                        double weight_low, weight_high;
                        find_between(k_grid, Nk, k03[ic], ind_low1, ind_high1, weight_low, weight_high);
                        double newmass = mass3[ic]*surv_prob[ia];

                        if (ia < G - 1){
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
                        }           
                        
                        double value_def_weight = vd3[(ia*Nk)+ind_low1]*weight_low + vd3[(ia*Nk)+ind_high1]*weight_high;
                        double value_nodef_weight = vnd3[(ia*Nk)+ind_low1]*weight_low + vnd3[(ia*Nk)+ind_high1]*weight_high;                           
                        if (value_def_weight > value_nodef_weight){
                            default_pop += mass3[ic];              
                            if ((iy == y_slice) & (iK == k_slice) & (ib == b_slice))
                                default_slice3[startadd+add3[ic]] = mass3[ic];
                        }
                        else{
                            nodefault_pop += mass3[ic];
                            if ((iy == y_slice) & (iK == k_slice) & (ib == b_slice))
                                nd_slice3[startadd+add3[ic]] = mass3[ic];
                        }
                    }//high
                }//closes ia
                
                sumAge(cap_dist0, cap_dist1, cap_dist2, cap_dist3, sum0, sum1, sum2, sum3, popSum, G, N_finek);
                                
                bigK_realized[(iy*NK*Nb)+(iK*Nb)+ib] = ddot(popSum, finer_grid, N_finek);
                default_mat[(iy*NK*Nb)+(iK*Nb)+ib] = default_pop > nodefault_pop ? 1.0 : 0.0;
                vote_nd[(iy*NK*Nb)+(iK*Nb)+ib] = nodefault_pop;
                vote_d[(iy*NK*Nb)+(iK*Nb)+ib] = default_pop;
                
                if ((iy == y_slice) & (iK == k_slice) & (ib == b_slice)){
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
                mxFree(add0);
                mxFree(add1);
                mxFree(add2);
                mxFree(add3);                
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
                mxFree(vd0);
                mxFree(vd1);
                mxFree(vd2);
                mxFree(vd3);              
                mxFree(vnd0);
                mxFree(vnd1);
                mxFree(vnd2);
                mxFree(vnd3);              
            } //closes iy
        } //closes iK
    }//closes ib

    std::cout << "distribution no default done " << std::endl;

    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], bigK_realized);
    mxSetM(plhs[0], NK*Nb*Ny);
    mxSetN(plhs[0], 1);
    
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[1], countMat);
    mxSetM(plhs[1], NK*Nb*Ny);
    mxSetN(plhs[1], 1);
    
    plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[2], default_mat);
    mxSetM(plhs[2], NK*Nb*Ny);
    mxSetN(plhs[2], 1);
    
    plhs[3] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[3], default_slice0);
    mxSetM(plhs[3], G*N_finek);
    mxSetN(plhs[3], 1);
      
    plhs[4] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[4], default_slice1);
    mxSetM(plhs[4], G*N_finek);
    mxSetN(plhs[4], 1);
    
    plhs[5] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[5], default_slice2);
    mxSetM(plhs[5], G*N_finek);
    mxSetN(plhs[5], 1);
    
    plhs[6] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[6], default_slice3);
    mxSetM(plhs[6], G*N_finek);
    mxSetN(plhs[6], 1);
    
    plhs[7] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[7], nd_slice0);
    mxSetM(plhs[7], G*N_finek);
    mxSetN(plhs[7], 1);
      
    plhs[8] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[8], nd_slice1);
    mxSetM(plhs[8], G*N_finek);
    mxSetN(plhs[8], 1);
    
    plhs[9] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[9], nd_slice2);
    mxSetM(plhs[9], G*N_finek);
    mxSetN(plhs[9], 1);    
    
    plhs[10] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[10], nd_slice3);
    mxSetM(plhs[10], G*N_finek);
    mxSetN(plhs[10], 1);
    
    plhs[11] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[11], capital_slice0);
    mxSetM(plhs[11], G*N_finek);
    mxSetN(plhs[11], 1);
    
    plhs[12] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[12], capital_slice1);
    mxSetM(plhs[12], G*N_finek);
    mxSetN(plhs[12], 1);
    
    plhs[13] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[13], capital_slice2);
    mxSetM(plhs[13], G*N_finek);
    mxSetN(plhs[13], 1);
    
    plhs[14] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[14], capital_slice3);
    mxSetM(plhs[14], G*N_finek);
    mxSetN(plhs[14], 1);
    
    plhs[15] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[15], vote_nd);
    mxSetM(plhs[15], NK*Nb*Ny);
    mxSetN(plhs[15], 1);
    
    plhs[16] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[16], vote_d);
    mxSetM(plhs[16], NK*Nb*Ny);
    mxSetN(plhs[16], 1);
}