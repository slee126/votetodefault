#include "mex.h"
#include <math.h>
#include "matrix.h"
#include<stdio.h>
#include <algorithm>
#include <iostream>
#include <omp.h>


// dot product of arr1 and arr2---------------------------------------------------------------------------------------------
double ddot(double inarr1[], double inarr2[], int N){
	double temp = 0;
	for (int i = 0; i < N; i++)
		temp += inarr1[i]*inarr2[i];
	return temp;
}

// Gets row(index) max col = (Ncol-1) of NrowxNcol Matrix---------------------------------------------------------------------------------------------
void getRow(double inarr[], double outarr[], int Nrow, int Ncol, int row){
	for (int i = 0; i < Ncol; i++){
		outarr[i] = inarr[row*Ncol + i];
	}
}

// Gets addressing for getElements---------------------------------------------------------------------------------------------
void evenSpace(int outarr[], int start, int stride, int N){
	for (int i = 0; i < N; i++)
		outarr[i] = start + stride*i;
}

//pulls out elements according rowaddress---------------------------------------------------------------------------------------------
void getElements(double inarr1[], double outarr1[], int rowAddress[], int Nrow){
	// Gets col(index) max col = (Ncol-1) of NrowxNcol Matrix
	int address;
	for (int i = 0; i < Nrow; i++){
		address = rowAddress[i];
		outarr1[i] = inarr1[address];
	}
}

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

//sums a vector---------------------------------------------------------------------------------------------
double getSum(double inarr[], int Nrow, int Ncol){
	double temp = 0;
	for (int i = 0; i < Nrow*Ncol; i++)
		temp += inarr[i];
	return temp;
}

// references the maximum maxval row num, col num of 2 dimensional array---------------------------------------------------------------------------------------------
double getNorm(double inarr1[], double inarr2[], int Nrow, int Ncol){
	double maxVal = abs(inarr1[0] - inarr2[0]);
	for (int i = 0; i < Nrow*Ncol; i++){
		if (abs(inarr1[i] - inarr2[i]) > maxVal){
			maxVal = abs(inarr1[i] - inarr2[i]);
		}
	}
	return maxVal;
}

//arr1 is the weighted matrix of arr1 and arr2 w/ weights damp and (1-damp)---------------------------------------------------------------------------------------------
void getWeighted(double inarr1[], double inarr2[], int Nrow, int Ncol, double damp){
	for (int i = 0; i < Nrow*Ncol; i++){
		inarr1[i] = inarr1[i] *(1-damp) + inarr2[i] * damp;
	}
}

//copies array no pointers---------------------------------------------------------------------------------------------
void copyArray(double inarr[], double outarr[], int Nrow, int Ncol){
	for (int i = 0; i < Nrow*Ncol; i++)
		outarr[i] = inarr[i];
}

//initializes array to value-----------------------------------------------------------------------------------------------------
void initArray(double outarr[], int Nrow, int Ncol, double val){
	for (int i = 0; i < Nrow*Ncol; i++)
		outarr[i] = val;
}

//utility
double util(double cons){
    return log(cons);
}


//MAIN----------------------------------------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){    
    //Needed Inputs 
    double *k_grid, *prob, *Z, *par, *b_grid, *vnd_init, *vd_init;
    k_grid = mxGetPr(prhs[0]);
    prob = mxGetPr(prhs[1]);
    Z = mxGetPr(prhs[2]);
    par = mxGetPr(prhs[3]);
    b_grid = mxGetPr(prhs[4]);
    vnd_init = mxGetPr(prhs[5]);
    vd_init = mxGetPr(prhs[6]);
    
    const double alpha = par[0];
    const double delta = par[1];
    const double tauk = par[2];
    const double tauw = par[3];
    const double taup = par[4];
    const double beta = par[5];
    const double loss = par[6];
    const double rf = par[7];
    const double theta = par[8];
    
    const int NK = 52;
    const int Ny = 15;
    const int Nb = 25;

    double *vo, *vo_temp, *vnd_temp, *kopt, *copt, *bopt, *default_decision, *price_mat;
    
    vo = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    copyArray(vnd_init, vo, NK*Ny, Nb);
    
    default_decision = (double *)mxCalloc(NK*Ny*Nb, sizeof(double));
    price_mat = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    initArray(price_mat, NK*Ny, Nb, 1/(1 + rf));
    
    vo_temp = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    vnd_temp = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    kopt = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    copt = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    bopt = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));

    
    double *cdef, *val_diff, *vd, *vd_temp, *kopt_d; 
    vd = (double *)mxMalloc(NK*Ny * sizeof(double));
    copyArray(vd_init, vd, NK, Ny);

    cdef = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));
    vd_temp = (double *)mxMalloc(NK*Ny * sizeof(double));
    kopt_d = (double *)mxMalloc(NK*Ny * sizeof(double));
    val_diff = (double *)mxMalloc(NK*Ny*Nb * sizeof(double));

    double damp = .5;
    double crit = 1;
    double tol = 1e-3;
    int iter = 0;

    while (crit  > tol){
		for (int iy = 0; iy < Ny; iy++){
			double z = Z[iy];
            double prob_vec[Ny];
            getRow(prob, prob_vec, Ny, Ny, iy);
            
			for (int ik = 0; ik < NK; ik++){
				 double k0 = k_grid[ik];
				 double r = (1 - tauk)*(z*alpha*pow((k0*3), alpha - 1) - delta);
				 double r_d = loss*r;
				 double w = (1 - taup)*(1- tauw)*z*(1 - alpha)*pow((k0*3), alpha);
				 double w_d = loss*w;

 #pragma omp parallel for
                for (int ib = 0; ib < Nb; ib++){
                    double b0 = b_grid[ib];
                    
                    int shock_add[Ny];
					double default_decision_vec[Ny];
					evenSpace(shock_add, ik*Nb + ib, Nb*NK, Ny);
					getElements(default_decision, default_decision_vec, shock_add, Ny);
					price_mat[(iy*NK*Nb)+(ik*Nb)+ib] = (1.0 - ddot(default_decision_vec, prob_vec, Ny)) / (1 + rf);

                    double d_colMaxVal = -1000.0;
                    int d_colMaxId = 0;
                
                    double nd_rowMaxVal = -1000.0;
                    int nd_rowMaxId = 0;
                    int nd_global_colMaxId = 0;
                    
                    for (int ikp = 0; ikp < NK; ikp++){
                        double kp = k_grid[ikp];
                        
                        int shock_add2[Ny];
                        double vd_vec[Ny];
                        int shock_add1[Ny];
                        double vo_vec[Ny];
                            
                        evenSpace(shock_add2, ikp, NK, Ny);
                        getElements(vd, vd_vec, shock_add2, Ny);
                        double temp0 = ddot(vd_vec, prob_vec, Ny);                          

                        evenSpace(shock_add1, ikp*Nb + (Nb-1), NK*Nb, Ny);
                        getElements(vo, vo_vec, shock_add1, Ny);
                        double temp1 = ddot(vo_vec, prob_vec, Ny);                          
                        double c_d = (1 + r_d)*k0 + w_d*.33 - kp;
                        double current_d = c_d > 0 ? util(c_d) + theta*beta*temp0 + (1 - theta)*beta*temp1: -1000.0;

                        if (current_d > d_colMaxVal){
                            d_colMaxVal = current_d;
                            d_colMaxId = ikp;
                        }
                        
                        double nd_colMaxVal = -1000.0;
                        int nd_colMaxId = 0;
                        
                        for(int ibp = 0; ibp < Nb; ibp++){
                            double bp = b_grid[ibp];

                            evenSpace(shock_add1, ikp*Nb + ibp, NK*Nb, Ny);
                            getElements(vo, vo_vec, shock_add1, Ny);
                            double temp1 = ddot(vo_vec, prob_vec, Ny);        
                            double price = price_mat[(iy*Nb*NK)+(ikp*Nb)+ibp];
                            double c = (1 + r)*k0 + w*.33 - kp + b0 - bp*price;
                            double current_nd = c > 0 ? util(c) + beta*temp1 : -1000.0;
                      
                            if (current_nd > nd_colMaxVal){
                                nd_colMaxVal = current_nd;
                                nd_colMaxId = ibp;
                            }
                        }//ibp
                        
                        if (nd_colMaxVal > nd_rowMaxVal){
							nd_rowMaxVal = nd_colMaxVal;
							nd_rowMaxId = ikp;
							nd_global_colMaxId = nd_colMaxId;
						}                       
                    }

                    int ind_d = iy*NK + ik;
                    vd_temp[ind_d] = d_colMaxVal;
                    kopt_d[ind_d] = k_grid[d_colMaxId];
               
                    int ind_nd = iy*NK*Nb + ik*Nb + ib;
                    vnd_temp[ind_nd] = nd_rowMaxVal;
                    bopt[ind_nd] = b_grid[nd_global_colMaxId];
                    kopt[ind_nd] = k_grid[nd_rowMaxId];
                    copt[ind_nd] = (1 + r)*k0 + w*.33 - k_grid[nd_rowMaxId] + \
                            b0 - b_grid[nd_global_colMaxId]*price_mat[(iy*Nb*NK)+(nd_rowMaxId*Nb)+nd_global_colMaxId];   
                    
                    val_diff[ind_nd] = vnd_temp[ind_nd] - vd_temp[ind_d];
                    cdef[ind_nd] = (1 + r_d)*k0 + w_d*.33 - k_grid[d_colMaxId] ;
                    if(vnd_temp[ind_nd] >= vd_temp[ind_d]){
                        default_decision[ind_nd] = 0;
                        vo_temp[ind_nd] = vnd_temp[ind_nd]; 
                    }
                    else{
                        default_decision[ind_nd] = 1.0;
                        vo_temp[ind_nd] = vd_temp[ind_d];
                    }                 
                
                }//ib END PARALLEL             
            }//ik

        }//iy
        copyArray(vd_temp, vd, NK, Ny);
        crit = getNorm(vo, vo_temp, NK, Ny);

		iter++;
        double defaults = getSum(default_decision, Nb*NK, Ny);
        std:: cout << "iter: " << iter << " crit is " << crit  << "  default num " << defaults << std::endl;
        if(crit > 600)
            break;
        getWeighted(vo, vo_temp, NK*Nb, Ny, damp);

    }

    std::cout << "gov problem done" << std::endl;
    
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], vo);
    mxSetM(plhs[0], NK*Nb*Ny);
    mxSetN(plhs[0], 1);
    
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[1], vnd_temp);
    mxSetM(plhs[1], NK*Nb*Ny);
    mxSetN(plhs[1], 1);

    plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[2], vd);
    mxSetM(plhs[2], NK*Ny);
    mxSetN(plhs[2], 1);
    
    plhs[3] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[3], kopt);
    mxSetM(plhs[3], NK*Nb*Ny);
    mxSetN(plhs[3], 1);
     
    plhs[4] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[4], copt);
    mxSetM(plhs[4], NK*Nb*Ny);
    mxSetN(plhs[4], 1);
   
    plhs[5] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[5], bopt);
    mxSetM(plhs[5], NK*Nb*Ny);
    mxSetN(plhs[5], 1);
   
    plhs[6] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[6], default_decision);
    mxSetM(plhs[6], NK*Nb*Ny);
    mxSetN(plhs[6], 1);
    
    plhs[7] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[7], price_mat);
    mxSetM(plhs[7], NK*Nb*Ny);
    mxSetN(plhs[7], 1);
    
    plhs[8] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[8], kopt_d);
    mxSetM(plhs[8], NK*Ny);
    mxSetN(plhs[8], 1);

    plhs[9] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[9], cdef);
    mxSetM(plhs[9], NK*Ny*Nb);
    mxSetN(plhs[9], 1);
}