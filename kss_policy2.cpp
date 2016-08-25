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
		outarr[i] = inarr[row*Ncol+i];
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

//finds numerical index---------------------------------------------------------------------------------------------
int find_in(double inarr1[], int N, double val){
	int ind_low = 0;
	int ind_high = 1;
    int weight_low = 1.0;
	int weight_high = 0;

    if (inarr1[0] < val)
        return 0;
    
	for (int i = 0; i < N; i++){
		if (inarr1[i] > val){
			ind_high = i;
			ind_low = i - 1;
			weight_high = (val - inarr1[ind_low])/(inarr1[ind_high] - inarr1[ind_low]);
			weight_low = 1 - weight_high;
            if(weight_high > weight_low)
                return ind_high;
            else
                return ind_low;
			break;
		}
		if (i == N-1){
            ind_high = i;
            return ind_high;
		}
	}
}

//dots the betahat and created design vector include debt---------------------------------------------------------------------------------------------
double getKprime_nd(double bigK, double aggShock, double debt, int iy, double betahat[], int N_reg, int NK, double aggregateK[]){    
    double temp1 = log(bigK);
    double temp2 = log(aggShock);
    double temp3 = debt;
    double *X;
    X = (double *)mxMalloc(N_reg * sizeof(double));
    X[0] = 1;
    X[1] = temp1;
    X[2] = temp2;
    X[3] = temp1*temp1;
    X[4] = iy < 4 ? 1.0 : 0.0;
    X[5] = X[4]*temp1;
    X[6] = iy == 10 ? 1.0 : 0.0;
    X[7] = X[6]*temp1;
    X[8] = iy == 11 ? 1.0 : 0.0;
    X[9] =  X[8]*temp1;
    X[10] = iy == 12 ? 1.0 : 0.0;
    X[11] = X[10]*temp1;
    X[12] = iy == 6 ? 1.0 : 0.0;
    X[13] = X[12]*temp1;
    X[14] = iy == 7 ? 1.0 : 0.0;
    X[15] = X[14]*temp1;
    X[16] = iy == 8 ? 1.0 : 0.0;
    X[17] = X[16]*temp1;
    X[18] = iy == 9 ? 1.0 : 0.0;
    X[19] = X[18]*temp1;
    X[20] = iy == 13 ? 1.0 : 0.0;
    X[21] = X[20]*temp1;
    X[22] = debt;
    double bigKprime = exp(ddot(X, betahat, N_reg));
    bigKprime = std::min(bigKprime, aggregateK[NK-1]);
    bigKprime = std::max(bigKprime, aggregateK[0]);         
    return bigKprime;
}

//dots the betahat and created design vector ---------------------------------------------------------------------------------------------
double getKprime_d(double bigK, double aggShock, int iy, double betahat[], int N_reg, int NK, double aggregateK[]){    
    double temp1 = log(bigK);
    double temp2 = log(aggShock);
    double *X;
    X = (double *)mxMalloc(N_reg * sizeof(double));
    X[0] = 1;
    X[1] = temp1;
    X[2] = temp2;
    X[3] = temp1*temp1;
    X[4] = iy < 4 ? 1.0 : 0.0;
    X[5] = X[4]*temp1;
    X[6] = iy == 10 ? 1.0 : 0.0;
    X[7] = X[6]*temp1;
    X[8] = iy == 11 ? 1.0 : 0.0;
    X[9] =  X[8]*temp1;
    X[10] = iy == 12 ? 1.0 : 0.0;
    X[11] = X[10]*temp1;
    X[12] = iy == 6 ? 1.0 : 0.0;
    X[13] = X[12]*temp1;
    X[14] = iy == 7 ? 1.0 : 0.0;
    X[15] = X[14]*temp1;
    X[16] = iy == 8 ? 1.0 : 0.0;
    X[17] = X[16]*temp1;
    X[18] = iy == 9 ? 1.0 : 0.0;
    X[19] = X[18]*temp1;
    X[20] = iy == 13 ? 1.0 : 0.0;
    X[21] = X[20]*temp1;
    double bigKprime = exp(ddot(X, betahat, N_reg));
    bigKprime = std::min(bigKprime, aggregateK[NK-1]);
    bigKprime = std::max(bigKprime, aggregateK[0]);             
    return bigKprime;
}

//interpolation bilinear with valueslice between k (xdim) and b (ydim)---------------------------------------------------------------------------------------------
void bilinear_Val(double valslice[], double kprime, double bprime, double aggregateK[], double b_grid[], double Value[], \
        int prod_type, int ia, int NK, int Nb, int Ny, int Nk, int G){      
    int ind_low_k, ind_high_k, ind_low_b, ind_high_b;
    double weight_low_k, weight_high_k, weight_low_b, weight_high_b;
    
    find_between(aggregateK, NK, kprime, ind_low_k, ind_high_k, weight_low_k, weight_high_k);
    find_between(b_grid, Nb, bprime, ind_low_b, ind_high_b, weight_low_b, weight_high_b);
    
    double denom1 = (aggregateK[ind_high_k] - aggregateK[ind_low_k])*(b_grid[ind_high_b] - b_grid[ind_low_b]);
    double xdev2 = aggregateK[ind_high_k] - kprime;
    double xdev1 = kprime - aggregateK[ind_low_k];
    double ydev2 = b_grid[ind_high_b] - bprime;
    double ydev1 = bprime - b_grid[ind_low_b];
    double tempy[2] = {ydev2, ydev1};
    
#pragma omp parallel for 
    for (int ii = 0; ii < Ny; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int ind11 = (ind_low_b*NK*Ny*4*G*Nk) + (ind_low_k*Ny*4*G*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;
            long long int ind21 = (ind_low_b*NK*Ny*4*G*Nk) + (ind_high_k*Ny*4*G*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;           
            long long int ind12 = (ind_high_b*NK*Ny*4*G*Nk) + (ind_low_k*Ny*4*G*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;
            long long int ind22 = (ind_high_b*NK*Ny*4*G*Nk) + (ind_high_k*Ny*4*G*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;
            
            double temp1[2] = {xdev2*Value[ind11] + xdev1*Value[ind21], xdev2*Value[ind12] + xdev1*Value[ind22]};            
            double num = ddot(temp1, tempy, 2);
            valslice[ii*Nk + jj] = num/denom1;
        }
    }
}

//get default individual from vnd and vd iterpolation---------------------------------------------------------------------------------------------
void get_value_slice_nd(double val_slice[], int iK, int ib, double value[], int prod_type, int ia, int NK, int Nb, int Ny, int Nk, int G){      
#pragma omp parallel for 
    for (int ii = 0; ii < Ny; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int add_nd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;
            val_slice[ii*Nk+jj] = value[add_nd];
        }//jj
    }//ii
}

//get default individual from default and vd iterpolation---------------------------------------------------------------------------------------------
void get_value_slice_d(double val_slice[], int iK,  double value[], int prod_type, int ia, int NK, int Ny, int Nk, int G){      
#pragma omp parallel for 
    for (int ii = 0; ii < Ny; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int add_d = (iK*Ny*G*4*Nk) + (ii*4*G*Nk) + (ia*4*Nk) + (prod_type*Nk) + jj;
            val_slice[ii*Nk+jj] = value[add_d];
        }//jj
    }//ii
}

//interpolation 1d with valueslice---------------------------------------------------------------------------------------------
void interp_Val(double valslice[], double kprime, double aggregateK[], double Value[], int prod_type, int ia, int  NK, int  Ny, int  Nk, int  G){      
    int ind_low, ind_high;
    double weight_low, weight_high;
    
    find_between(aggregateK, NK, kprime, ind_low, ind_high, weight_low, weight_high);
    
#pragma omp parallel for 
    for (int ii = 0; ii < Ny; ii++){
        for (int jj = 0; jj < Nk; jj++){
            long long int lowadd = (ind_low*Ny*G*4*Nk) + (ii*G*4*Nk) + (ia*Nk*4) + (prod_type*Nk) +  jj;
            long long int highadd = (ind_high*Ny*G*4*Nk) + (ii*G*4*Nk) + (ia*Nk*4) + (prod_type*Nk) + jj;
            valslice[ii*Nk+jj] = weight_low*Value[lowadd] + weight_high*Value[highadd];
        }//jj
    }//ii
}

//interpolation bprime using bdec with states b0, k0 and iy---------------------------------------------------------------------------------------------
double interp_bprime(double kprime, double aggregateK[], double bprime_mat[], int iy, int ib, int NK, int Nb){        
    int ind_low_k, ind_high_k;
    double weight_low_k, weight_high_k;
    
    find_between(aggregateK, NK, kprime, ind_low_k, ind_high_k, weight_low_k, weight_high_k);
    double bprime = weight_low_k*bprime_mat[(iy*NK*Nb)+(ind_low_k*Nb)+ib] + 
            weight_high_k*bprime_mat[(iy*NK*Nb)+(ind_high_k*Nb)+ib];
    return bprime;
}

//copies array ---------------------------------------------------------------------------------------------
void copyArray(double inarr[], double outarr[], int Nrow, int Ncol){
	for (int i = 0; i < Nrow*Ncol; i++)
		outarr[i] = inarr[i];
}


//sets valueo to max of nd and d---------------------------------------------------------------------------------------------
void setValue_o(double Value_nd[], double Value_d[], double Value_o[], int NK, int Ny, int Nb, int Nk, int G, int endAge, int startAge){
    for (int ia = endAge; ia > startAge; ia --){
        for (int ib = 0; ib < Nb; ib++){
            for(int iK = 0; iK < NK; iK++){
                for(int iy = 0; iy < Ny; iy++){
#pragma omp parallel for                    
                    for(int k = 0; k < Nk; k++){
                        long long int currentadd_nd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;	
                        long long int currentadd_d = (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;	
                        
                        Value_o[currentadd_nd] = Value_nd[currentadd_nd] >= Value_d[currentadd_d] ? Value_nd[currentadd_nd] : Value_d[currentadd_d];
                        Value_o[currentadd_nd+Nk] = Value_nd[currentadd_nd+Nk] >= Value_d[currentadd_d+Nk] ? Value_nd[currentadd_nd+Nk] : Value_d[currentadd_d+Nk];
                        Value_o[currentadd_nd+2*Nk] = Value_nd[currentadd_nd+2*Nk] >= Value_d[currentadd_d+2*Nk] ? Value_nd[currentadd_nd+2*Nk] : Value_d[currentadd_d+2*Nk];
                        Value_o[currentadd_nd+3*Nk] = Value_nd[currentadd_nd+3*Nk] >= Value_d[currentadd_d+3*Nk] ? Value_nd[currentadd_nd+3*Nk] : Value_d[currentadd_d+3*Nk];                           
                    }//k
                } //iy
            }//iK
        }//ib       
    } //ia
}


//MAIN-------------------------------------------------------------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){    
    //Parameters
    const double eps[2] = {.68, 1.32};
    const double probz[4] = {.9608, .0392, .0392, .9608};
	const double skill[2] = {.57, 1.43};
    
//     Needed Inputs 
    double *k_grid, *aggregateK, *prob, *Z, *mass, *surv_prob, *prod_vec, \
            *par, *b_grid, *b_net, *default_decision, *KprimeMat_nd, *KprimeMat_d, *bprimeMat;
    k_grid = mxGetPr(prhs[0]);
    prob = mxGetPr(prhs[1]);
    Z = mxGetPr(prhs[2]);
    mass = mxGetPr(prhs[3]);
    surv_prob = mxGetPr(prhs[4]);
    prod_vec = mxGetPr(prhs[5]);
    aggregateK = mxGetPr(prhs[6]);
    par = mxGetPr(prhs[7]);
    b_grid = mxGetPr(prhs[8]);
    bprimeMat = mxGetPr(prhs[9]);
    b_net = mxGetPr(prhs[10]);
    default_decision = mxGetPr(prhs[11]);
    KprimeMat_nd = mxGetPr(prhs[12]);
    KprimeMat_d = mxGetPr(prhs[13]);
    
    const double alfa = par[0];
    const double delta = par[1];
    const double tauk = par[2];
    const double tauw = par[3];
    const double taup = par[4];
    const double beta = par[5];
    const double loss = par[6];
    const double rf = par[7];
    const double theta = par[8];
    const double fullout = par[9];
    
    const int G = 30;
    const int R = 20;
    const int Nk = 200;
    const int Ny = 15;
    const int NK = 52;
    const int Nb = 25;
    const int N_reg_d = 22;
    const int N_reg_nd = 23;

    
    double mid_k = k_grid[(int)(Nk/2)];

    double *Value_nd, *kopt_nd, *copt_nd, *bopt_nd, *Value_d, *kopt_d, *copt_d, *Value_o;
    mwSize ELEMENTS_nd = NK*Ny*G*4*Nk*Nb;
    mwSize COL = 1;
    Value_nd = (double *)mxMalloc(ELEMENTS_nd * sizeof(double));
    kopt_nd = (double *)mxMalloc(ELEMENTS_nd * sizeof(double));
    copt_nd = (double *)mxMalloc(ELEMENTS_nd * sizeof(double));
    bopt_nd = (double *)mxMalloc(ELEMENTS_nd * sizeof(double));
    Value_o = (double *)mxMalloc(ELEMENTS_nd * sizeof(double));
    
    mwSize ELEMENTS_d = NK*Ny*G*4*Nk;
    Value_d = (double *)mxMalloc(ELEMENTS_d * sizeof(double));
    kopt_d = (double *)mxMalloc(ELEMENTS_d * sizeof(double));
    copt_d = (double *)mxMalloc(ELEMENTS_d * sizeof(double));

    double working = 0.0, retired = 0.0;
    
    for (int i = 0; i < R - 1; i++)
		working += mass[2*i]*2;
    for (int i = R - 1; i < G; i++)
		retired += mass[2*i]*2;

    double prob_vec[Ny];
    
    double ValueSlice0_o[Nk*Ny];
    double ValueSlice1_o[Nk*Ny];
    double ValueSlice2_o[Nk*Ny];
    double ValueSlice3_o[Nk*Ny];
    
    double ValueSlice0_d[Nk*Ny];
    double ValueSlice1_d[Nk*Ny];
    double ValueSlice2_d[Nk*Ny];
    double ValueSlice3_d[Nk*Ny];
    
    double r, w, receipts, pen, bond, r_d, w_d, receipts_d, pen_d, bigK, aggShock;
    double dd, bigKprime_nd, bprime, bigKprime_d;
    int kprime_ind, bprime_ind;
    
////Last age---------------------------------------------------------------------------------------------------------
    for (int ib = 0; ib < Nb; ib++){
        for (int iK = 0; iK < NK; iK++){
            for (int iy = 0; iy < Ny; iy++){
                bigK = aggregateK[iK];
                aggShock = Z[iy];

                r = (1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta);
                receipts = working*.33*taup*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                pen = (1 - tauk)*receipts/retired;                                              
                bond = default_decision[(iy*NK*Nb)+(iK*Nb)+ib]  == 0 ? b_net[(iy*NK*Nb)+(iK*Nb)+ib] : 0.0;              
                
#pragma omp parallel for
                for (int k = 0; k < Nk; k++){  
                    double c = (1.0 + r)*k_grid[k] + pen + bond;
                    long long int currentadd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((G - 1)*4*Nk) + k;
                    double util = log(c);
                    Value_nd[currentadd] = util;
                    Value_nd[currentadd+Nk] = util;
                    Value_nd[currentadd+Nk*2] = util;
                    Value_nd[currentadd+Nk*3] = util;

                    copt_nd[currentadd] = c;
                    copt_nd[currentadd+Nk] = c;
                    copt_nd[currentadd+Nk*2] = c;
                    copt_nd[currentadd+Nk*3] = c;
                }//k END PARALLEL
                
//----------------------------------------------------------------------------------------------------
                //default region
                if (ib == Nb - 1){
                    double r_d = loss*r;
                    double pen_d = loss*pen;
#pragma omp parallel for
                    for (int k = 0; k < Nk; k++){  
                        double c_d = (1.0 + r_d)*k_grid[k] + pen_d;
                        long long int currentadd = (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((G - 1)*4*Nk) + k;
                        
                        double util_d = log(c_d);
                        Value_d[currentadd] = util_d;
                        Value_d[currentadd+Nk] = util_d;
                        Value_d[currentadd+Nk*2] = util_d;
                        Value_d[currentadd+Nk*3] = util_d;

                        copt_d[currentadd] = c_d;
                        copt_d[currentadd+Nk] = c_d;
                        copt_d[currentadd+Nk*2] = c_d;
                        copt_d[currentadd+Nk*3] = c_d;
                    }//k END PARALLEL                    
                } //end default region
//----------------------------------------------------------------------------------------------------
            }//iy
        }//iK
    }//ib
    
    setValue_o(Value_nd, Value_d, Value_o, NK, Ny, Nb, Nk, G, G, G-1);


//Retirement------------------------------------------------------------------------------------------------
    for (int ia = G - 1; ia > R - 1; ia--){
        for (int ib = 0; ib < Nb; ib++){
            for (int iK = 0; iK < NK; iK++){
                double bigK = aggregateK[iK];
                for (int iy = 0; iy < Ny; iy++){
                    aggShock = Z[iy];                                    
                    
                    getRow(prob, prob_vec, Ny, Ny, iy);
                    dd = default_decision[(iy*NK*Nb)+(iK*Nb)+ib];

                    if (dd == 0){
                        bond = b_net[(iy*NK*Nb)+ib];
                        r = (1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta);
                        receipts = working*.33*taup*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                        
                        bigKprime_nd = KprimeMat_nd[(iy*NK*Nb)+(iK*Nb)+ib];
                        bprime = bprimeMat[(iy*NK*Nb)+(iK*Nb)+ib];
                        
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);
                        bprime_ind = find_in(b_grid, Nb, bprime);
    
                        get_value_slice_nd(ValueSlice0_o, kprime_ind, bprime_ind, Value_o, 0, ia, NK, Nb, Ny, Nk, G);                                  
                    }
                    else{
                        bond = 0.0;
                        r = loss*((1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta));
                        receipts = loss*working*.33*taup*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                        
                        bigKprime_nd = KprimeMat_d[(iy*NK)+iK];
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);
                        
                        get_value_slice_d(ValueSlice0_o, kprime_ind, Value_d, 0, ia, NK, Ny, Nk, G);                                  
                    }
                    pen = (1 - tauk)*receipts/retired;

#pragma omp parallel for 
                    for (int k = 0; k < Nk; k++){
                        double k0 = k_grid[k];
                        double currentMax = -1000.0;
                        double currentCons = 0.0;
                        int currentMax_ind = 0;

                        for (int kp = 0; kp < Nk; kp++){
                            int ind_kp = k0 > mid_k ? Nk - 1 - kp : kp;
                            double kprime = k_grid[ind_kp];
                            double c = (1.0 + r)*k0 + pen - kprime + bond;

                            if (c < 0) continue;
                            int tempind[Ny];
                            double temparr_nd[Ny];
                            evenSpace(tempind, ind_kp, Nk, Ny);
                            getElements(ValueSlice0_o, temparr_nd, tempind, Ny);

                            double temp = log(c) + beta*surv_prob[ia-1]*ddot(temparr_nd, prob_vec, Ny);
                            if (temp >= currentMax){
                                currentMax = temp;
                                currentMax_ind = ind_kp;
                                currentCons = c;
                            }
                            else if(kp > 25)
                                break;
                        } //kp

                        long long int currentadd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;		
                        Value_nd[currentadd] = currentMax;
                        kopt_nd[currentadd] = k_grid[currentMax_ind];
                        copt_nd[currentadd] = currentCons;
                        Value_nd[currentadd+Nk] = currentMax;
                        kopt_nd[currentadd+Nk] = k_grid[currentMax_ind];
                        copt_nd[currentadd+Nk] = currentCons;
                        Value_nd[currentadd+Nk*2] = currentMax;
                        kopt_nd[currentadd+Nk*2] = k_grid[currentMax_ind];
                        copt_nd[currentadd+Nk*2] = currentCons;
                        Value_nd[currentadd+Nk*3] = currentMax;
                        kopt_nd[currentadd+Nk*3] = k_grid[currentMax_ind];
                        copt_nd[currentadd+Nk*3] = currentCons;
                    } //k  END PARALLEL
                    
//-----------------------------------------------------------------------------------------------------------------------------
                    //Default Region
                    if(ib == Nb - 1){
                        r_d = loss*((1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta));
                        receipts_d = loss*working*.33*taup*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                        pen_d = (1 - tauk)*receipts_d/retired;
                        
                        bigKprime_d = KprimeMat_d[(iy*NK)+iK];
                        kprime_ind = find_in(aggregateK, NK, bigKprime_d);
                        get_value_slice_d(ValueSlice0_d, kprime_ind, Value_d, 0, ia, NK, Ny, Nk, G);                                  

                        bigKprime_nd = KprimeMat_nd[(iy*NK*Nb)+(iK*Nb)+ib];                        
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);
                        bprime_ind = ib;
                        get_value_slice_nd(ValueSlice0_o, kprime_ind, bprime_ind, Value_o, 0, ia, NK, Nb, Ny, Nk, G);       
                                                
#pragma omp parallel for 
                        for (int k = 0; k < Nk; k++){
                            double k0 = k_grid[k];
                            double currentMax = -1000.0;
                            double currentCons = 0.0;
                            int currentMax_ind = 0;

                            for (int kp = 0; kp < Nk; kp++){
                                int ind_kp = k0 > mid_k ? Nk - 1 - kp : kp;
                                double kprime = k_grid[ind_kp];
                                double c_d = (1.0 + r_d)*k0 + pen_d - kprime;

                                if (c_d < 0) continue;
                                int tempind[Ny];
                                double temparr_d[Ny];
                                double temparr_nd[Ny];
                                evenSpace(tempind, ind_kp, Nk, Ny);
                                getElements(ValueSlice0_d, temparr_d, tempind, Ny);
                                getElements(ValueSlice0_o, temparr_nd, tempind, Ny);

                                double temp = log(c_d) + beta*surv_prob[ia-1]*((1 - theta)*ddot(temparr_d, prob_vec, Ny) + theta*ddot(temparr_nd, prob_vec, Ny));
                                if (temp >= currentMax){
                                    currentMax = temp;
                                    currentMax_ind = ind_kp;
                                    currentCons = c_d;
                                }
                                else if(kp > 25)
                                    break;
                            } //kp

                            int currentadd = (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;		

                            Value_d[currentadd] = currentMax;
                            kopt_d[currentadd] = k_grid[currentMax_ind];
                            copt_d[currentadd] = currentCons;
                            Value_d[currentadd+Nk] = currentMax;
                            kopt_d[currentadd+Nk] = k_grid[currentMax_ind];
                            copt_d[currentadd+Nk] = currentCons;
                            Value_d[currentadd+Nk*2] = currentMax;
                            kopt_d[currentadd+Nk*2] = k_grid[currentMax_ind];
                            copt_d[currentadd+Nk*2] = currentCons;
                            Value_d[currentadd+Nk*3] = currentMax;
                            kopt_d[currentadd+Nk*3] = k_grid[currentMax_ind];
                            copt_d[currentadd+Nk*3] = currentCons;
                        } //k  END PARALLEL
                    }  // end default region          
//-----------------------------------------------------------------------------------------------------------------------------              
                } //iy                
           } //iK
        } //ib
        setValue_o(Value_nd, Value_d, Value_o, NK, Ny, Nb, Nk, G, G-1, R-1);
    }  //age

    std::cout << "retired policy done " << std::endl;

    
//WorkingAge-----------------------------------------------------------------------------------------------------------------------
    for (int ia = R - 1; ia > 0; ia--){
        for(int ib = 0; ib < Nb; ib++){
            for (int iK = 0; iK < NK; iK++){
                bigK = aggregateK[iK];
                for (int iy = 0; iy < Ny; iy++){
                    aggShock = Z[iy];
                    
                    getRow(prob, prob_vec, Ny, Ny, iy);
                    dd = default_decision[(iy*NK*Nb)+(iK*Nb)+ib];

                    if (dd == 0){
                        bond = b_net[(iy*NK*Nb)+ib];
                        r = (1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta);
                        w = (1 - taup)*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                        
                        bigKprime_nd = KprimeMat_nd[(iy*NK*Nb)+(iK*Nb)+ib];
                        bprime = bprimeMat[(iy*NK*Nb)+(iK*Nb)+ib];                       
                        
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);
                        bprime_ind = find_in(b_grid, Nb, bprime);

                        get_value_slice_nd(ValueSlice0_o, kprime_ind, bprime_ind, Value_o, 0, ia, NK, Nb, Ny, Nk, G);      
                        get_value_slice_nd(ValueSlice1_o, kprime_ind, bprime_ind, Value_o, 1, ia, NK, Nb, Ny, Nk, G);      
                        get_value_slice_nd(ValueSlice2_o, kprime_ind, bprime_ind, Value_o, 2, ia, NK, Nb, Ny, Nk, G);      
                        get_value_slice_nd(ValueSlice3_o, kprime_ind, bprime_ind, Value_o, 3, ia, NK, Nb, Ny, Nk, G);      
                    }
                    else{
                        bond = 0;
                        r = loss*((1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta));
                        w = loss*(1 - taup)*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                        
                        bigKprime_nd = KprimeMat_d[(iy*NK)+iK];
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);

                        get_value_slice_d(ValueSlice0_o, kprime_ind, Value_d, 0, ia, NK, Ny, Nk, G);          
                        get_value_slice_d(ValueSlice1_o, kprime_ind, Value_d, 1, ia, NK, Ny, Nk, G);          
                        get_value_slice_d(ValueSlice2_o, kprime_ind, Value_d, 2, ia, NK, Ny, Nk, G);          
                        get_value_slice_d(ValueSlice3_o, kprime_ind, Value_d, 3, ia, NK, Ny, Nk, G);          
                    }
                    
#pragma omp parallel for 
                    for (int k = 0; k < Nk; k++){
                        double k0 = k_grid[k];
                        double currentMax0 = -1000.0, currentMax1 = -1000.0, currentMax2 = -1000.0, currentMax3 = -1000.0;
                        double currentCon0 = 0.0, currentCon1 = 0.0, currentCon2 = 0.0, currentCon3 = 0.0;

                        int currentMax_ind0 = 0, currentMax_ind1 = 0, currentMax_ind2 = 0, currentMax_ind3 = 0;
                        int found0 = 0, found1 = 0, found2 = 0, found3 = 0;

                        for (int kp = 0; kp < Nk; kp++){
                            int tempind[Ny];
                            double temparr0_nd[Ny];
                            double temparr1_nd[Ny];
                            double temparr2_nd[Ny];
                            double temparr3_nd[Ny];

                            int ind_kp = k0 > mid_k ? Nk - 1 - kp : kp;
                            double kprime = k_grid[ind_kp];

                            evenSpace(tempind, ind_kp, Nk, Ny);
                            getElements(ValueSlice0_o, temparr0_nd, tempind, Ny);
                            getElements(ValueSlice1_o, temparr1_nd, tempind, Ny);
                            getElements(ValueSlice2_o, temparr2_nd, tempind, Ny);
                            getElements(ValueSlice3_o, temparr3_nd, tempind, Ny);

                            //low skill low shock
                            if (found0 == 0){         
                                double c0 = (1.0 + r)*k0 - kprime + w*eps[0]*skill[0]*.33*(1 - tauw)*prod_vec[ia-1] + bond;      
                                double temp0 = c0 > 0 ? log(c0) + beta*surv_prob[ia-1]*(probz[0]*ddot(temparr0_nd, prob_vec, Ny) + \
                                    probz[1]*ddot(temparr1_nd, prob_vec, Ny)) : -1000;
                                if (temp0 >= currentMax0){
                                    currentMax0 = temp0;
                                    currentMax_ind0 = ind_kp;
                                    currentCon0 = c0;
                                }
                                else if(kp > 25)
                                    found0 = 1;
                            }

                            //low skill high shock
                            if(found1 == 0){
                                double c1 = (1.0 + r)*k0 - kprime + w*eps[1]*skill[0]*.33*(1 - tauw)*prod_vec[ia-1] + bond;
                                double temp1 = c1 > 0 ? log(c1) + beta*surv_prob[ia-1]*(probz[2]*ddot(temparr0_nd, prob_vec, Ny) + \
                                    probz[3]*ddot(temparr1_nd, prob_vec, Ny)) : -1000;
                                if (temp1 >= currentMax1){
                                    currentMax1 = temp1;
                                    currentMax_ind1 = ind_kp;
                                    currentCon1 = c1;
                                }
                                else if(kp > 25)
                                    found1 = 1;
                            }

                            //high skill low shock
                            if (found2 == 0){
                                double c2 = (1.0 + r)*k0 - kprime + w*eps[0]*skill[1]*.33*(1 - tauw)*prod_vec[ia-1] + bond;
                                double temp2 = c2 > 0 ? log(c2) + beta*surv_prob[ia-1]*(probz[0]*ddot(temparr2_nd, prob_vec, Ny) + \
                                    probz[1]*ddot(temparr3_nd, prob_vec, Ny)) : -1000;
                                if (temp2 >= currentMax2){
                                    currentMax2 = temp2;
                                    currentMax_ind2 = ind_kp;
                                    currentCon2 = c2;
                                }
                               else if(kp > 25)
                                    found2 = 1;                  
                            }

                            //high skill high shock
                            if (found3 == 0){
                                double c3 = (1.0 + r)*k0 - kprime + w*eps[1]*skill[1]*.33*(1 - tauw)*prod_vec[ia-1] + bond;
                                double temp3 = c3 > 0 ? log(c3) + beta*surv_prob[ia-1]*(probz[2]*ddot(temparr2_nd, prob_vec, Ny) + \
                                    probz[3]*ddot(temparr3_nd, prob_vec, Ny)) : -1000;
                                if (temp3 >= currentMax3){
                                    currentMax3 = temp3;
                                    currentMax_ind3 = ind_kp;
                                    currentCon3 = c3;
                                }
                                else if(kp > 25)
                                    found3 = 1;
                            }

                            if ((found0 == 1) & (found1 == 1) & (found2 == 1) & (found3 == 1))
                                break;
                        }//kp

                        int currentadd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;
                        Value_nd[currentadd] = currentMax0;
                        kopt_nd[currentadd] = k_grid[currentMax_ind0];
                        copt_nd[currentadd] = currentCon0;
                        Value_nd[currentadd+Nk] = currentMax1;
                        kopt_nd[currentadd+Nk] = k_grid[currentMax_ind1];
                        copt_nd[currentadd+Nk] = currentCon1;
                        Value_nd[currentadd+2*Nk] = currentMax2;
                        kopt_nd[currentadd+2*Nk] = k_grid[currentMax_ind2];
                        copt_nd[currentadd+2*Nk] = currentCon2;
                        Value_nd[currentadd+3*Nk] = currentMax3;
                        kopt_nd[currentadd+3*Nk] = k_grid[currentMax_ind3];
                        copt_nd[currentadd+3*Nk] = currentCon3;
                    }//k  END PARALLEL

//-----------------------------------------------------------------------------------------------------------------------------
                    //Default Region
                    if(ib == Nb - 1){
                        r_d = loss*((1 - tauk)*(aggShock*alfa*pow(bigK*3, (alfa - 1)) - delta));
                        w_d = loss*(1 - taup)*aggShock*(1 - alfa)*pow(bigK*3, alfa);
                       
                        bigKprime_d = KprimeMat_d[(iy*NK)+iK];
                        kprime_ind = find_in(aggregateK, NK, bigKprime_d);
                        
                        get_value_slice_d(ValueSlice0_d, kprime_ind, Value_d, 0, ia, NK, Ny, Nk, G);                                  
                        get_value_slice_d(ValueSlice1_d, kprime_ind, Value_d, 1, ia, NK, Ny, Nk, G);                                  
                        get_value_slice_d(ValueSlice2_d, kprime_ind, Value_d, 2, ia, NK, Ny, Nk, G);                                  
                        get_value_slice_d(ValueSlice3_d, kprime_ind, Value_d, 3, ia, NK, Ny, Nk, G);                                  

                        bigKprime_nd = KprimeMat_nd[(iy*NK*Nb)+(iK*Nb)+ib];                        
                        kprime_ind = find_in(aggregateK, NK, bigKprime_nd);
                        bprime_ind = ib;
                        
                        get_value_slice_nd(ValueSlice0_o, kprime_ind, bprime_ind, Value_o, 0, ia, NK, Nb, Ny, Nk, G);       
                        get_value_slice_nd(ValueSlice1_o, kprime_ind, bprime_ind, Value_o, 1, ia, NK, Nb, Ny, Nk, G);       
                        get_value_slice_nd(ValueSlice2_o, kprime_ind, bprime_ind, Value_o, 2, ia, NK, Nb, Ny, Nk, G);       
                        get_value_slice_nd(ValueSlice3_o, kprime_ind, bprime_ind, Value_o, 3, ia, NK, Nb, Ny, Nk, G);       
                                                                  
#pragma omp parallel for 
                        for (int k = 0; k < Nk; k++){
                            double k0 = k_grid[k];
                            double currentMax0 = -1000.0, currentMax1 = -1000.0, currentMax2 = -1000.0, currentMax3 = -1000.0;
                            double currentCon0 = 0.0, currentCon1 = 0.0, currentCon2 = 0.0, currentCon3 = 0.0;

                            int currentMax_ind0 = 0, currentMax_ind1 = 0, currentMax_ind2 = 0, currentMax_ind3 = 0;
                            int found0 = 0, found1 = 0, found2 = 0, found3 = 0;

                            for (int kp = 0; kp < Nk; kp++){
                                int tempind[Ny];
                                double temparr0_nd[Ny];
                                double temparr1_nd[Ny];
                                double temparr2_nd[Ny];
                                double temparr3_nd[Ny];

                                double temparr0_d[Ny];
                                double temparr1_d[Ny];
                                double temparr2_d[Ny];
                                double temparr3_d[Ny];
                                
                                int ind_kp = k0 > mid_k ? Nk - 1 - kp : kp;
                                double kprime = k_grid[ind_kp];

                                evenSpace(tempind, ind_kp, Nk, Ny);
                                getElements(ValueSlice0_o, temparr0_nd, tempind, Ny);
                                getElements(ValueSlice1_o, temparr1_nd, tempind, Ny);
                                getElements(ValueSlice2_o, temparr2_nd, tempind, Ny);
                                getElements(ValueSlice3_o, temparr3_nd, tempind, Ny);
                                
                                getElements(ValueSlice0_d, temparr0_d, tempind, Ny);
                                getElements(ValueSlice1_d, temparr1_d, tempind, Ny);
                                getElements(ValueSlice2_d, temparr2_d, tempind, Ny);
                                getElements(ValueSlice3_d, temparr3_d, tempind, Ny);
                                
                                //low skill low shock
                                if (found0 == 0){         
                                    double c0_d = (1.0 + r_d)*k0 - kprime + w_d*eps[0]*skill[0]*.33*(1 - tauw)*prod_vec[ia-1];      
                                    double val_weight0 = (1 - theta)*(probz[0]*ddot(temparr0_d, prob_vec, Ny) + probz[1]*ddot(temparr1_d, prob_vec, Ny)) + 
                                            theta*(probz[0]*ddot(temparr0_nd, prob_vec, Ny) + probz[1]*ddot(temparr1_nd, prob_vec, Ny));                                    
                                    double temp0 = c0_d > 0 ? log(c0_d) + beta*surv_prob[ia-1]*val_weight0 : -1000;
                                    if (temp0 >= currentMax0){
                                        currentMax0 = temp0;
                                        currentMax_ind0 = ind_kp;
                                        currentCon0 = c0_d;
                                    }
                                    else if(kp > 25)
                                        found0 = 1;
                                }

                                //low skill high shock
                                if(found1 == 0){
                                    double c1_d = (1.0 + r_d)*k0 - kprime + w_d*eps[1]*skill[0]*.33*(1 - tauw)*prod_vec[ia-1];
                                    double val_weight1 = (1 - theta)*(probz[2]*ddot(temparr0_d, prob_vec, Ny) + probz[3]*ddot(temparr1_d, prob_vec, Ny)) + 
                                            theta*(probz[2]*ddot(temparr0_nd, prob_vec, Ny) + probz[3]*ddot(temparr1_nd, prob_vec, Ny));                                       
                                    double temp1 = c1_d > 0 ? log(c1_d) + beta*surv_prob[ia-1]*val_weight1 : -1000;
                                    if (temp1 >= currentMax1){  
                                        currentMax1 = temp1;
                                        currentMax_ind1 = ind_kp;
                                        currentCon1 = c1_d;
                                    }
                                    else if(kp > 25)
                                        found1 = 1;
                                }

                                //high skill low shock
                                if (found2 == 0){
                                    double c2_d = (1.0 + r_d)*k0 - kprime + w_d*eps[0]*skill[1]*.33*(1 - tauw)*prod_vec[ia-1];
                                    double val_weight2 = (1 - theta)*(probz[0]*ddot(temparr2_d, prob_vec, Ny) + probz[1]*ddot(temparr3_d, prob_vec, Ny)) + 
                                        theta*(probz[0]*ddot(temparr2_nd, prob_vec, Ny) + probz[1]*ddot(temparr3_nd, prob_vec, Ny));                                                                                
                                    double temp2 = c2_d > 0 ? log(c2_d) + beta*surv_prob[ia-1]*val_weight2 : -1000;
                                    if (temp2 >= currentMax2){
                                        currentMax2 = temp2;
                                        currentMax_ind2 = ind_kp;
                                        currentCon2 = c2_d;
                                    }
                                   else if(kp > 25)
                                        found2 = 1;                  
                                }

                                //high skill high shock
                                if (found3 == 0){
                                    double c3_d = (1.0 + r_d)*k0 - kprime + w_d*eps[1]*skill[1]*.33*(1 - tauw)*prod_vec[ia-1];
                                    double val_weight3 = (1 - theta)*(probz[2]*ddot(temparr2_d, prob_vec, Ny) + probz[3]*ddot(temparr3_d, prob_vec, Ny)) + 
                                            theta*(probz[2]*ddot(temparr2_nd, prob_vec, Ny) + probz[3]*ddot(temparr3_nd, prob_vec, Ny));     
                                    double temp3 = c3_d > 0 ? log(c3_d) + beta*surv_prob[ia-1]*val_weight3 : -1000;
                                    if (temp3 >= currentMax3){
                                        currentMax3 = temp3;
                                        currentMax_ind3 = ind_kp;
                                        currentCon3 = c3_d;
                                    }
                                    else if(kp > 25)
                                        found3 = 1;
                                }

                                if ((found0 == 1) & (found1 == 1) &  (found2 == 1) & (found3 ==1))
                                    break;
                            }//kp

                            int currentadd = (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;
                            
                            Value_d[currentadd] = currentMax0;
                            kopt_d[currentadd] = k_grid[currentMax_ind0];
                            copt_d[currentadd] = currentCon0;
                            Value_d[currentadd+Nk] = currentMax1;
                            kopt_d[currentadd+Nk] = k_grid[currentMax_ind1];
                            copt_d[currentadd+Nk] = currentCon1;
                            Value_d[currentadd+2*Nk] = currentMax2;
                            kopt_d[currentadd+2*Nk] = k_grid[currentMax_ind2];
                            copt_d[currentadd+2*Nk] = currentCon2;
                            Value_d[currentadd+3*Nk] = currentMax3;
                            kopt_d[currentadd+3*Nk] = k_grid[currentMax_ind3];
                            copt_d[currentadd+3*Nk] = currentCon3;
                        } //k  END PARALLEL
                    }  // end default region    
//-----------------------------------------------------------------------------------------------------------------------------
                }  //iy
            }  //ik
        }  //ib
        setValue_o(Value_nd, Value_d, Value_o, NK, Ny, Nb, Nk, G, R-1, 0);
    } //age
    
    std::cout << "working policy done " << std::endl;

    
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], kopt_nd);
    mxSetM(plhs[0], ELEMENTS_nd);
    mxSetN(plhs[0], COL);
    
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[1], kopt_d);
    mxSetM(plhs[1], ELEMENTS_d);
    mxSetN(plhs[1], COL);
    
    plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[2], Value_nd);
    mxSetM(plhs[2], ELEMENTS_nd);
    mxSetN(plhs[2], COL);
    
    plhs[3] = mxCreateNumericMatrix(0, 0,  mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[3], Value_d);
    mxSetM(plhs[3], ELEMENTS_d);
    mxSetN(plhs[3], COL);

    if (fullout > 0)
    {
        plhs[4] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mxSetPr(plhs[4], Value_o);
        mxSetM(plhs[4], ELEMENTS_nd);
        mxSetN(plhs[4], COL);
        
        plhs[5] = mxCreateNumericMatrix(0, 0,  mxDOUBLE_CLASS, mxREAL);
        mxSetPr(plhs[5], copt_nd);
        mxSetM(plhs[5], ELEMENTS_nd);
        mxSetN(plhs[5], COL);
    
        plhs[6] = mxCreateNumericMatrix(0, 0,  mxDOUBLE_CLASS, mxREAL);
        mxSetPr(plhs[6], copt_d);
        mxSetM(plhs[6], ELEMENTS_d);
        mxSetN(plhs[6], COL);
    }   
}