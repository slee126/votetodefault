% clear; close all; clc;

need_guess = 1;
compile = 0;

if compile == 1
    mex gov_prob_guess.cpp COMPFLAGS="/openmp -O2 -Wall $COMPFLAGS";
    mex gov_prob.cpp COMPFLAGS="/openmp -O2 -Wall $COMPFLAGS";
    mex kss_policy.cpp COMPFLAGS="/openmp -O2 -Wall $COMPFLAGS";
    mex kss_policy2.cpp COMPFLAGS="/openmp -O2 -Wall $COMPFLAGS";
    mex kss_dist_d.cpp COMPFLAGS="-O2 -Wall $COMPFLAGS";
    mex kss_dist_nd.cpp COMPFLAGS="-O2 -Wall $COMPFLAGS";
end


% gov problem -------------------------------------------------------------
load mexPar;
prob2 = reshape(prob', Ny*Ny, 1);

delta = .1;
alpha = .33;
beta = .95;
loss = .96;
rf = .04;
theta = .3;

pargov = zeros(8, 1);
pargov(1) = alpha;
pargov(2) = delta;
pargov(3) = .30;
pargov(4) = .3;
pargov(5) = .2;
pargov(6) = beta;
pargov(7) = loss;
pargov(8) = rf;
pargov(9) = theta;

mass2 = reshape(mass', G*2, 1);
Nk = 200;
NK = 52;
Ny = 15;
Nb = 25;
kmin = 0;
kmax = 2.5;
kgrid = linspace(kmin, kmax, Nk)';
finer_grid = linspace(kmin, kmax, 500)';

iter = 1;
par = pargov;
par(6) = 1.02;
par(10) = 0;
grid_slice = [0; 5; 17];

aggregateK = transpose(.25:.01:.76);
%----------------------------------------------------------------------------------------------------------------------
if (need_guess == 1)
    vd = zeros(NbigK*Ny, 1);
    for i = 1:Ny
       for j = 1:NbigK 
          temp = Z(i)*aggregateK(j)^alpha;
          vd((i-1)*Nk + j)  = log(temp);
       end
    end

    vo  = zeros(NbigK*Nb*Ny, 1);
    for k = 1:Ny
        for i = 1:NbigK
           for j = 1:Nb 
              temp = Z(k)*aggregateK(i)^alpha + b_grid(j);
              vo((k-1)*Nb*NbigK + (i-1)*Nb + j)  = log(temp);
           end
        end
    end
    
    tic;
    [vo1, vnd, vd1, kopt, copt, bopt, default_decision_avg, price_mat, kopt_def, c_def] = ...
            gov_prob_guess(aggregateK, prob2, Z', pargov, b_grid, vo, vd);
    toc;

    vo = vo1;
    vd = vd1;
    inflow = repmat(b_grid, NbigK*Ny, 1) - bopt.*price_mat;
    inflow(default_decision_avg == 1) = 0;
    save guess vo vd bopt inflow default_decision_avg;
     
    [kopt_nd, kopt_d, Value_nd, Value_d] = kss_policy(kgrid, prob2, Z', mass2, surv, betahat_d, prod_profile, ...
        aggregateK, par, b_grid, bopt, inflow, default_decision_avg, betahat_nd);
     
    [bigK_realized_d, countMat_d, cs0_def, cs1_def, cs2_def, cs3_def] = kss_dist_d(kgrid, mass2, surv, finer_grid, kopt_d, grid_slice);
     
    [bigK_realized_nd, countMat_nd, default_mat, ds0, ds1, ds3, ds4, ...
        nds0, nds1, nds2, nds3, cs0, cs1, cs2, cs3, vote_nd, vote_d] = kss_dist_nd(kgrid, mass2, surv, finer_grid, kopt_nd, Value_nd, Value_d, grid_slice);
    
    save -v7.3 initrun;
end


%---------------------------------------------------------------------------------------------------------------------------
% 
% Full Problem
% 
%---------------------------------------------------------------------------------------------------------------------------

crit = 1;
tol = 1e-1;
iter = 1;
% load initrun;
tic;
diary on;
while(iter < 10 && crit > tol)

    [vo1, vnd, vd1, kopt_nd_avg, copt_nd_avg, bopt, price_mat, kopt_def_avg, c_def_avg] = ...
        gov_prob(aggregateK, prob2, Z', pargov, b_grid, vo, vd, default_mat);
    
    vo = vo1;
    vd = vd1;
    inflow = repmat(b_grid, NbigK*Ny, 1) - bopt.*price_mat;
    inflow(default_mat == 1) = 0;
    
    [kopt_nd, kopt_d, Value_nd_new, Value_d_new] = kss_policy2(kgrid, prob2, Z', mass2, surv, prod_profile, ...
        aggregateK, par, b_grid, bopt, inflow, default_mat, bigK_realized_nd, bigK_realized_d);
             
    [bigK_realized_d, countMat_d, cs0_def, cs1_def, cs2_def, cs3_def] = kss_dist_d(kgrid, mass2, surv, finer_grid, kopt_d, grid_slice);
     
    [bigK_realized_nd, countMat_nd, default_mat, ds0, ds1, ds3, ds4, ...
        nds0, nds1, nds2, nds3, cs0, cs1, cs2, cs3, vote_nd, vote_d] = kss_dist_nd(kgrid, mass2, surv, finer_grid, kopt_nd, Value_nd_new, Value_d_new, grid_slice);
    
    crit = max(abs(Value_nd_new - Value_nd));
    iter = iter + 1;
    
    Value_nd = Value_nd_new;    
    fprintf('iter is %i crit is %1.6f defaultcases is %1.2f\n', iter, crit, sum(default_mat));
    toc;
    
%     if (mod(iter, 3) == 0)
%         save -v7.3 run1;
%     end   
end

save -v7.3 run1;
diary off;

grid_slice = [0; 35; 14]

