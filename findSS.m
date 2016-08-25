
clear; close all; clc;
tic;

eps = [0.68; 1.32];
probz = [0.98 0.02; 0.02 0.98]^2;
skill = [.57; 1.43];

beta = 1.02;
Nk = 1000;
kmax = 4;
kmin = 0;

kgrid = linspace(kmin, kmax, Nk)';
kp = kgrid;

finer_grid = linspace(kmin, kmax, 2*Nk)';
step = finer_grid(3) - finer_grid(2);

G = 30;
R = 20;
load prodpop;
surv = zeros(G, 1);
for i = 1:G
   surv(i) = surv_prob((i-1)*3 + 1); 
end

ages = (21:2:60);
prod_profile = interp1(20:74, prod_vec, ages)';

mass = ones(G, 1);
for i = 1:G-1
    mass(i+1) = mass(i)*surv(i); 
end
mass = mass/sum(mass);
mass = repmat(mass, 1, 2)/2;

working = sum(mass(1:R-1, 1))*2;
retired = sum(mass(R:end, 1))*2;

alfa = .33;
bigK = .45;
delta = .1;
taup = .2;
tauk = .3;
tauw = .3;
loss = .96;

bigKold = bigK;
Value = zeros(4*G, Nk);
kopt = zeros(4*G, Nk);

iter = 0;
crit = 1;


while(crit > 1e-2 && iter < 100 )
    % Last age -----------------------------------------------------------------
    r = loss*(1 - tauk)*(alfa*(bigK*3)^(alfa - 1) - delta);
    w = loss*(1 - taup)*(1 - alfa)*(bigK*3)^alfa;
    receipts = loss*working*.33*taup*(1 - alfa)*(bigK*3)^alfa;
    pen = (1 - tauk)*receipts/retired;
    
    %current indices            
    c = (1 + r)*kgrid + pen;
    Value(237:240, :) = repmat(log(c), 1, 4)'; 

    % -------------------------------------------------------------------------
    % 
    % 
    % Retirement-------------------------------------------------------------
    
    for age = G-1:-1:R        
        for k = 1:Nk
            k0 = kgrid(k);
            c = (1 + r)*k0 + pen - kp;
            ind = c<0;
            temp = log(c) + beta*surv(age)*Value((age*4+1), :)'; 
            temp(ind) = -1000;
            [Value((age-1)*4 + 1:(age-1)*4 + 4, k), ind] = max(temp);
            kopt((age-1)*4 + 1:(age-1)*4 + 4, k) = kp(ind);
        end %closes k (capital grid)
    end  % closes age
    
    % --------------------------------------------------------------------------
    % 
    % Working------------------------------------------------------------------
    for age = R-1:-1:1

        for k = 1:Nk
            k0 = kgrid(k);

            % low prod low shock
            idio_prob_vec = probz(1, :);
            idio_shock = eps(1);
            idio_skill = skill(1);
            c = (1 + r)*k0 + idio_skill*idio_shock*w*.33*(1 - tauw)*prod_profile(age) - kp;
            ind = c<0;
            temp1 = idio_prob_vec(1)*Value(age*4 + 1, :)'; 
            temp2 = idio_prob_vec(2)*Value(age*4 + 2, :)'; 

            vtemp = (temp1 + temp2);
            temp = log(c) + beta*surv(age)*vtemp; 
            temp(ind) = -1000;
            [Value((age-1)*4 + 1, k), ind] = max(temp);
            kopt((age-1)*4 + 1, k) = kp(ind);

            % low prod high shock
            idio_prob_vec = probz(2, :);
            idio_shock = eps(2);
            idio_skill = skill(1);
            c = (1 + r)*k0 + idio_skill*idio_shock*w*.33*(1 - tauw)*prod_profile(age) - kp;
            ind = c<0;
            temp1 = idio_prob_vec(1)*Value(age*4 + 1, :)'; 
            temp2 = idio_prob_vec(2)*Value(age*4 + 2, :)'; 

            vtemp = (temp1 + temp2);
            temp = log(c) + beta*surv(age)*vtemp; 
            temp(ind) = -1000;
            [Value((age-1)*4 + 2, k), ind] = max(temp);
            kopt((age-1)*4 + 2, k) = kp(ind);

            % high prod low shock
            idio_prob_vec = probz(1, :);
            idio_shock = eps(1);
            idio_skill = skill(2);
            c = (1 + r)*k0 + idio_skill*idio_shock*w*.33*(1 - tauw)*prod_profile(age) - kp;
            ind = c<0;
            temp1 = idio_prob_vec(1)*Value(age*4 + 3, :)'; 
            temp2 = idio_prob_vec(2)*Value(age*4 + 4, :)'; 

            vtemp = (temp1+temp2);
            temp = log(c) + beta*surv(age)*vtemp; 
            temp(ind) = -1000;
            [Value((age-1)*4 + 3, k), ind] = max(temp);
            kopt((age-1)*4 + 3, k) = kp(ind);
            
            % high prod high shock
            idio_prob_vec = probz(2, :);
            idio_shock = eps(2);
            idio_skill = skill(2);
            c = (1 + r)*k0 + idio_skill*idio_shock*w*.33*(1 - tauw)*prod_profile(age) - kp;
            ind = c<0;
            temp1 = idio_prob_vec(1)*Value(age*4 + 3, :)'; 
            temp2 = idio_prob_vec(2)*Value(age*4 + 4, :)'; 

            vtemp = (temp1+temp2);
            temp = log(c) + beta*surv(age)*vtemp; 
            temp(ind) = -1000;
            [Value((age-1)*4 + 4, k), ind] = max(temp);
            kopt((age-1)*4 + 4, k) = kp(ind);
        end 
    end %age 
   
    % --------------------------------------------------------------------------
    fprintf('policy done\n');

    %--------------------------------------------------------------------------
    
   
    %distributional grid 2 x finer than policy
    cap_dist_low_low = zeros(2*Nk, G);
    cap_dist_low_high = zeros(2*Nk, G);
    cap_dist_high_low = zeros(2*Nk, G);
    cap_dist_high_high = zeros(2*Nk, G);

    cap_dist_low_low(1, 1) = mass(1, 1)/2;
    cap_dist_low_high(1, 1) = mass(1, 1)/2;
    cap_dist_high_low(1, 1) = mass(1, 2)/2;
    cap_dist_high_high(1, 1) = mass(1, 2)/2;

    for i = 1:G-1
        k0_q_vec_low_low = finer_grid(cap_dist_low_low(:, i) > 0);
        k0_q_vec_low_high = finer_grid(cap_dist_low_high(:, i) > 0);
        k0_q_vec_high_low = finer_grid(cap_dist_high_low(:, i) > 0);
        k0_q_vec_high_high = finer_grid(cap_dist_high_high(:, i) > 0);

        mass_low_low = cap_dist_low_low(cap_dist_low_low(:, i) > 0, i)*surv(i);
        mass_low_high = cap_dist_low_high(cap_dist_low_high(:, i) > 0, i)*surv(i);
        mass_high_low = cap_dist_high_low(cap_dist_high_low(:, i) > 0, i)*surv(i);
        mass_high_high = cap_dist_high_high(cap_dist_high_high(:, i) > 0, i)*surv(i);

        %mass for age 
        tot_mass = sum(mass_low_low) + sum(mass_low_high) + sum(mass_high_low) + sum(mass_high_high);

        kp_q_vec_low_low = interp1(kgrid, kopt((i-1)*4 + 1, :), k0_q_vec_low_low);
        kp_q_vec_low_high = interp1(kgrid, kopt((i-1)*4 + 2, :), k0_q_vec_low_high);
        kp_q_vec_high_low = interp1(kgrid, kopt((i-1)*4 + 3, :), k0_q_vec_high_low);
        kp_q_vec_high_high = interp1(kgrid, kopt((i-1)*4 + 4, :), k0_q_vec_high_high);

        %low type and low shock
        for j = 1:length(k0_q_vec_low_low)
            kp_q_low_low = min(kp_q_vec_low_low(j), finer_grid(end));
            ind0_low_low = floor(kp_q_low_low/step) + 1;                  
            w1_low_low = (kp_q_low_low - finer_grid(ind0_low_low))/step;

            ind_h1 = ind0_low_low;
            ind_h2 = min(ind0_low_low + 1, 2*Nk);

            cap_dist_low_low(ind_h2, i+1) =  cap_dist_low_low(ind_h2, i+1) +  w1_low_low*mass_low_low(j)*probz(1, 1);
            cap_dist_low_low(ind_h1, i+1) = cap_dist_low_low(ind_h1, i+1) +  (1 - w1_low_low)*mass_low_low(j)*probz(1, 1);

            cap_dist_low_high(ind_h2, i+1) = cap_dist_low_high(ind_h2, i+1) + w1_low_low*mass_low_low(j)*probz(1, 2);
            cap_dist_low_high(ind_h1, i+1) = cap_dist_low_high(ind_h1, i+1) + (1 - w1_low_low)*mass_low_low(j)*probz(1, 2);
        end

        %low type and high shock
        for j = 1:length(k0_q_vec_low_high)
            kp_q_low_high = min(kp_q_vec_low_high(j), finer_grid(end));
            ind0_low_high = floor(kp_q_low_high/step) + 1;                            

            w1_low_high = (kp_q_low_high - finer_grid(ind0_low_high))/step;
            ind_h1 = ind0_low_high;
            ind_h2 = min(ind0_low_high+1, 2*Nk);

            cap_dist_low_high(ind_h2, i+1) =  cap_dist_low_high(ind_h2, i+1) + w1_low_high*mass_low_high(j)*probz(2, 2);
            cap_dist_low_high(ind_h1, i+1) = cap_dist_low_high(ind_h1, i+1) + (1 - w1_low_high)*mass_low_high(j)*probz(2, 2);

            cap_dist_low_low(ind_h2, i+1) = cap_dist_low_low(ind_h2, i+1) + w1_low_high*mass_low_high(j)*probz(2, 1);
            cap_dist_low_low(ind_h1, i+1) = cap_dist_low_low(ind_h1, i+1) + (1 - w1_low_high)*mass_low_high(j)*probz(2, 1);
        end

        %high type and low shock
        for j = 1:length(k0_q_vec_high_low)
            kp_q_high_low = min(kp_q_vec_high_low(j), finer_grid(end));

            ind0_high_low = floor(kp_q_high_low/step) + 1;             
            w1_high_low = mod(kp_q_high_low, step)/step; 

            ind_h1 = ind0_high_low;
            ind_h2 = min(ind0_high_low+1, 2*Nk);

            cap_dist_high_low(ind_h2, i+1) =   cap_dist_high_low(ind_h2, i+1) + w1_high_low*mass_high_low(j)*probz(1, 1);
            cap_dist_high_low(ind_h1, i+1) = cap_dist_high_low(ind_h1, i+1) + (1 - w1_high_low)*mass_high_low(j)*probz(1, 1);

            cap_dist_high_high(ind_h2, i+1) = cap_dist_high_high(ind_h2, i+1) + w1_high_low*mass_high_low(j)*probz(1, 2);
            cap_dist_high_high(ind_h1, i+1) = cap_dist_high_high(ind_h1, i+1) + (1 - w1_high_low)*mass_high_low(j)*probz(1, 2);
        end

        %high type and high shock
        for j = 1:length(k0_q_vec_high_high)
            kp_q_high_high = min(kp_q_vec_high_high(j), finer_grid(end));
            ind0_high_high = floor(kp_q_high_high/step) + 1;                  
            w1_high_high = mod(kp_q_high_high, step)/step;           

            ind_h1 = ind0_high_high;
            ind_h2 = min(ind0_high_high+1, 2*Nk);

            cap_dist_high_high(ind_h2, i+1) =  cap_dist_high_high(ind_h2, i+1)  + w1_high_high*mass_high_high(j)*probz(2, 2);
            cap_dist_high_high(ind_h1, i+1) = cap_dist_high_high(ind_h1, i+1) + (1 - w1_high_high)*mass_high_high(j)*probz(2, 2);

            cap_dist_high_low(ind_h2, i+1) = cap_dist_high_low(ind_h2, i+1) + w1_high_high*mass_high_high(j)*probz(2, 1);
            cap_dist_high_low(ind_h1, i+1) = cap_dist_high_low(ind_h1, i+1)  + (1 - w1_high_high)*mass_high_high(j)*probz(2, 1);
        end
    end
    all_ = cap_dist_high_high + cap_dist_high_low + cap_dist_low_high + cap_dist_low_low;
    bigKnew = sum(all_, 2)'*finer_grid;
    crit = abs(bigKnew-bigK);
    fprintf('bigK: %1.8f r: %1.8f crit: %1.8f iter: %d\n', bigKnew, r, crit, iter);
    iter = iter + 1;
    bigK = bigKnew*.1 + .9*bigK;
end
    
bigKold = bigK;
minK = .4*bigK;
maxK = 1.2*bigK;

aggregateK = transpose(minK:.01:maxK);
NbigK = length(aggregateK);

sd = .02;
rho = .8;
Ny = 15;
[prob, s] = markovappr(rho, sd, 3, Ny);
Z = exp(s);
    
kss = bigK;

N = 1000;
temp = randsample(6:NbigK-5, N, true);
bigK_series = aggregateK(temp);

Z_series = zeros(N, 1);
Z_series(1) = 8; 
for i = 2:N
    Z_series(i) = randsample(Ny, 1, true, prob(Z_series(i-1), :));
end

Nb = 25;
b_min = -.05;
b_max = 0.0;
b_grid = linspace(b_min, b_max, Nb)';
temp = randsample(Nb, N, true);
B_series = b_grid(temp);

load betahat;

save mexPar;




