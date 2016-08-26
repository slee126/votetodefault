% Generate plots for endowment and production economy

clear; close all; clc;

need_guess = 0;
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
rf = .05;
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
iter = 1;
par = pargov;
par(6) = 1.02;
par(10) = 0;
grid_slice = [0; 0; 0];

aggregateK = transpose(.25:.01:.76);
r = (1-tauk)*(alpha*(3*aggregateK).^(alpha-1) - delta);

%----------------------------------------------------------------------------------------------------------------------
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

[vo1, vnd, vd1, kopt, copt, bopt, default_decision_avg, price_mat, kopt_def, c_def] = ...
        gov_prob_guess(aggregateK, prob2, Z', pargov, b_grid, vo, vd);

sd = .02;
rho = .9;
Ny = 15;
[prob, s] = markovappr(rho, sd, 3, Ny);
Z = exp(s);
prob2 = reshape(prob', Ny*Ny, 1);

load guess;
[vo1, vnd, vd1, kopt, copt, bopt, default_decision_avg, price_mat, kopt_def, c_def] = ...
        gov_prob_guess(aggregateK, prob2, Z', pargov, b_grid, vo, vd);
        
default_mat = reshape(default_decision_avg, Nb, NbigK*Ny);
default_slice = reshape(default_mat(13, :), NbigK, Ny);
save iter1 default_mat default_slice price_mat bopt s


sd = .02;
rho = .8;
Ny = 15;
[prob, s] = markovappr(rho, sd, 3, Ny);
Z = exp(s);
prob2 = reshape(prob', Ny*Ny, 1);

load guess;
[vo1, vnd, vd1, kopt, copt, bopt, default_decision_avg, price_mat, kopt_def, c_def] = ...
        gov_prob_guess(aggregateK, prob2, Z', pargov, b_grid, vo, vd);
        
default_mat = reshape(default_decision_avg, Nb, NbigK*Ny);
default_slice = reshape(default_mat(13, :), NbigK, Ny);
save iter2 default_mat default_slice price_mat bopt s

sd = .02;
rho = 0;
Ny = 15;
[prob, s] = markovappr(rho, sd, 3, Ny);
Z = exp(s);
prob2 = reshape(prob', Ny*Ny, 1);

load guess;
[vo1, vnd, vd1, kopt, copt, bopt, default_decision_avg, price_mat, kopt_def, c_def] = ...
        gov_prob_guess(aggregateK, prob2, Z', pargov, b_grid, vo, vd);
        
default_mat = reshape(default_decision_avg, Nb, NbigK*Ny);
default_slice = reshape(default_mat(13, :), NbigK, Ny);
save iter3 default_mat default_slice price_mat bopt s


%-------------------------------------------------------------------------
figure(1); clf;
set(1, 'defaulttextinterpreter', 'latex');

load iter1;
subplot(331)
default_slice = reshape(default_mat(14, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
ylabel('Capital');
title('$\rho$=.9, $\sigma$=.02, bond state=-.022');

subplot(334)
default_slice = reshape(default_mat(13, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
ylabel('Capital');
title('$\rho$=.9, $\sigma$=.02, bond state=-.025');

subplot(337)
default_slice = reshape(default_mat(12, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
xlabel('Z');
ylabel('Capital');
title('$\rho$=.9, $\sigma$=.02, bond state=-.027');

load iter2;
subplot(332)
default_slice = reshape(default_mat(14, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
title('$\rho$=.8, $\sigma$=.02, bond state=-.022');

subplot(335)
default_slice = reshape(default_mat(13, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
title('$\rho$=.8, $\sigma$=.02, bond state=-.025');

subplot(338)
default_slice = reshape(default_mat(12, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
xlabel('Z');
title('$\rho$=.8, $\sigma$=.02, bond state=-.027');

load iter3;
default_slice = reshape(default_mat(14, :), NbigK, Ny);
subplot(333)
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
title('$\rho$=0, $\sigma$=.02, bond state=-.022');

subplot(336)
default_slice = reshape(default_mat(13, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
title('$\rho$=0, $\sigma$=.02, bond state=-.025');

subplot(339)
default_slice = reshape(default_mat(12, :), NbigK, Ny);
imagesc(exp(s), aggregateK, ~default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
xlabel('Z');
title('$\rho$=0, $\sigma$=.02, bond state=-.027');

figtitle('Default Set Production Economy', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/defaultset_Prod.pdf
%--------------------------------------------------------------------------
% 
% 
%--------------------------------------------------------------------------
figure(2); clf ;
set(2, 'defaulttextinterpreter', 'latex');

load iter1;
subplot(331)
pmat = reshape(price_mat, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
ylabel('Bond Price');
title('$\rho$=.9, $\sigma$=.02, bond state=-.022');

subplot(334)
b_ind = 13;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
ylabel('Bond Price');
title('$\rho$=.9, $\sigma$=.02, bond state=-.025');

subplot(337)
b_ind = 12;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
ylabel('Bond Price');
title('$\rho$=.9, $\sigma$=.02, bond state=-.027');

load iter2;
subplot(332)
pmat = reshape(price_mat, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.022');

subplot(335)
b_ind = 13;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.025');

subplot(338)
b_ind = 12;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
xlabel('Capital');
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.027');

load iter3;
subplot(333)
pmat = reshape(price_mat, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'all Z'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=0, $\sigma$=.02, bond state=-.022');

subplot(336)
b_ind = 13;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'all Z'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=0, $\sigma$=.02, bond state=-.025');

subplot(339)
b_ind = 12;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'all Z'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
title('$\rho$=0, $\sigma$=.02, bond state=-.027');

figtitle('Bond Price Production Economy', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/bondPrice_prod.pdf
%--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
figure(3); clf ;
set(3, 'defaulttextinterpreter', 'latex');

load iter1;
subplot(331)
bmat = reshape(bopt, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
ylabel('Bond Choice');
title('$\rho$=.9, $\sigma$=.02, bond state=-.022');

subplot(334)
b_ind = 13;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
ylabel('Bond Choice');
title('$\rho$=.9, $\sigma$=.02, bond state=-.025');

subplot(337)
b_ind = 12;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
ylabel('Bond Choice');
title('$\rho$=.9, $\sigma$=.02, bond state=-.027');

load iter2;
subplot(332)
bmat = reshape(bopt, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.022');

subplot(335)
b_ind = 13;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.025');

subplot(338)
b_ind = 12;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
title('$\rho$=.8, $\sigma$ =.02, bond state=-.027');

load iter3;
subplot(333)
bmat = reshape(bopt, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.022');

subplot(336)
b_ind = 13;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('$\rho$=.8, $\sigma$=.02, bond state=-.025');

subplot(339)
b_ind = 12;
y_ind = 3;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 9;1
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
y_ind = 12;
plot(aggregateK, bmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
grid on; axis tight;
legend({'Z=.9064', 'Z=1.02', 'Z=1.08'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
title('$\rho$=0, $\sigma$=.02, bond state=-.027');

figtitle('Bond Chocie Production Economy', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/bondChoice_prod.pdf
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
figure(4); clf ;
set(4, 'defaulttextinterpreter', 'latex');

load iter1;
subplot(121)
pmat = reshape(price_mat, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 3;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
r = (1-tauk)*(exp(s(y_ind))*.33*(aggregateK*3).^(.33-1) - delta);
bp = 1./(1+r);
plot(aggregateK, bp);
grid on; axis tight;
legend({'bond price', 'implied net rental return'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
ylabel('Bond Price');
title('bond state=-.022, Z=.90');

subplot(122)
pmat = reshape(price_mat, 25, NbigK*Ny)';
b_ind = 14;
y_ind = 12;
plot(aggregateK, pmat((y_ind-1)*NbigK+1:y_ind*NbigK, b_ind));  hold on;
r = (1-tauk)*(exp(s(y_ind))*.33*(aggregateK*3).^(.33-1) - delta);
bp = 1./(1+r);
plot(aggregateK, bp);
grid on; axis tight;
legend({'bond price', 'implied net rental return'}, 'location', 'best', 'Color', 'none', 'box', 'off');
xlabel('Capital');
ylabel('Bond Price');
title('bond state=-.022, Z=1.08');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
figtitle('Implied Rental Return Bond Price', 'FontSize', 16);
export_fig plots/implied.pdf;











