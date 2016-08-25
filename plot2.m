% load run1;
% 
default_mat = reshape(default_mat, 25, 15*52);
figure(1); clf;
set(1, 'defaulttextinterpreter', 'latex');
default_slice = reshape(default_mat(18, :), NbigK, Ny);
imagesc(exp(s), aggregateK, default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
xlabel('Technology Shock');
ylabel('Aggregate Capital');
title('Default set for bond state=-.0125', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig defaultset_KS.pdf

price_mat = reshape(price_mat, 25, 15*52);
figure(2); clf;
set(2, 'defaulttextinterpreter', 'latex');
price_slice = reshape(price_mat(18, :), NbigK, Ny)';
plot(Z, price_slice(:, 10)); hold on;
plot(Z, price_slice(:, 20)); hold on;
plot(Z, price_slice(:, 30)); hold on;
plot(Z, price_slice(:, 40)); hold on;
plot(Z, price_slice(:, 50)); hold on;
xlabel('Technology Shock');
ylabel('Bond Price');
grid on; axis tight;
legend({'K=.34', 'K=.44', 'K=.54'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('Prices for bond state=-.0125', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig bondPrice_KS.pdf


b_mat = reshape(bopt, 25, 15*52);
figure(3); clf;
set(3, 'defaulttextinterpreter', 'latex');
b_slice = reshape(b_mat(14, :), NbigK, Ny)';
subplot(131)
plot(Z, b_slice(:, 10));
grid on;
xlabel('Technology Shock');
ylabel('Bond Decision');
title('aggregateK=.34'); 
subplot(132)
plot(Z, b_slice(:, 15));
grid on;
xlabel('Technology Shock');
ylabel('Bond Decision');
title('aggregateK=.39'); 
subplot(133)
plot(Z, b_slice(:, 20)); hold on;
grid on;
xlabel('Technology Shock');
ylabel('Bond Decision');
title('aggregateK=.44'); 
figtitle('bond policy for bond state=-.0125', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig bondPolicy_KS.pdf

% 
% 
Value_nd_new = reshape(Value_nd_new, 200, 15*25*52*120);
Value_d_new = reshape(Value_d_new, 200, 15*52*120);

% 
% 
% iy = 3;
% iK = 41;
% ib = 19;
ia = 30;
% 
% currentadd_nd = (ib*NK*Ny*G*4*Nk) + (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;	
% currentadd_d = (iK*Ny*G*4*Nk) + (iy*G*4*Nk) + ((ia - 1)*4*Nk) + k;	
%                         
% 
% Value_d_new(currentadd_d:currentadd_d+20)
% 
% [Value_d_new(currentadd_d:currentadd_d+20), Value_nd_new(currentadd_nd:currentadd_nd+20)]
% 

counter = 1;
figure(4); clf; 
set(4, 'defaulttextinterpreter', 'latex');

ia = 30;
ikset = [12; 16; 20; 27];
for i = 1:4
    iK = ikset(i);
    counter = 1;
    label_ = {};
    subplot(2,2,i)
    for ib = 17:17
        for iy = 0:2:14
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);

            surplus = value_nd_slice - value_d_slice;
            if (sum(surplus < 0))
                iyset{i}(counter) = iy;
                inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
                fprintf('surplus sum  %i ib %i iK %i iy %i outflow %1.6f\n', sum(surplus < 0), ib, iK, iy, inflow_);
                plot(kgrid, surplus); hold on;  
                temp = sprintf('Z=%1.3f\n',Z(iy+1));
                label_{counter} = temp;
                counter = counter + 1;
                plot(kgrid, surplus); hold on;  
            end

        end
    end

    title(sprintf('K=%1.2f', aggregateK(iK+1)), 'FontSize', 16);
    grid on;
    legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
    xlabel('$$k_t$$');
    ylabel('surplus');
end

figtitle('$$V^{R, G} - V^{D, G}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value30_KS.pdf


figure(5); clf;
set(5, 'defaulttextinterpreter', 'latex');
ia = 29;
for h = 1:4
    iK = ikset(h);
    subplot(2,2,h)
    counter = 1;
    clear label;
    for i = 1:length(iyset{h})
        iy = iyset{h}(i);
        currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
        currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
        value_nd_slice = Value_nd_new(:, currentadd_nd);
        value_d_slice = Value_d_new(:, currentadd_d);
        surplus = value_nd_slice - value_d_slice;
        plot(kgrid, surplus); hold on;
        temp = sprintf('Z=%1.3f\n',Z(iy+1));
        label{counter} = temp;
        counter = counter + 1;
    end
    grid on;   
    title(sprintf('K=%1.2f', aggregateK(iK+1)), 'FontSize', 16);
    legend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
    ylabel('surplus');
end
figtitle('$$V^{R, G-1} - V^{D, G-1}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value29_KS.pdf


figure(6); clf;
set(6, 'defaulttextinterpreter', 'latex');
ia = 21;
for h = 1:4
    iK = ikset(h);
    subplot(2,2,h)
    counter = 1;
    clear label;
    for i = 1:length(iyset{h})
        iy = iyset{h}(i);
        currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
        currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
        value_nd_slice = Value_nd_new(:, currentadd_nd);
        value_d_slice = Value_d_new(:, currentadd_d);
        surplus = value_nd_slice - value_d_slice;
        plot(kgrid, surplus); hold on;
        temp = sprintf('Z=%1.3f\n',Z(iy+1));
        label{counter} = temp;
        counter = counter + 1;
    end
    grid on;   
    title(sprintf('K=%1.2f', aggregateK(iK+1)), 'FontSize', 16);
    legend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
    xlabel('$$k_t$$');
    ylabel('surplus');
end
figtitle('$$V^{R, R} - V^{D, R}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value20_KS.pdf


titles_ = {'productive type low and low idiosyncratic shock', 'productive type low and high idiosyncratic shock', ...
    'productive type high and low idiosyncratic shock', 'productive type high and high idiosyncratic shock'};
figure(8); clf;
set(8, 'defaulttextinterpreter', 'latex');
ia = 19;
iK = 42;
for h = 1:4
    subplot(2,2,h)
    counter = 1;
    clear label;
    for i = 1:length(iyset{h})
        iy = iyset{h}(i);
        currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        value_nd_slice = Value_nd_new(:, currentadd_nd);
        value_d_slice = Value_d_new(:, currentadd_d);
        surplus = value_nd_slice - value_d_slice;
        plot(kgrid, surplus); hold on;
        temp = sprintf('Z=%1.3f\n',Z(iy+1));
        label{counter} = temp;
        counter = counter + 1;
    end
    grid on;   
    title(titles_{h}, 'FontSize', 16);
    legend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
    xlabel('$$k_t$$');
    ylabel('surplus');
end
figtitle('$$V^{R, R-1} - V^{D, R-1}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value19_KS.pdf


figure(9); clf;
set(9, 'defaulttextinterpreter', 'latex');
ia = 10;
iK = 42;
for h = 1:4
    subplot(2,2,h)
    counter = 1;
    clear label;
    for i = 1:length(iyset{h})
        iy = iyset{h}(i);
        currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        value_nd_slice = Value_nd_new(:, currentadd_nd);
        value_d_slice = Value_d_new(:, currentadd_d);
        surplus = value_nd_slice - value_d_slice;
        plot(kgrid, surplus); hold on;
        temp = sprintf('Z=%1.3f\n',Z(iy+1));
        label{counter} = temp;
        counter = counter + 1;
    end
    grid on;   
    title(titles_{h}, 'FontSize', 16);
    legend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
    xlabel('$$k_t$$');
    ylabel('surplus');
end
figtitle('$$V^{R, 10} - V^{D, 10}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value10_KS.pdf



figure(10); clf;
set(10, 'defaulttextinterpreter', 'latex');
ia = 5;
iK = 42;
for h = 1:4
    subplot(2,2,h)
    counter = 1;
    clear label;
    for i = 1:length(iyset{h})
        iy = iyset{h}(i);
        currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
        value_nd_slice = Value_nd_new(:, currentadd_nd);
        value_d_slice = Value_d_new(:, currentadd_d);
        surplus = value_nd_slice - value_d_slice;
        plot(kgrid, surplus); hold on;
        temp = sprintf('Z=%1.3f\n',Z(iy+1));
        label{counter} = temp;
        counter = counter + 1;
    end
    grid on;   
    title(titles_{h}, 'FontSize', 16);
    legend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
    xlabel('$$k_t$$');
    ylabel('surplus');
end
figtitle('$$V^{R, 5} - V^{D, 5}$$', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig Value5_KS.pdf


mass_by_age = mass(1:30, 1)*2;
figure(11); clf;
set(11, 'defaulttextinterpreter', 'latex');
ds0_new = reshape(ds0, 500, 30);
ds1_new = reshape(ds1, 500, 30);
ds2_new = reshape(ds3, 500, 30);
ds3_new = reshape(ds4, 500, 30);
dsall = ds0_new + ds1_new + ds2_new + ds3_new;
sub0 = sum(ds0_new);
sub1 = sum(ds1_new);
sub2 = sum(ds2_new);
sub3 = sum(ds3_new);
suball = sum(dsall);
temp = suball./mass_by_age';
% temp(end) = 1;
plot(1:30, temp);
xlabel('age');
ylabel('Proportion of Default Vote');
grid on;
title('Vote by Age for Low Capital Stock and Low Aggregate Shock', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig votedef_byage.pdf


nds0_new = reshape(nds0, 500, 30);
nds1_new = reshape(nds1, 500, 30);
nds2_new = reshape(nds2, 500, 30);
nds3_new = reshape(nds3, 500, 30);
ndsall = nds0_new + nds1_new + nds2_new + nds3_new;


csall = cs0 + cs1 + cs2 + cs3;
csall = reshape(csall, 500, 30);
sub_csall = sum(csall, 2);
suball2 = sum(dsall, 2);
figure(12); clf;
set(12, 'defaulttextinterpreter', 'latex');
plot(finer_grid, suball2./sub_csall);
grid on;
xlabel('$$k_t$$');
ylabel('Proportion of Default Vote');
title('Vote by wealth for Low Capital Stock and Low Aggregate Shock', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig votedef_bycapital.pdf

figure(12); clf;
set(12, 'defaulttextinterpreter', 'latex');
plot(1:30, mass_by_age);
grid on; axis tight;
title('Population density by age', 'FontSize', 16);
xlabel('age');
ylabel('Proportion of population');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig popdensity.pdf


temp = 0;
mu = sub_csall'*finer_grid;
for i = 1:500
   for j = 1:500
      temp = temp + sub_csall(i)*sub_csall(j)*abs(finer_grid(i) - finer_grid(j));
   end
end
G = temp/(2*mu)
figure(13); clf;
set(13, 'defaulttextinterpreter', 'latex');
plot(finer_grid, sub_csall);
grid on; axis tight;
title('Population density by wealth', 'FontSize', 16);
xlabel('$$k_t$$');
ylabel('Proportion of population');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig gini.pdf




figure(7);
set(7, 'defaulttextinterpreter', 'latex');
clf;
iy = 3;
iK = 10;
ib = 17;
subplot(221)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('Low Capital Low Shock');
legend({'No Default', 'Default'}, 'Location', 'best', 'Color', 'none', 'box', 'off');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
iy = 14;
iK = 10;
ib = 19;
subplot(222)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('Low Capital High Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
iy = 3;
iK = 40;
ib = 19;
subplot(223)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('High Capital Low Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
iy = 14;
iK = 40;
ib = 19;
subplot(224)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('High Capital High Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
figtitle('Aggregate Capital Stock Comparison', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig AggCapital.pdf



plot(Z, bigK_realized_d(20:52:end)); hold on;
plot(Z, bigK_realized_d(30:52:end)); hold on;
plot(Z, bigK_realized_d(40:52:end)); hold on;
plot(Z, bigK_realized_d(50:52:end)); hold on;
legend({'10', '20', '30', '40', '50'}, 'Location', 'best', 'box', 'off', 'Color', 'none');
% 
% 
% 
% inflow_new = reshape(inflow, 25, 15*52);
% figure(5); clf;
% set(5, 'defaulttextinterpreter', 'latex');
% 
% iy = 3;
% plot(aggregateK, inflow_new(19, (iy-1)*52+1:iy*52)); hold on;
% 
% % 
% % grid on;
% % clickableLegend(label, 'Location', 'best', 'Color', 'none', 'Box', 'off');
% 

% 
% 
% grid on;
% 
% 
% legend(leg, 'Location', 'best', 'box', 'off', 'Color', 'none');
% 
% leg = {};
% for i = 1:15
%     leg{i}  = num2str(2*i); 
% end
% 
% figure(4); clf;
% set(4, 'defaulttextinterpreter', 'latex');
% plot(kgrid, value_nd_slice); hold on;
% plot(kgrid, value_d_slice, '-r'); hold on;
% 
% working = sum(mass(1:R-1, 1))*2;
% ret = 1-working;
% 
%     
% [kopt_nd, kopt_d, Value_nd_new, Value_d_new] = kss_policy2test(kgrid, prob2, Z', mass2, surv, prod_profile, ...
%     aggregateK, par, b_grid, bopt, inflow, default_mat, bigK_realized_nd, bigK_realized_d);
% 
%     
%     
% 
