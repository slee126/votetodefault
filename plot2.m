load run1;
% % 
default_mat1 = reshape(default_mat, 25, 15*52);
figure(1); clf;
set(1, 'defaulttextinterpreter', 'latex');
default_slice = reshape(default_mat1(14, :), NbigK, Ny);
imagesc(exp(s), aggregateK, default_slice', 'AlphaData', .5);           
colormap(gray);                            
grid on;
xlabel('Technology Shock');
ylabel('Aggregate Capital');
title('Default set for bond state=-.023', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/defaultset_KS.pdf

price_mat1 = reshape(price_mat, 25, 15*52);
figure(2); clf;
set(2, 'defaulttextinterpreter', 'latex');
price_slice = reshape(price_mat1(14, :), NbigK, Ny)';
plot(Z, price_slice(:, 10)); hold on;
plot(Z, price_slice(:, 20)); hold on;
plot(Z, price_slice(:, 40)); hold on;
xlabel('Technology Shock');
ylabel('Bond Price');
grid on;
legend({'K=.34', 'K=.44', 'K=.64'}, 'location', 'best', 'Color', 'none', 'box', 'off');
title('Prices for bond state=-.023', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/bondPrice_KS.pdf


b_mat1 = reshape(bopt, 25, 15*52);
figure(3); clf;
set(3, 'defaulttextinterpreter', 'latex');
b_slice = reshape(b_mat1(14, :), NbigK, Ny)';
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
figtitle('bond policy for bond state=-.023', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/bondPolicy_KS.pdf


Value_nd_new = reshape(Value_nd_new, 200, 15*25*52*120);
Value_d_new = reshape(Value_d_new, 200, 15*52*120);

iy = 3;
iK = 41;
ib = 14;
ia = 30;

ibset = [10; 13];
ia = 30;
ikset1 = [30; 35; 40; 45];
for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(4); clf; 
    set(4, 'defaulttextinterpreter', 'latex');
    for i = 1:4
        iK = ikset1(i);
        counter = 1;
        label_ = {};
        subplot(2,2,i)

        for iy = 0:4:14
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);

            surplus = value_nd_slice - value_d_slice;
            iyset{i}(counter) = iy;
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            fprintf('surplus sum  %i ib %i iK %i iy %i outflow %1.6f\n', sum(surplus < 0), ib, iK, iy, inflow_);
            if(sum(surplus<0))
                plot(kgrid, surplus, 'LineWidth', 2); hold on;
            else
                plot(kgrid, surplus); hold on;
            end
            temp = sprintf('Z=%1.3f B=%1.3f', Z(iy+1), inflow_);
            label_{counter} = temp;
            counter = counter + 1;
        end

        title(sprintf('K=%1.2f D=-%1.2f', aggregateK(iK+1), b_grid(ib)), 'FontSize', 16);
        grid on;
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, G} - V^{D, G}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value30_KS_%i.pdf', ib);
    export_fig(temp);
end

for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(5); clf;
    set(5, 'defaulttextinterpreter', 'latex');
    ia = 29;
    for h = 1:4
        iK = ikset1(h);
        subplot(2,2,h)
        counter = 1;
        clear label_;
        for iy = 0:4:14
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);
            surplus = value_nd_slice - value_d_slice;
            if(sum(surplus<0))
                plot(kgrid, surplus, 'LineWidth', 2); hold on;
            else
                plot(kgrid, surplus); hold on;
            end
            temp = sprintf('Z=%1.3f B=%1.3f', Z(iy+1), inflow_);
            label_{counter} = temp;
            counter = counter + 1;
        end
        grid on;   
        title(sprintf('K=%1.2f', aggregateK(iK+1)), 'FontSize', 16);
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, G-1} - V^{D, G-1}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value29_KS_%i.pdf', ib);
    export_fig(temp);
end

for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(6); clf;
    set(6, 'defaulttextinterpreter', 'latex');
    ia = 21;
    for h = 1:4
        iK = ikset1(h);
        subplot(2,2,h)
        counter = 1;
        clear label_;
        for iy = 0:4:14
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + 1;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);
            surplus = value_nd_slice - value_d_slice;
            if(sum(surplus<0))
                plot(kgrid, surplus, 'LineWidth', 2); hold on;
            else
                plot(kgrid, surplus); hold on;
            end
            temp = sprintf('Z=%1.3f B=%1.3f', Z(iy+1), inflow_);
            label_{counter} = temp;
            counter = counter + 1;
        end
        grid on;   
        title(sprintf('K=%1.2f', aggregateK(iK+1)), 'FontSize', 16);
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, R+1} - V^{D, R+1}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value21_KS_%i.pdf', ib);
    export_fig(temp);
end

titles_ = {'productive type low and low idiosyncratic shock', 'productive type low and high idiosyncratic shock', ...
    'productive type high and low idiosyncratic shock', 'productive type high and high idiosyncratic shock'};
for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(8); clf;
    set(8, 'defaulttextinterpreter', 'latex');
    ia = 19;
    iK = 42;
    for h = 1:4
        subplot(2,2,h)
        counter = 1;
        clear label_;
        for iy = 0:4:14
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);
            surplus = value_nd_slice - value_d_slice;
            plot(kgrid, surplus); hold on;
            temp = sprintf('Z=%1.3f\n',Z(iy+1));
            label_{counter} = temp;
            counter = counter + 1;
        end
        grid on;   
        title(titles_{h}, 'FontSize', 16);
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, R-1} - V^{D, R-1}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value19_KS_%i.pdf', ib);
    export_fig(temp);
end

for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(9); clf;
    set(9, 'defaulttextinterpreter', 'latex');
    ia = 10;
    iK = 42;
    for h = 1:4
        subplot(2,2,h)
        counter = 1;
        clear label_;
        for iy = 0:4:14
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);
            surplus = value_nd_slice - value_d_slice;
            plot(kgrid, surplus); hold on;
            temp = sprintf('Z=%1.3f\n',Z(iy+1));
            label_{counter} = temp;
            counter = counter + 1;
        end
        grid on;   
        title(titles_{h}, 'FontSize', 16);
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, 10} - V^{D, 10}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value10_KS_%i.pdf', ib);
    export_fig(temp);
end

for ib_ind = 1:2
    ib = ibset(ib_ind);
    figure(10); clf;
    set(10, 'defaulttextinterpreter', 'latex');
    ia = 5;
    iK = 42;
    for h = 1:4
        subplot(2,2,h)
        counter = 1;
        clear label_;
        for iy = 0:4:14
            inflow_ = inflow(iy*NbigK*Nb + iK*Nb + ib + 1);
            currentadd_nd = (ib*NbigK*Ny*G*4) + (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            currentadd_d = (iK*Ny*G*4) + (iy*G*4) + (ia - 1)*4 + h;
            value_nd_slice = Value_nd_new(:, currentadd_nd);
            value_d_slice = Value_d_new(:, currentadd_d);
            surplus = value_nd_slice - value_d_slice;
            plot(kgrid, surplus); hold on;
            temp = sprintf('Z=%1.3f\n',Z(iy+1));
            label_{counter} = temp;
            counter = counter + 1;
        end
        grid on;   
        title(titles_{h}, 'FontSize', 16);
        legend(label_, 'Location', 'best', 'Color', 'none', 'Box', 'off');
        xlabel('$$k_t$$');
        ylabel('surplus');
    end
    temp = sprintf('$$V^{R, 5} - V^{D, 5}$$ D=%1.2f', b_grid(ib));
    figtitle(temp, 'FontSize', 16);
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
    temp = sprintf('plots/Value5_KS_%i.pdf', ib);
    export_fig(temp);
end

find(vote_d > .43 & vote_d < .6)

% 
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
export_fig plots/votedef_byage.pdf


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
export_fig plots/votedef_bycapital.pdf

figure(13); clf;
set(13, 'defaulttextinterpreter', 'latex');
plot(1:30, mass_by_age);
grid on; axis tight;
title('Population density by age', 'FontSize', 16);
xlabel('age');
ylabel('Proportion of population');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/popdensity.pdf
% 
% 
temp = 0;
mu = sub_csall'*finer_grid;
for i = 1:500
   for j = 1:500
      temp = temp + sub_csall(i)*sub_csall(j)*abs(finer_grid(i) - finer_grid(j));
   end
end
G = temp/(2*mu)
figure(14); clf;
set(14, 'defaulttextinterpreter', 'latex');
plot(finer_grid, sub_csall);
grid on; axis tight;
title('Population density by wealth', 'FontSize', 16);
xlabel('$$k_t$$');
ylabel('Proportion of population');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/gini.pdf
% 
% 
% 
% 
figure(7);
set(7, 'defaulttextinterpreter', 'latex');
clf;
iy = 3;
iK = 10;
ib = 13;
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
ib = 13;
subplot(222)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('Low Capital High Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
iy = 3;
iK = 40;
ib = 12;
subplot(223)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('High Capital Low Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
iy = 14;
iK = 40;
ib = 13;
subplot(224)
plot(b_grid, bigK_realized_nd(iy*NbigK*Nb + iK*Nb + 1:iy*NbigK*Nb + iK*Nb + Nb)); hold on;
plot([b_grid(1), b_grid(end)], [bigK_realized_d(iy*NbigK + iK) bigK_realized_d(iy*NbigK + iK)]);
grid on;
title('High Capital High Shock');
ylabel('$$K_{t+1}$$');
xlabel('Debt $$D_t$$');
figtitle('Aggregate Capital Stock Comparison', 'FontSize', 16);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1200 500]);
export_fig plots/AggCapital.pdf


