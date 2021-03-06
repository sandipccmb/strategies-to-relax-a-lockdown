% v4: Plot illustrative figures for main manuscript

clear all; %load Model_setup_India;

Setup_model_India2;

p.q = 1; r.q = 0;

R0 = 3;
r.beta = r.beta*R0;

hosp_capac = 57194;

% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.(gps.geo{1}).ad) = seed; init(i.S.(gps.geo{1}).ad) = init(i.S.(gps.geo{1}).ad) - seed;


% --- Baseline (no intervention)
tf = 1e3;
M0 = make_model4(p, r, i, s, gps, prm);

geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
[t,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
% Find timeseries of mortality
tmort   = sum(soln0(:,i.aux.mort),2);
t_ldown = find(tmort>10, 1, 'first');

cbeta_list = [0.25 0.6];

for li = 1:length(cbeta_list)
    fprintf('%0.5g ', li);
    
    % --- With lockdown
    r1 = r; p1 = p;
    r1.beta = r.beta*(1-cbeta_list(li));
    M1 = make_model4(p1, r1, i, s, gps, prm);
    
    geq = @(t,in) goveqs_scaleup(t, in, M0, M1, i, s, p1, r, agg, sel, t_ldown+[0 10]);
    [t, soln1] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
    
    % Find when the level of hospitalisations has come down to 10% of
    % capacity
    hosp = sum(soln1(:,s.H),2);
    inds = find(hosp>0.1*hosp_capac);
    t_relax = inds(end);
    
    % Relaxing lockdown
    r2 = r1; p2 = p1;
    %r2.beta = r.beta*(1-cbeta_list(li)/2);
    r2.beta = r.beta;
    r2.q = 1/4; p2.q = 1;
    M2 = make_model4(p2, r2, i, s, gps, prm);
    
    geq = @(t,in) goveqs_scaleup(t, in, M1, M2, i, s, p2, r, agg, sel, t_relax+[0 10]);
    [t, tmp] = ode15s(geq, [t_relax:1:tf], soln1(t_relax,:), odeset('Nonnegative',1:i.nx));
    soln2 = soln1; soln2(t_relax+1:end,:) = tmp;

    % Relaxing lockdown
    r3 = r2; p3 = p2;
    r3.beta = r.beta;
    M3 = make_model4(p3, r3, i, s, gps, prm);
    
    hosp = sum(soln2(:,s.H),2);
    inds = find(hosp>0.1*hosp_capac);
    t_relax2 = inds(end);
    geq = @(t,in) goveqs_scaleup(t, in, M2, M3, i, s, p3, r, agg, sel, t_relax2+[0 10]);
    [t, tmp] = ode15s(geq, [t_relax2:1:tf], soln2(t_relax2,:), odeset('Nonnegative',1:i.nx));
    soln3 = soln2; soln3(t_relax2+1:end,:) = tmp;
    
    
    allsol1(:,:,li) = soln1; 
    allsol2(:,:,li) = soln2; 
    allsol3(:,:,li) = soln3; 
end
fprintf('\n');

allsol = cat(3, soln0, allsol2);

% Multiplier to put cases per 1,000 population
mult = 1e3/sum(prm.N);


Ht = squeeze(sum(allsol(:,s.H,:),2));

ff=figure; fs = 14; lw = 2;

% --- Plot numbers needing hospitalisation
subplot(2,1,1);
hold on;
% plot(Ht(:,1)/1e3, '--', 'linewidth', lw);
plot(Ht(:,1:end)*mult, 'linewidth', lw);
line(xlim, hosp_capac*[1 1]*mult, 'linestyle', ':', 'Color', 0.2*[1 1 1], 'linewidth', lw);
xlim([0 700]);

xlabel('Days');
ylabel({'Persons needing hospitalisation','(per thousand popn)'});
legend('No lockdown','25% reduction in transmission','60% reduction in transmission','Hospital capacity');
set(gca,'fontsize',fs);

xl = xlim; yl = ylim;
h = fill([t_ldown, t_ldown, t_relax, t_relax],[yl(1) yl(2) yl(2) yl(1)],'r');
set(h,'facealpha',.1,'EdgeColor','none');
%set(0,'DefaultLegendAutoUpdate','off');

h = fill([t_relax, t_relax, xl(2), xl(2)],[yl(1) yl(2) yl(2) yl(1)],'g');
set(h,'facealpha',.1,'EdgeColor','none');

% --- Plot seroprevalence over time
subplot(2,1,2);
hold on;
Rt = squeeze(sum(allsol(:,s.R,:),2))/sum(prm.N(1,:));
plot(Rt, 'linewidth', lw);
set(gca,'fontsize',fs);
xlim([0 700]);
set(ff,'Position',[680   351   462   627]);

ll = line(xlim, (1 - 1/R0)*[1 1], 'linestyle', ':', 'Color', 0.2*[1 1 1], 'linewidth', lw);

legend(ll,'Approx herd immunity threshold');
ylabel('Proportion of population immune');
xlabel('Days');

xl = xlim; yl = ylim;
h = fill([t_ldown, t_ldown, t_relax, t_relax],[yl(1) yl(2) yl(2) yl(1)],'r');
% h = fill([t_ldown, t_ldown, t_relax, t_relax],[yl(1) yl(2) yl(2) yl(1)],0*[1 1 1]);
set(h,'facealpha',.1,'EdgeColor','none');

h = fill([t_relax, t_relax, xl(2), xl(2)],[yl(1) yl(2) yl(2) yl(1)],'g');
% h = fill([t_relax, t_relax, xl(2), xl(2)],[yl(1) yl(2) yl(2) yl(1)],0*[1 1 1]);
set(h,'facealpha',.1,'EdgeColor','none');
plot(t_relax*[1 1], Rt(t_relax,2:3), '*');

% % Plot hosp and cumulative mortality one on top of the other, to check
% % appropriate thresholds for initiating lockdown
% figure;
% subplot(2,1,1);
% plot(Ht(:,1)); xlim([0 250]); 
% subplot(2,1,2);
% y = sum(soln0(:,i.aux.mort),2);
% plot(y); xlim([0 250]);

