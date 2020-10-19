% --- Use the code 'Make_dots' to create the inputs for this script

clear all;

load scanres_asym0_R0_2p5_v2

% Find a point having cumul incidence closest to xref%
xref = 0.2;
dif  = (tR-xref).^2; ind = find(dif==min(dif));

vec = samples(ind,:);

qdur_list = [4,11];

for ii = 1:length(qdur_list)
    
    if mod(isam,mk)==0; fprintf('%0.5g ', isam/mk); end
    
    cbeta     = vec(2);
    p.sympto  = vec(3);
    p.c       = vec(4);
    r.beta    = beta0;
    
    r.waning  = vec(5);
    p.imm     = vec(6);
    
    
    p.q = 1; r.q = 0;
    
    M = make_model4(p, r, i, s, gps, prm);
    tmp = prm.N';
    I_inds = intersect([s.E, s.IA, s.IP, s.IN1, s.IN2, s.IS1, s.IS2],s.del);
    R = find_R0(M, tmp(:), I_inds, s.E);
    
    r.beta = r.beta*R0/R;
    
    % --- Baseline (no intervention)
    M0 = make_model4(p, r, i, s, gps, prm);
    geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
    [t,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
    
    % Find date of first death
    tmort   = sum(soln0(:,i.aux.mort),2);
    t_ldown = find(tmort>vec(1), 1, 'first');
    
    
    % --- With lockdown
    r1 = r; p1 = p;
    r1.beta = r.beta*(1-cbeta);
    M1 = make_model4(p1, r1, i, s, gps, prm);
    
    geq = @(t,in) goveqs_scaleup(t, in, M0, M1, i, s, p1, r, agg, sel, t_ldown+[0 10]);
    [t, soln1] = ode15s(geq, [0:1:tf], init, opts);
    
    % Find when the hospitalisations have come down to 10% of capacity again
    hosp = sum(soln1(:,s.H),2);
    if max(hosp)>0.1*hosp_capac
        inds = find(hosp>0.1*hosp_capac);
        tmp = inds(end);
    else
        tmp = find(hosp==max(hosp))+14;
    end
    t_relax = max(tmp, t_ldown + 21);
    
    % What proportion are immune by the time of relaxing the lockdown
    tR(isam) = sum(soln1(t_relax,s.R))/sum(prm.N(1,:));
    
    % --- Relaxing the lockdown
    r2 = r1; p2 = p1;
    r2.beta = r.beta;
    r2.q_asym = r2.gamma;
    
    r2.q = 1/qdur_list(ii);
    
    M2 = make_model4(p2, r2, i, s, gps, prm);
    geq = @(t,in) goveqs_scaleup(t, in, M1, M2, i, s, p2, r2, agg, sel, t_relax+[0 15]);
    [t,tmp] = ode15s(geq, [t_relax:1:tf], soln1(t_relax,:), odeset('Nonnegative',1:i.nx));
    soln2 = soln1; soln2(t_relax+1:end,:) = tmp;
    
    allsto(:,:,:,ii) = cat(3,soln0,soln1,soln2);
    
end

% Multiplier to put cases per 1,000 population
mult = 1e3/sum(prm.N);

y = squeeze(sum(allsto(:,s.H,:,:),2));

ff=figure; lw = 1.5; fs = 14;
%set(0,'DefaultLegendAutoUpdate','off');

subplot(2,2,[1,3]); hold on;
plot(tR, dursto, '.', 'markersize', 15);
xlabel('Proportion immune when lifting lockdown');
ylabel({'Timeliness (duration of symptoms before quarantine)','to avoid second wave exceeding hospital capacity'});
xlim([0 0.4]);
ylim([0 20]);
set(gca,'fontsize',fs);
plot(xref*[1 1], qdur_list, '*','Color',[0.2 0.2 0.2],'Markersize',10,'linewidth',2);
text(xref-0.03,qdur_list(2)-0.5,'P1','fontsize',14);
text(xref-0.03,qdur_list(1)-0.5,'P2','fontsize',14);

subplot(2,2,2); hold on;
plot(y(:,[1,3],2)*mult, 'linewidth', lw);
xlim([50 600]);
line(xlim, hosp_capac*[1 1]*mult, 'linestyle', '--');
title('P1 scenario','fontsize',3);
set(gca,'fontsize',fs);
legend('No lockdown','Lockdown','Hospital bed capacity');
yl = ylim; xl = xlim;
h1 = fill([t_ldown, t_ldown, t_relax, t_relax], [yl, fliplr(yl)], 'r');
h2 = fill([t_relax, t_relax, xl(2), xl(2)], [yl, fliplr(yl)], 'g');
set(h1,'facealpha',.1,'EdgeColor','none');
set(h2,'facealpha',.1,'EdgeColor','none');
ylabel({'Numbers needing hospitalisation','(per thousand popn)'},'fontsize',13);


subplot(2,2,4); hold on;
plot(y(:,[1,3],1)*mult, 'linewidth', lw);
xlim([50 600]);
line(xlim, hosp_capac*[1 1]*mult, 'linestyle', '--');
yl = ylim; xl = xlim;
h1 = fill([t_ldown, t_ldown, t_relax, t_relax], [yl, fliplr(yl)], 'r');
h2 = fill([t_relax, t_relax, xl(2), xl(2)], [yl, fliplr(yl)], 'g');
set(h1,'facealpha',.1,'EdgeColor','none');
set(h2,'facealpha',.1,'EdgeColor','none');
title('P2 scenario','fontsize',10);
set(gca,'fontsize',fs);
ylabel({'Numbers needing hospitalisation','(per thousand popn)'},'fontsize',13);

set(ff,'Position',[573   454   961   483])

% load scanres_R0_3
% % load scanres_asym_R0_3
% plot(tR, dursto, 'r.');
%
% ylim([0 20]);