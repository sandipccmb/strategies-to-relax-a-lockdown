% Code used to create the dots needed for Figures 3, 4

clear all; 

Setup_model_India2;

tf         = 2000;
hosp_capac = 57194;
beta0      = r.beta;

R0         = 2.5; 

% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.(gps.geo{1}).ad) = seed; init(i.S.(gps.geo{1}).ad) = init(i.S.(gps.geo{1}).ad) - seed;


opts = odeset('Nonnegative',1:i.nx,'RelTol',1e-4,'AbsTol',1e-4); %The odeset function lets you adjust the integration parameters of the following ODE solvers.

% Deaths before lockdown; Effectiveness of lockdown; Proportion sympto;
% Realtive infectiousness of asympto; r.waning;  p.imm;
ranges  = [1 1500; 0.15 0.85; 2/3 0.9; 0.1 2/3; 0 1/(365/12*6); 0 0.1]';  % Deaths before lockdown; Effectiveness of lockdown; Proportion sympto; Realtive infectiousness of sympto
nsam    = 1000;
samples = repmat(ranges(1,:),nsam,1) + repmat(diff(ranges),nsam,1).*lhsdesign(nsam, size(ranges,2));

asym = 1;

mk = round(nsam/25);
tic
for isam = 1:size(samples,1)
    
    if mod(isam,mk)==0; fprintf('%0.5g ', isam/mk); end
        
    cbeta     = samples(isam,2); % reduction of transmission (0.15 0.85)
    p.sympto  = samples(isam,3); % proportion of symptomatic (2/3 0.9)
    p.c       = samples(isam,4); % Realtive infectiousness of asympto (0.1 2/3)
    r.beta    = beta0;
    
    r.waning  = samples(isam,5);
    p.imm     = samples(isam,6);
    
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
    t_ldown = find(tmort>samples(isam,1), 1, 'first');


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
    r2.q_asym = r2.gamma*asym;
    
    qdur0 = 1e-3; qdur1 = 100;
    
    ii = 1;

    while (ii < 100) && (qdur1-qdur0)>1e-1
        qdur = mean([qdur0,qdur1]);
        r2.q    = 1/qdur;
        M2      = make_model4(p2, r2, i, s, gps, prm);
        
        geq = @(t,in) goveqs_scaleup(t, in, M1, M2, i, s, p2, r2, agg, sel, t_relax+[0 15]);
        [t,tmp] = ode15s(geq, [t_relax:1:tf], soln1(t_relax,:), odeset('Nonnegative',1:i.nx));
        soln2 = soln1; soln2(t_relax+1:end,:) = tmp;
        
        % Check whether hospitalisation exceeds capacity
        hosp2   = sum(soln2(t_relax+15:end, s.H), 2);
        exceedH = (max(hosp2)>hosp_capac);
        % And whether there is a decreasing hosp trend
        vec = diff(hosp2);
        declining = (max(vec)<1e-3);
        
        if exceedH
            qdur1 = qdur;
        else
            qdur0 = qdur;
        end
        ii = ii+1;
    end
    
    dursto(isam) = qdur;
    
    
    % --- Pulling together data for outcomes
    ldown_sto = [t_ldown, t_relax, p.c];                                % Lockdown start and stop times, and effectiveness
    
    % Get the symptomatic incidence during lockdown
    sympto1 = -sum(diff(soln1([1,t_relax],s.S),[],1))*p.sympto;
    % Symptomatic incidence after lockdown
    sympto2 = -sum(diff(soln2([t_relax,end],s.S),[],1))*p.sympto;
    inc_sympto = [sympto1, sympto2];
    
    % symptomatic prevalence at time of lifting
    prev_sympto = sum(soln1(t_relax, s.prevalent));
    
    % Peak sympto prevalence after lockdown
    prev_peak = max(sum(soln2(t_relax:end,[s.IN1, s.IN2, s.IS1, s.IS2]),2),[],1);
    
    outmat(isam,:) = [ldown_sto, inc_sympto, prev_sympto, prev_peak];
    
end
toc

readme = 'v2 is new simulation done to include Hemanshu numbers';

if R0==2.5
    rc = '2p5';
elseif R0==3
    rc = '3';
end
fname = ['scanres_asym',num2str(asym),'_R0_',rc,'_v2'];
save(fname);

figure; plot(tR, dursto, '.'); ylim([0 20]);

fprintf('\n');
% save(['scan_declining_R',strrep(num2str(R0),'.','p')]);
% save(['scan_exceedH_R',strrep(num2str(R0),'.','p')]);



% --- Plot the overall result figures
% prev = squeeze(sum(allsol(:,s.prevalent,:),2));
% hosp = squeeze(sum(allsol(:,s.H,:),2));
% figure;
% subplot(1,2,1); plot(prev); xlim([0 1e3]);
% subplot(1,2,2); plot(hosp); xlim([0 1e3]);

allsol = cat(3,soln0,soln1,soln2);
hosp = squeeze(sum(allsol(:,s.H,:),2));
figure; hold on;
plot(hosp); 
line(t_ldown*[1 1], ylim, 'linestyle', '--');
line(t_relax*[1 1], ylim, 'linestyle', '--');
line(xlim, hosp_capac*[1 1], 'linestyle', '--');
xlim([0 500]);

