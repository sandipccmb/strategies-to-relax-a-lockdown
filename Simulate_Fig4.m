% --- Use the code 'Make_dots' to create the inputs for this script

clear all;

ff=figure; fs = 14; ms = 15;

Rs = {'2p5', '3'};

cols = linspecer(2);

tis = {'R_0 = 2.5','R_0 = 3'};

for ir = 1:length(Rs)
    subplot(1,2,ir); hold on;

    fname = ['scanres_asym0_R0_',Rs{ir},'_v2'];
    load(fname);    
    scatter1 = scatter(tR,dursto,'MarkerFaceColor', cols(ir,:), 'MarkerEdgeColor','None');
    %scatter1.MarkerFaceAlpha = .5;

    fname = ['scanres_asym1_R0_',Rs{ir},'_v2'];
    load(fname);    
    scatter1 = scatter(tR,dursto,'MarkerFaceColor', cols(ir,:), 'MarkerEdgeColor','None');
    scatter1.MarkerFaceAlpha = .5;
    
    xlim([0 0.5]);
    ylim([0 20]);
    
    xlabel('Seroprevalence when lifting lockdown');
    if ir == 1
        ylabel({'Timeliness (duration of symptoms before quarantine)','to avoid second wave exceeding hospital capacity'});
    end
    set(gca,'fontsize',fs);
    title(tis{ir});
    
end

set(ff,'Position',[680   574   940   404]);
