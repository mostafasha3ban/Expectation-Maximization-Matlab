clear all; clc; close all

load data

meanQest = mean(Qest, 1);
meanRest = mean(Rest, 1);
meanFest = mean(Fest, 1);
meanHest = mean(Hest, 1);


h = figure(1);
subplot(221); hold on; box on;  plot(sig_w_range, Qtrue, 'linewidth', 2);
xlabel('\sigma_w'); ylabel('Q'); 
subplot(222); hold on; box on;  plot(sig_w_range, sqrt(Rtrue), 'linewidth', 2);
xlabel('\sigma_w'); ylabel('R'); 
subplot(223); hold on; box on;  plot(sig_w_range, Ftrue, 'linewidth', 2);
xlabel('\sigma_w'); ylabel('F'); 
subplot(224); hold on; box on;  plot(sig_w_range, Htrue, 'linewidth', 2);
xlabel('\sigma_w'); ylabel('H'); 
for iSim = 1:4
    subplot(221); plot(sig_w_range, squeeze(meanQest(:,iSim,:)), '-*'); 
    subplot(222); plot(sig_w_range, squeeze(meanRest(:,iSim,:)), '-*');
    subplot(223); plot(sig_w_range, squeeze(meanFest(:,iSim,:)), '-*');
    subplot(224); plot(sig_w_range, squeeze(meanHest(:,iSim,:)), '-*');
end

lgnd = {'True','Scen-1','Scen-2 (est. H)','Scen-3 (est. F)','Scen-4 (est.F&H)'};
subplot(221); legend(lgnd, 'location', 'best'); 
subplot(222); legend(lgnd, 'location', 'best'); 
subplot(223); legend(lgnd, 'location', 'best');
subplot(224); legend(lgnd, 'location', 'best'); 

