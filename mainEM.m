clear all; clc; close all
N = 250;
dbstop('error')

modelID = 1;
sig_v   = 1; % kept always a constant
sig_w   = 1;
F_range = [0.1:0.1:0.9];
H_range = [1];
nMCruns = 1;

Fest = zeros(nMCruns, 4, length(F_range), length(H_range));
Hest = zeros(nMCruns, 4, length(F_range), length(H_range));
Qest = zeros(nMCruns, 4, length(F_range), length(H_range));
Rest = zeros(nMCruns, 4, length(F_range), length(H_range));
LLHest = zeros(nMCruns, 4, length(F_range), length(H_range));

Ftrue = zeros(1, length(F_range), length(H_range));
Htrue = zeros(1, length(F_range), length(H_range));
Qtrue = zeros(1, length(F_range), length(H_range));
Rtrue = zeros(1, length(F_range), length(H_range));
LLHtrue = zeros(nMCruns, length(F_range), length(H_range));

tic
for iH = 1:length(H_range)
for iF = 1:length(F_range)  % simulation scenario
    F = F_range(iF); 
    Q = sig_v^2;
    H   = H_range(iH);
    R = sig_w^2;
    x = zeros(1, N);x(1,1) = 0;  
    z = zeros(1,N); z(:,1)  = x(1,1);
    
    Ftrue(iF,iH) = F;
    Htrue(iF,iH) = H;
    Qtrue(iF,iH) = Q;
    Rtrue(iF,iH) = R;
    
    for iMC = 1:nMCruns
                
            for k = 2:N
                x(:,k) = F*x(:,k-1) + sqrt(Q)*randn(size(Q,1),1);
                z(:,k) = H*x(:,k) + sqrt(R)*randn(size(R,1),1);
            end
            
            sum1 = x(:,1);
            sum2 = (z(:,1) - H*x(:,1)).^2;
            
                for i=2:N
                 sum1 = sum1 + (x(:,i) - F*x(:,i-1)).^2;
                 sum2 = sum2 + (z(:,i) - H*x(:,i)).^2;
                end

        c = - (2 + 4*N)*log(1/sqrt(2*pi));
        LLHtrue(iMC,iF,iH) = -(N*log(Q) + (Q^-1*sum1) + N*log(R) + (R^-1*sum2) + c)/2;
    
        for iSim = 1:4
            iSim;
        switch iSim
        case 1
            options.F = 'knownF'; options.H = 'knownH';
            Fem = F; Hem = H;
        case 2
            options.F = 'knownF'; options.H = 'estimateH';
            Fem = F; Hem =  5*rand*H;
        case 3
            options.F = 'estimateF'; options.H = 'knownH';
            Fem = rand*F; Hem = H;
        case 4
            options.F = 'estimateF'; options.H = 'estimateH';
            Fem = rand*F; Hem =  5*rand*H;
        end

        x0 = 0;
        P0 = 0.1;
        xkk = zeros(1,N); x(1) = x0;
        Pkk = zeros(1,N); Pkk(1) = P0;
                    
        Qem = 5*rand*Q; % initialize away from the trut (other means is better)
        Rem = 5*rand*R;  % initialize away from the trut (other means is better)
        Q1  = Qem;
        R1  = Rem;
            
            % the actual EM
            x0n = x0;
            P0n = P0;
            i=0;
            while i<50
                [xkkEM, PkkEM, x_predEM, P_predEM, kgainEM] = KF(x0n, P0n,Fem, Qem, Hem, Rem, z);
                [xkn, Pkn, x0n, P0n, PL1, PL10]= ...
                    KFsmooth(x0, P0, xkkEM, PkkEM, x_predEM, P_predEM, Fem, Hem, kgainEM);
                [Qem, Rem, Fem, Hem] =  EMmax(x0n, P0n, xkn, Pkn, Hem, Fem, PL1, PL10, z, options);
                i=i+1;
            end

            sum1 = x(:,1);
            sum2 = (z(:,1) - H*x(:,1)).^2;
                for i=2:N
                   sum1 = sum1 + (x(:,i) - Fem*x(:,i-1)).^2;
                   sum2 = sum2 + (z(:,i) - Hem*x(:,i)).^2;
                end


          llh1 = -( N*log(Qem) + (Qem^-1*sum1) + N*log(Rem) + (Rem^-1*sum2) + c)/2;
         
            Qest(iMC, iSim, iF, iH) = Qem;
            Rest(iMC, iSim, iF, iH) = Rem;
            Fest(iMC, iSim, iF, iH) = Fem;
            Hest(iMC, iSim, iF, iH) = Hem;
            LLHest(iMC,iSim,iF, iH) = llh1;

        end
    end
end
end

meanLLH = mean(LLHest,1);

save data1000badINI

toc
h= figure; hold on; grid on; box on;
plot(F_range, mean(LLHtrue,1),'--','linewidth', 1)
plot(F_range,squeeze(meanLLH(1,1,:)),'-*','linewidth', 1)
plot(F_range,squeeze(meanLLH(1,2,:)),'-s','linewidth', 1)
plot(F_range,squeeze(meanLLH(1,3,:)),'-d','linewidth', 1)
plot(F_range,squeeze(meanLLH(1,4,:)),'-+','linewidth', 1)
legend({'True LLh','Case 1', 'Case 2', 'Case 3', 'Case 4'}, 'location', 'best')
xlabel('F')
ylabel('Likelihood')
set(gca, 'fontsize', 16)

fig_name = ['./Fig/Nruns=' num2str(nMCruns) 'badINI'];
print(h, fig_name, '-depsc')
saveas(h, fig_name, 'fig')

Qerror = Qest - 1;
Qerror = Qerror.^2;
meanQerror = squeeze(mean(Qerror,1));
rmseQ = sqrt(meanQerror);

Rerror = Rest - 1;
Rerror = Rerror.^2;
meanRerror = squeeze(mean(Rerror,1));
rmseR = sqrt(meanRerror);

Herror = Hest - 1;
Herror = Herror.^2;
meanHerror = squeeze(mean(Herror,1));
rmseH = sqrt(meanHerror);


Ferror = zeros(nMCruns, 4, length(F_range));
for i=1:length(F_range)
    Ferror(:,:,i) = Fest(:,:,i) - F_range(i);
end
Ferror = Ferror.^2;
meanFerror = squeeze(mean(Ferror,1));
rmseF = sqrt(meanFerror);

h= figure; hold on
box on; grid on
plot(F_range,rmseQ(1,:),'-*')
plot(F_range,rmseR(1,:),'-+')
legend({'Q','R'}, 'location', 'best')
xlabel('F')
ylabel('RMSE')
set(gca, 'fontsize', 14)
title('Case 1')
fig_name = ['./Fig/RMSECase1' num2str(nMCruns) 'badINI'];
print(h, fig_name, '-depsc')
saveas(h, fig_name, 'fig')

h= figure; hold on
plot(F_range,rmseQ(2,:),'-*')
plot(F_range,rmseR(2,:),'-+')
plot(F_range,rmseH(2,:),'-d')
box on; grid on
xlabel('F')
ylabel('RMSE')
set(gca, 'fontsize', 14)
title('Case 2')
legend({'Q','R','H'}, 'location', 'best')
fig_name = ['./Fig/RMSECase2' num2str(nMCruns) 'badINI'];
print(h, fig_name, '-depsc')
saveas(h, fig_name, 'fig')

h= figure; hold on
plot(F_range,rmseQ(3,:),'-*')
plot(F_range,rmseR(3,:),'-+')
plot(F_range,rmseF(3,:),'-s')
box on; grid on
xlabel('F')
ylabel('RMSE')
set(gca, 'fontsize', 14)
title('Case 3')
legend({'Q','R','F'}, 'location', 'best')
fig_name = ['./Fig/RMSECase3' num2str(nMCruns) 'badINI'];
print(h, fig_name, '-depsc')
saveas(h, fig_name, 'fig')

h= figure; hold on
plot(F_range,rmseQ(4,:),'-*')
plot(F_range,rmseR(4,:),'-+')
plot(F_range,rmseF(4,:),'-d')
plot(F_range,rmseH(4,:),'-s')
box on; grid on
xlabel('F')
ylabel('RMSE')
set(gca, 'fontsize', 14)
title('Case 4')
legend({'Q','R','H','F'}, 'location', 'best')
fig_name = ['./Fig/RMSECase4' num2str(nMCruns) 'badINI'];
print(h, fig_name, '-depsc')
saveas(h, fig_name, 'fig')