function [Q, R, F, H] = EMmax(x0n, P0n, xkn, Pkn, H, F, PL1, PL10, z, options )
    n = length(xkn);
    S11 = xkn(:,1)*xkn(:,1)' + Pkn(:,:,1);
    S10 = xkn(:,1)*x0n' + PL10;
    S00 = x0n*x0n' + P0n;
    for i = 2:n
        S11 =  S11 + xkn(:,i)*xkn(:,i)' + Pkn(:,:,i);
        S10 =  S10 + xkn(:,i)*xkn(:,i-1)' + PL1(:,:,i-1);
        S00 =  S00 + xkn(:,i-1)*xkn(:,i-1)' + Pkn(:,:,i-1);
    end 
    
    switch options.F
        case 'estimateF'
            Q =  (S11 - (S10*inv(S00)*S10'))/n;
            F = S10*inv(S00); %Maximizing for F
        case 'knownF'
            Q = (S11-(2*F*S10)+(F*F*S00))/n; % Assume F is known
    end
    
    switch options.H
        case 'estimateH'
            zk = z(:,1)*z(:,1);
            zkx = z(:,1)*xkn(:,1);
            xkp = xkn(:,1)*xkn(:,1) + Pkn(:,:,1);
            for i = 2:n
                 zk = zk + z(:,i)*z(:,i);
                 zkx = zkx + z(:,i)*xkn(:,i);
                 xkp = xkp + xkn(:,i)*xkn(:,i) + Pkn(:,:,i);
            end
             H = zkx/xkp; % maximizing for Hend
             R = (zk-(zkx^2/xkp))/n; % Using H Max
        case 'knownH'
            tmp = 0;
            for i = 1:n
                tmp = tmp  +  (z(:,i) - H*xkn(:,i))* (z(:,i) - H*xkn(:,i))' ...
                           + H*Pkn(:,:,i)*H';
            end
            R = tmp/n; 
    end
           
end