function [xkn, Pkn, x0n, P0n, PL1, PL10]= KFsmooth(x0n, P0n, xkk, Pkk, xk_pred, Pk_pred, F, H, Kgain)
    Ns = size(xkk,2);
    Nx = length(x0n);
    xkn  = xkk;
    Pkn  = Pkk;
    C    = zeros(Nx, Nx, Ns-1);
    for k = Ns-1:-1:1
        C(:,:,k)    = Pkk(:,:,k)*F'*inv(Pk_pred(:,:,k+1));
        xkn(:,k)    = xkk(:,k)+C(:,:,k)*(xkn(:,k+1)-xk_pred(:,k+1));
        Pkn(:,:,k)  = Pkk(:,:,k)+C(:,:,k)*(Pkn(:,:,k+1)-Pk_pred(:,:,k+1))*C(:,:,k)'; 
    end
    C0 = P0n*F'/Pk_pred(1);
    x0n = x0n + C0*(xkn(:,1)-xk_pred(:,1));
    P0n = P0n + C0*(Pkn(:,:,1)-Pk_pred(:,:,1))*C0'; 
    
    % Lag-one Cov Smoother
    PL1     = zeros(Nx, Nx, Ns-1);
            W = Kgain(:,:,end);
    for k = Ns:-1:1
        if k==Ns
            PL1(:,:,k-1) = (eye(Nx) - W*H)*F*Pkk(:,:,end-1);
        elseif k>=2
            PL1(:,:,k-1) = Pkk(:,:,k-1)*C(:,:,k-1)' + ...
                        C(:,:,k)*(PL1(:,:,k) - F*Pkk(:,:,k-1))*C(:,:,k-1)';     
        elseif k==1
            PL10 = P0n*C0' + C(:,:,k)*(PL1(:,:,k) - F*P0n)*C0'; 
        end
    end
end