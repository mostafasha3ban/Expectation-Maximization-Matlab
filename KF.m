function [xkk, Pkk, x_pred, P_pred, Kgain] = KF(x0, P0, F, Q, H, R, zk) 
    if size(zk,2) > 1 % batch filter
        N = size(zk,2);
        Nx      = length(x0);
        xkk     = zeros(Nx, N);
        Pkk     = zeros(Nx, Nx, N);
        x_pred  = zeros(size(xkk)); 
        P_pred  = zeros(size(Pkk));
        [xkk(:,1), Pkk(:,:,1), x_pred(:,1), P_pred(:,:,1), Kgain(:,:,1)] =  ...
                                            KFiter(x0, P0, F, Q, H, R, zk(1));
        for i = 2:N
            [xkk(:,i), Pkk(:,:,i), x_pred(:,i), P_pred(:,:,i), Kgain(:,:,i)] =  ...
                      KFiter(xkk(:,i-1), Pkk(:,:,i-1), F, Q, H, R, zk(i));
        end 
    elseif size(zk,2)==1 % single run filter
        [xkk, Pkk, x_pred, P_pred, Kgain] = KFiter(x0, P0, F, Q, H, R, zk);
    end
end

function [xkk, Pkk, x_pred, P_pred, G] = KFiter(x0, P0, F, Q, H, R, zk)

    x_pred = F*x0;          %state prediction
    P_pred = F*P0*F' + Q ;   %state prediction covariance
    z_hat = H*x_pred;        %measurement prediction
    inov  = zk - z_hat;     %inovation
    S     = R + H*P_pred*H'; %innovation covariance
    G     = P_pred*H'*inv(S);%kalman gain
    xkk   = x_pred + G*inov; %state update
    Pkk   = P_pred - G*S*G'; %covariance update
    %Pkk   = (eye(size(Q,1)) - G*H)*P_pred*(eye(size(Q,1)) - G*H)' + G*R*G';
end