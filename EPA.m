function  H = EPA(Y, M,  coefficients,sequence_Phi,F)%输入参数为Y为接收信号矩阵、X为发射信号矩阵、N为信号长度、B为一个参数、tau_ind_real真正时间延迟索引
% M: number of antennas
% J: number of grid points
% P: number of pilots
% L: number of delay taps
% G: number of frames
% etc: number of active grid points
%% grid parameters
[P, G] = size(Y);
D=P;
L=30;
J=100;
norm_Y = norm(Y, 'fro') / sqrt(P*G);
Y = Y / norm_Y;
max_doa = 1;
theta = ( -max_doa: 2*max_doa/J: max_doa - 2*max_doa/J ).';
theld_theta =  theta(2) - theta(1) ;
K = sqrt(M);
A = exp( -(0: M-1).'*1j*pi * (theta.') )/K ;
dA_theta=(-(0: M-1).'*1j*pi).*A;
S=coefficients;
B=A.'*S;

%% initialization
alpha = 1;
v_Apri = 1;
Z_Apri = zeros(L, G);
Gamma = ones(L, J);
U = zeros(L, J);
etc = 10;
iter = 1;
converged = false;
maxiter=300;
SSH=S*S';
XHX=F'* F/D;  % When sequence_Phi is an orthogonal pilot, XHX=F'* F/D.
%% beginning
while ~converged

    ZT_Sigma = zeros(G, G, L);
    ZT_post=zeros(L,G);
    % EPA
    %% A
    Sigma = zeros(L, L,G);
    Z_Apost = zeros(L, G);
    TraceSigma=zeros(G,1);
    resid=zeros(G,1);

    Sigma_g = pinv( 1/v_Apri * eye(L) + 1/alpha * XHX );
    for g=1:G
        X=diag(sequence_Phi(:,g))*F;
        X=X/sqrt(D);
        XHY= X'*Y(:,g);
        Sigma(:,:,g)=Sigma_g;
        Z_Apost(:,g) = Sigma_g * ( 1/v_Apri * Z_Apri(:,g) + 1/alpha * XHY);
        v_Apost=real(diag(Sigma_g));
        v_Apost2 = real( trace(Sigma_g) / L );
        TraceSigma(g)=real( trace( XHX * Sigma_g) );
        resid(g)=norm( Y(:,g) - X * Z_Apost(:,g), 'fro' )^2;
    end
    v_Aext= 1./ (1./v_Apost - 1/v_Apri);
    v_Aext2 = 1 / ( 1/v_Apost2 - 1/v_Apri);
    Z_Aext = v_Aext2 * ( 1/v_Apost2 * Z_Apost - 1/v_Apri * Z_Apri);

    alpha_old = alpha;
    rho_alpha = 0.98;
    alpha = sum(  resid + G*TraceSigma  ) / P / G;
    alpha = rho_alpha * alpha_old + (1 - rho_alpha) * alpha;

    %% B
    U_old=U;
    power = abs(  sum( U .* conj(U) )  );
    [~, active_ind] = sort(power);
    active_ind = active_ind( end: -1: end - etc + 1 );
    SS=0;
    for ll=1:L

        Dinv = diag(Gamma(ll,:));
        B_H = conj(B);
        B_T = B.';
        BTD=B_T * Dinv;
        MM = v_Aext(ll) * eye(G) + (BTD) * B_H;
        mid = MM \ BTD;
        U_Sigma_ll = Dinv - BTD' * mid;
        U_ll = Z_Aext(ll, :) * B' * (U_Sigma_ll.') / v_Aext(ll);
        SS = U_Sigma_ll(active_ind, active_ind)+SS;
        U(ll,:) = U_ll;
        Gamma(ll,:) = real( U_ll .* conj( U_ll ) ) + real( diag( U_Sigma_ll ) ).';
    end
    rho_v = 0.8;
    U = rho_v * U_old + (1-rho_v) * U;

    power = abs(U).^2;
    power_vet=power(:);
    [~,ind_power] = sort(power_vet);
    ind_shrange=ind_power(1:end-etc+1);
    if iter>100
        Gamma(ind_shrange)=  1e-8;
    end
    %% ZT
    for ll=1:L
        C=(B.') .* (ones(G,1)*Gamma(ll,:)) * conj(B);
        ZT_Sigma_ll = C * pinv(eye(G) + C/v_Aext(ll));
        ZT_post_ll=  Z_Aext(ll,:)*(ZT_Sigma_ll.')/ v_Aext(ll);
        ZT_Sigma(:,:,ll)= ZT_Sigma_ll;
        ZT_post(ll,:) = ZT_post_ll;
    end
    v_Tpos= real(sum(diag(sum(ZT_Sigma,3)))/L/G);
    v_Apri = 1. / ( 1/v_Tpos - 1./v_Aext2);
    Z_Apri = v_Apri * ( 1/v_Tpos * ZT_post - 1./v_Aext2 * Z_Aext);
    %% off-grid
    UU = U(:, active_ind).' * conj( U(:, active_ind) );
    UU_SS = UU + sum(SS, 3);
    CC_theta = dA_theta(:, active_ind).'*SSH * conj( dA_theta(:, active_ind) );
    P_theta = real( CC_theta .* UU_SS );
    term5 = 0;
    for jj = 1: L
        term5 = term5 + diag(   U(jj,active_ind)'  ) * dA_theta(:, active_ind)' *conj(S)* (( Z_Aext(jj,:).' - S.'*A(:, active_ind) * U(jj,active_ind).') );
    end
    v_theta = real(  term5 - diag( dA_theta(:, active_ind)'*conj(S)*S.' *  A(:, active_ind) * sum(SS, 3) )  );
    beta_theta = pinv(P_theta) * v_theta;


    Theld_theta=theld_theta/50*0.95^(iter);
    theta_small=find(abs(beta_theta)<Theld_theta);
    beta_theta(theta_small)=sign(beta_theta(theta_small))*Theld_theta;
    theta_max=find (abs(beta_theta)>theld_theta/5);
    beta_theta(theta_max)=sign(beta_theta(theta_max)) * theld_theta/5;

    theta(active_ind)=theta(active_ind)+beta_theta;

    A(:, active_ind)=exp( -(0: M-1).'*1j*pi*(theta(active_ind).')) / K;
    B(active_ind,:)  =A(:, active_ind).'*S;
    dA_theta(:, active_ind)=(-(0: M-1).'*1j*pi).*A(:, active_ind);
    %%
    if  iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
end
H = U * A.' * norm_Y/sqrt(D);


