%% solve the problem 
%  Original Problem:
%      min_{U,V,E}  beta*||E||_2 + lambda* ||S||_1 + tau*||U||_{TV}
%                  s.t. Y = UV' + E + S, rank(U)=r,  V'*V=I
%  Auxiliary Problem:
%      min_{U,V,E} 1/2*||W.E||_2 + lambda* ||S||_1 + tau*(\sum_i ||G_i||_1)
%                  s.t. D_x(U) = G_1,D_y(U) = G_2,
%                       X = UV, rank(U)=r,  V'*V=I
%                          ===============================                      
%         D is difference operator
%  Jiangjun Peng, October 2022
%   If you have any question, please contact "andrew.pengjj@gmail.com" for
%   more information
%  ------------------------------------------------------------------------
function [DenoisedHSI,U,V] = RCTV(Nhsi, beta,lambda, tau, r)
%% Input Varaible
%   Nhsi   --- Polluted HSI data
%   beta   --- trade-off parameter of Gaussian noise
%   lambda --- trade-off parameter of sparse noise
%   rk     --- rank of needed-to-be recovered Tensor along mode-3
%% Output Variable
%   DenoisedHSI --- Denoised HSI
%   U           --- Representative Coefficient Image(RCI)
%   V           --- Orthogonal Basis
%% Initializing admm variables
tol     = 1e-5;
maxIter = 100;
rho     = 1.5;
max_mu  = 1e6;
[M,N,p] = size(Nhsi);
D       = zeros(M*N,p) ;
for i=1:p
    bandp = Nhsi(:,:,i);
    D(:,i)= bandp(:);
end

delta = 0.1;
if nargin < 2
    beta = 0.5;
end
if nargin < 3
    lambda = 1;
end
if nargin < 4
    tau = [0.01,0.01];
end
if nargin < 5
    r = 6;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf);
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
mu = 1/norm_two;
normD   = norm(D,'fro');
%% FFT setting
h               = M;
w               = N;
d               = r;
sizeU           = [h,w,d];
%% 
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
determ  =  Eny_x + Eny_y;
%% Initializing optimization variables

[u,s,v]= svd(D,'econ');
U              = u(:,1:r)*s(1:r,1:r);
V              = v(:,1:r);
X              = U*V';
S              = zeros(M*N,p);%Sparse
E              = zeros(M*N,p);%Gaussian

M1 = zeros(M*N*r,1);  % multiplier for Dx_U-G1
W1 = ones(M*N*r,1);
M2 = zeros(M*N*r,1);  % multiplier for Dy_U-G2
W2 = ones(M*N*r,1);
M3 = zeros([M*N,p]);  % multiplier for D-UV^T-E
W3 = ones(M*N,p);
% main loop
iter = 0;
eps =1e-4;
tic
while iter<maxIter
    iter          = iter + 1;  
    
    %% -Update G1,G2
    G1 = softthre_s(diff_x(U,sizeU)+M1/mu,tau(1)/mu,W1);
    G2 = softthre_s(diff_y(U,sizeU)+M2/mu,tau(2)/mu,W2);
    %% -Update U
    diffT_p  = diff_xT(G1-M1/mu,sizeU)+diff_yT(G2-M2/mu,sizeU);
    temp     = (D-E-S+M3/mu)*V;
    numer1   = reshape( diffT_p +temp(:), sizeU);
    x        = real( ifftn( fftn(numer1) ./ (determ + 1+eps) ) );
    U        = reshape(x,[M*N,r]);
    %% -Update V
    [u,~,v]     = svd((D-E-S+M3/mu)'*U,'econ');
    V           = u*v';
    %% -Update E
    E = mu*(D-U*V'-S+M3/mu)/(2*beta+mu);
    %% -Update S
    if lambda >100
        % 忽略掉高斯噪音影响
        S  = 0;
    else
        S = softthre_s(D-U*V'-E+M3/mu,lambda/mu,W3);
    end
    %% -Update Multiplier
    leq1 = diff_x(U,sizeU)-G1;
    leq2 = diff_y(U,sizeU)-G2;
    leq3 = D-U*V'-E-S;
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = norm(leq2,'fro')/normD;
    stopC3 = norm(leq3,'fro')/normD;
    if mod(iter,10)==0
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
                ',U_x rele = ' num2str(stopC1,'%2.3e') ',U_y rele = ' num2str(stopC2,'%2.3e')...
                ',X-UV = ' num2str(stopC3,'%2.3e')]);
    end
    if stopC1<tol && stopC2<tol && stopC3 <tol
        break;
    else
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        mu = min(max_mu,mu*rho); 
    end 
end
DenoisedHSI = reshape(U*V',[M,N,p]);
end
function out=softthre_s(a,tau,w)
out = sign(a).* max( abs(a) - tau*w, 0);
end
