%%INFO: 使用Bessel函数求解Gordeyev积分。
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   theta       - 自变量
%   psin        - 与中性成分碰撞频率 [Hz]
%   alpha       - 散射差矢与磁场夹角 [rad]
%   phi         - 回旋频率 [Hz]
%   N           - 求和最大个数
% Output:
%   J           - 积分结果
%%----------------------------------------------------------------------%%
% author: Washy [IGG]
% date: 2020/05/07
%%----------------------------------------------------------------------%%

function J = gordeyevBessel(theta,psin,alpha,phi,N)

if nargin < 5
    N = 8;
end

%%
Nm = 4096;
epsJN = 1e-6;

eta = 0.5*sin(alpha)^2/phi^2;

while 1
    JN = abs(besseli(N,eta,1));
%     JN = abs(besseli(N,eta,1)*sqrt(pi)/cos(alpha)*fadf(N*phi));
    if JN < epsJN || N >= Nm
        break
    end
    N = 2*N;
end

%%
n = -N:N;
chi = besseli(n,eta,1);

n = repmat(n', 1, length(theta));
theta = repmat(theta, 2*N+1, 1);

thetan = (theta - 1i*psin - n*phi)/cos(alpha);
G = sqrt(pi)/cos(alpha)*fadf(-thetan);

J = chi*G;

end