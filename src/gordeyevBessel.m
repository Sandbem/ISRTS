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

dN = 1;
n = -N:dN:N;

eta = sin(alpha)^2/2/phi^2;
chi = besseli(n,eta,1);

n = repmat(n', 1, length(theta));
theta = repmat(theta, 2*N/dN+1, 1);

thetan = (theta - 1i*psin - n*phi)/cos(alpha);
G = sqrt(pi)/cos(alpha)*fadf(-thetan);

J = chi*G*dN;

end