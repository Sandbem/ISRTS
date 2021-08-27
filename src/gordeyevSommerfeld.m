%%INFO: 使用Sommerfeld积分求解Gordeyev积分。
%%----------------------------------------------------------------------%%
% Needs: SingleParticleACF.m; SommerfeldIntegral.m; 
%%----------------------------------------------------------------------%%
% Inputs:
%   theta       - 自变量
%   psin        - 与中性成分碰撞频率 [Hz]
%   alpha       - 散射差矢与磁场夹角 [rad]
%   phi         - 回旋频率 [Hz]
%   psic_para   - 库仑碰撞频率 (平行于磁场) [Hz]
%   psic_perp   - 库仑碰撞频率 (垂直于磁场) [Hz]
% Output:
%   J           - 积分结果
%%----------------------------------------------------------------------%%
% author: Washy [IGG]
% date: 2020/05/08
%%----------------------------------------------------------------------%%

function J = gordeyevSommerfeld(theta,psin,alpha,phi,psic_para,psic_perp)

% 被积函数
func = @SingleParticleACF;
% 积分下限
a = 0;
% 积分上限
b = 1;
% 分段个数
N = 32;
% 约束条件
restriction = [1e-9,512,1e-6,32768];

funPar = [psin, phi, alpha, psic_para, psic_perp];

[J,~,~] = SommerfeldIntegral(func,a,b,N,theta,funPar,restriction);

end