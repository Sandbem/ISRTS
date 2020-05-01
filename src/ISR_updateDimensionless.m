%%INFO: ISR_main.m的子函数, 用来更新parameters中dimensionless部分。
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - 计算过程产生的中间变量 struct
% Output:
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%

function parameters = ISR_updateDimensionless(parameters)

k = parameters.radar.k;

%% 电子

vTe  = parameters.plasmas.vTe;
rece = 1/k/vTe;

% 电子多普勒频率
parameters.dimensionless.thetae     = rece*parameters.plasmas.omegae;
% 电子-中性碰撞频率
parameters.dimensionless.psien      = rece*parameters.factors.nuen;
% 电子归一化回旋频率
parameters.dimensionless.phie       = rece*parameters.plasmas.Oe;
% 电子归一化库仑碰撞频率(平行)
parameters.dimensionless.psiec_par  = rece*parameters.factors.nuei;
% 电子归一化库仑碰撞频率(垂直)
parameters.dimensionless.psiec_perp = rece*(parameters.factors.nuei+parameters.factors.nuee);

%% 离子

vTi  = parameters.plasmas.vTi;
reci = 1/k./vTi;

% 离子多普勒频率
parameters.dimensionless.thetai     = reci'*parameters.plasmas.omegai;
% 离子-中性碰撞频率
parameters.dimensionless.psiin      = reci.*parameters.factors.nuin;
% 离子回旋频率
parameters.dimensionless.phii       = reci.*parameters.plasmas.Oi;
% 离子库仑碰撞频率
parameters.dimensionless.psiic      = reci.*parameters.factors.nuii;

end