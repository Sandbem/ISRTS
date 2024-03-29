%%INFO: ISR_main.m的子函数, 用来更新parameters中dimensionless部分。
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - 计算过程产生的中间变量 struct
% Output:
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
%参考文献
% [1] Dougherty, J. P. & Farley, D. T. (1963). A Theory of Incoherent 
%   Scattering of Radio Waves by a Plasma .3. Scattering in a Partly 
%   Ionized Gas. Journal of Geophysical Research, 1963, 68, 5473-5486.
%   doi:10.1029/jz068i019p05473
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
parameters.dimensionless.thetae     = rece*parameters.plasmas.omegae; %[1] eq(7)
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
parameters.dimensionless.thetai     = reci'*parameters.plasmas.omegai; %[1] eq(7)
% 离子-中性碰撞频率
parameters.dimensionless.psiin      = reci.*parameters.factors.nuin;
% 离子回旋频率
parameters.dimensionless.phii       = reci.*parameters.plasmas.Oi;
% 离子库仑碰撞频率
parameters.dimensionless.psiic      = reci.*parameters.factors.nuii;

end