%%INFO: ISR_main.m的子函数, 用来更新parameters中plasmas部分。
%%----------------------------------------------------------------------%%
% Needs: constants.m; 
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - 计算过程产生的中间变量 struct
% Output:
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%更新内容[Washy 2020/05/01]
% 1. 去除输入: name, density, temperatur, constants
% 2. 简化计算过程, 去除if语句和for循环
%%----------------------------------------------------------------------%%

function parameters = ISR_updatePlasmas(parameters)

res  = constants({'kB', 'eps0'});
kB   = res.kB;
eps0 = res.eps0;

omega = parameters.omega;

k = parameters.radar.k;

B = parameters.factors.B;

%% 电子参数

ude = parameters.factors.ude;

me  = parameters.plasmas.me;
qe  = parameters.plasmas.qe;
ne  = parameters.plasmas.ne;
Te  = parameters.plasmas.Te;

% 电子热速度 [m/s]
parameters.plasmas.vTe      = sqrt(2*kB*Te/me);
% 电子等离子体频率 [Hz]
parameters.plasmas.wpe      = sqrt(ne*qe^2/eps0/me);
% 电子德拜半径 [m]
parameters.plasmas.he       = sqrt(eps0*kB*Te/ne/qe^2);
% 电子回旋频率 [Hz]
parameters.plasmas.Oe       = qe*B/me;
% 电子多普勒频率 [rad*Hz]
parameters.plasmas.omegae   = omega - k*ude;

%% 离子参数

udi = parameters.factors.udi;

mi  = parameters.plasmas.mi;
qi  = parameters.plasmas.qi;
ni  = parameters.plasmas.ni;
Ti  = parameters.plasmas.Ti;

% 离子热速度 [m/s]
parameters.plasmas.vTi      = sqrt(2*kB*Ti./mi);
% 离子等离子体频率 [Hz]
parameters.plasmas.wpi      = sqrt(ni.*qi.^2/eps0./mi);
% 离子德拜半径 [m]
parameters.plasmas.hi       = sqrt(eps0*kB.*Ti./ni./qi.^2);
% 离子回旋频率 [Hz]
parameters.plasmas.Oi       = qi*B./mi;
% 离子多普勒频率 [rad*Hz]
parameters.plasmas.omegai   = omega - k*udi;

end