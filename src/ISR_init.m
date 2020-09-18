%%INFO: ISR_main.m的子函数, 用来初始化parameters以及其中radar部分的计算。
%%----------------------------------------------------------------------%%
% Needs: constants.m; analysisIon.m; 
%%----------------------------------------------------------------------%%
% Inputs:
%   frequency   - 频率范围 [Hz]
%   fradar      - 雷达频率 [Hz]
%   theta       - 散射角 [°]
% Output:
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%%更新内容[Washy 2020/01/09]
% 1. 增加变量: parameters.dimensionless.He, parameters.dimensionless.Hi
%%----------------------------------------------------------------------%%
%%更新内容[Washy 2020/05/01]
% 1. 增加输入: ion, ne, Ti, Tr, percent
%%----------------------------------------------------------------------%%

function parameters = ISR_init(ion,ne,Ti,Tr,percent,frequency,fradar,theta_s)

res = constants({'c','me','e'});

parameters = struct;

% parameters.frequency = frequency;           % 多普勒频率 [Hz]
parameters.omega = 2*pi*frequency;          % 多普勒频率 [rad*Hz]

%% 雷达参数
parameters.radar = struct;

parameters.radar.frequency = fradar;        % 雷达频率 [Hz]
parameters.radar.theta     = theta_s;       % 散射角 [°]
parameters.radar.k         = 4*pi*fradar/res.c*sin(theta_s*pi/360); % 散射差矢 [1/m]

%% 额外因素
parameters.factors = struct;

parameters.factors.ude    = [];             % 电子漂移速度(视线方向)
parameters.factors.udi    = [];             % 离子漂移速度(视线方向)

parameters.factors.nuen   = [];             % 电子-中性碰撞频率 [Hz]
parameters.factors.nuin   = [];             % 离子-中性碰撞频率 [Hz]

parameters.factors.B      = [];             % 磁感应强度 [T]
parameters.factors.alpha  = [];             % 与磁场夹角 [°]

parameters.factors.nuii   = [];             % 离子-离子库仑碰撞频率 [Hz]
parameters.factors.nuei   = [];             % 电子-离子库仑碰撞频率 [Hz]
parameters.factors.nuee   = [];             % 电子-电子库仑碰撞频率 [Hz]

% parameters.factors.mode   = [];             % 控制开关: 0 关; 1 开
parameters.factors.gordmode = [];           % gordeyev积分: 1: Sommerfeld; 2: Bessel
parameters.factors.theomode = [];           % 理论谱计算公式: 1: Kudeki & Milla; 2: Farley et al.

%% 等离子体参数
parameters.plasmas = struct;

% 电子参数
parameters.plasmas.namee  = {'E-'};         % 电子种类名称 cell
parameters.plasmas.me     = res.me;         % 电子质量 [kg]
parameters.plasmas.qe     = res.e;          % 电子电荷量 [C]

parameters.plasmas.ne     = ne;             % 电子数密度 [m^-3]
parameters.plasmas.Te     = Ti*Tr;          % 电子温度 [K]

parameters.plasmas.vTe    = [];             % 电子热速度 [m/s]
parameters.plasmas.wpe    = [];             % 电子等离子体频率 [Hz]
parameters.plasmas.he     = [];             % 电子德拜半径 [m]
parameters.plasmas.Oe     = [];             % 电子回旋频率 [Hz]

parameters.plasmas.omegae = [];             % 电子多普勒频率 [rad*Hz]

% 离子参数
[mi, qi] = analysisIon(ion);
parameters.plasmas.namei  = ion;            % 离子种类名称 cell
parameters.plasmas.mi     = mi';            % 离子质量 [kg]
parameters.plasmas.qi     = qi';            % 离子电荷量 [C]

parameters.plasmas.ni     = ne*percent;     % 离子数密度 [m^-3]
parameters.plasmas.Ti     = Ti;             % 离子温度 [K]

parameters.plasmas.vTi    = [];             % 离子热速度 [m/s]
parameters.plasmas.wpi    = [];             % 离子等离子体频率 [Hz]
parameters.plasmas.hi     = [];             % 离子德拜半径 [m]
parameters.plasmas.Oi     = [];             % 离子回旋频率 [Hz]

parameters.plasmas.omegai = [];             % 离子多普勒频率 [rad*Hz]

%% 无量纲参量
parameters.dimensionless = struct;

parameters.dimensionless.thetae     = [];   % 电子无量纲自变量
parameters.dimensionless.thetai     = [];   % 离子无量纲自变量

parameters.dimensionless.psien      = [];   % 电子-中性无量纲碰撞频率
parameters.dimensionless.psiin      = [];   % 离子-中性无量纲碰撞频率

parameters.dimensionless.phie       = [];   % 电子无量纲回旋频率
parameters.dimensionless.phii       = [];   % 离子无量纲回旋频率

parameters.dimensionless.psiec_par  = [];   % 电子无量纲库仑碰撞频率(平行)
parameters.dimensionless.psiec_perp = [];   % 电子无量纲库仑碰撞频率(垂直)
parameters.dimensionless.psiic      = [];   % 离子无量纲库仑碰撞频率

% parameters.dimensionless.He         = [];   % 电子无量纲
% parameters.dimensionless.Hi         = [];   % 离子无量纲

%% 理论谱参量
parameters.spectrum = struct;

parameters.spectrum.Je    = [];             % 电子Gordeyev积分
parameters.spectrum.Ji    = [];             % 离子Gordeyev积分

parameters.spectrum.spece = [];             % 电子散射谱
parameters.spectrum.speci = [];             % 离子散射谱

end