%%INFO: 本函数为理论谱计算函数主程序。
%%----------------------------------------------------------------------%%
% Needs: ISR_init.m; 
%   ISR_updateFactors.m; ISR_updatePlasmas; ISR_updateDimensionless.m; 
%   ISR_gordeyev.m; ISR_spectrumKM.m; ISR_spectrumMac.m;
%%----------------------------------------------------------------------%%
% Inputs:
%   ion         - 离子成分 cell
%   ne          - 电子数密度 [m^-3]
%   Ti          - 离子温度 [K]
%   Tr          - 电子离子温度比
%   percent     - 离子成分比例
%   frequency   - 频率范围 [Hz]
%   fradar      - 雷达频率 [Hz]
%   theta       - 散射角 [°]
%   factors     - 漂移、碰撞、磁场、库仑、模式 struct
% Output:
%   spec        - 理论谱 [s]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/09,30
%%----------------------------------------------------------------------%%

function [spec,parameters]=ISR_main(ion,ne,Ti,Tr,percent,frequency,fradar,theta,factors)
%% 输入检查
if nargin == 8
    factors = struct;
elseif nargin == 9
    validateattributes(factors,{'struct'},{'nonempty'})
elseif nargin < 8
    error('输入参量至少为8个！')
else
    error('输入参量至多为9个！')
end

validateattributes(ion,{'cell'},{'row'})
validateattributes(ne,{'double'},{'real','>',0})
validateattributes(Ti,{'double'},{'real','>',0})
validateattributes(Tr,{'double'},{'real','>',0})
validateattributes(percent,{'double'},{'row','>=',0,'<=',1})
validateattributes(frequency,{'double'},{'row','real'})
validateattributes(fradar,{'double'},{'real','>',0})

%% 初始化parameters
parameters = ISR_init(ion,ne,Ti,Tr,percent,frequency,fradar,theta);

%% 额外因素
[parameters, mode, theomode] = ISR_updateFactors(parameters,factors);

%% 等离子体参数
parameters = ISR_updatePlasmas(parameters);

%% 无量纲参量
parameters = ISR_updateDimensionless(parameters);

%% Gordeyev积分
[Je,Ji] = ISR_gordeyev(parameters, mode, theomode);

%% 积分及理论谱所需参量
k   = parameters.radar.k;
ni  = parameters.plasmas.ni;
vTe = parameters.plasmas.vTe;
vTi = parameters.plasmas.vTi;
he  = parameters.plasmas.he;

%% 理论谱
if theomode == 1 % K&M
    
    thetae = parameters.dimensionless.thetae;
    thetai = parameters.dimensionless.thetai;
    
    hi   = parameters.plasmas.hi;
    
    spec = ISR_spectrumKM(Je,Ji,thetae,thetai,k,ne,ni,vTe,vTi,he,hi);
 
elseif theomode == 2 % Mac
    
    qi     = parameters.plasmas.qi;
    qe     = parameters.plasmas.qe;
    omegae = parameters.plasmas.omegae;
    omegai = parameters.plasmas.omegai;
    
    psien  = parameters.dimensionless.psien;
    psiin  = parameters.dimensionless.psiin;
    
    eta = ni.*qi.^2/(ne*qe^2);
    mu  = eta*Tr;

    spec = ISR_spectrumMac(Je,Ji,omegae,omegai,k,ne,he,eta,mu,vTe,vTi,psien,psiin);
    
end

end