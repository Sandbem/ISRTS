%%INFO: 本函数为理论谱计算函数主程序。
%%----------------------------------------------------------------------%%
% Needs: ISR_init.m; 
%   ISR_updateFactors.m; ISR_updatePlasmas; ISR_updateDimensionless.m; 
%   ISR_gordeyev.m; ISR_spectrumKM.m; ISR_spectrumMac.m;
%%----------------------------------------------------------------------%%
% Inputs:
%   ion         - 离子成分 cell
%   ne          - 电子数密度 [m^-3]
%   Te          - 电子温度 [K]
%   Ti          - 离子温度 [K]
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
%%更新内容[Washy 2021/04/02]
% 1. 修改变量: 将Tr改为Te
% 2. 调换顺序: 将Ti放在Te后面
%%----------------------------------------------------------------------%%

function [spec,parameters]=ISR_main(ion,ne,Te,Ti,percent,frequency,fradar,theta,factors)
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
validateattributes(Te,{'double'},{'real','>',0})
validateattributes(Ti,{'double'},{'real','>',0})
validateattributes(percent,{'double'},{'row','>=',0,'<=',1})
validateattributes(frequency,{'double'},{'row','real'})
validateattributes(fradar,{'double'},{'real','>',0})

%% 初始化parameters
parameters = ISR_init(ion,ne,Te,Ti,percent,frequency,fradar,theta);

%% 额外因素
[parameters, mode] = ISR_updateFactors(parameters,factors);

%% 等离子体参数
parameters = ISR_updatePlasmas(parameters);

%% 无量纲参量
parameters = ISR_updateDimensionless(parameters);

%% Gordeyev积分
parameters = ISR_gordeyev(parameters,mode);

%% 理论谱
if parameters.factors.theomode == 1 % K&M
    [spec,parameters] = ISR_spectrumKM(parameters);
elseif parameters.factors.theomode == 2 % Mac
    [spec,parameters] = ISR_spectrumMac(parameters);
else
    error('参数 factors.theomode 设置错误!')
end

end