%%INFO: ISR_main.m的子函数, 用来更新parameters中factors部分。
%%----------------------------------------------------------------------%%
% Needs: CollisionCoulomb.m
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - 计算过程产生的中间变量 struct
%   factors     - 漂移、碰撞、磁场、模式 struct
% Outputs:
%   parameters  - 计算过程产生的中间变量 struct
%   mode        - 模式控制(关: 0; 开: 1)
%   theomode    - 理论谱计算模式(Mac: 1; KM: 2)
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%更新内容[Washy 2019/09/24]
%%INFO: 将碰撞频率nu从单个数值输入改为行向量输入，其中第一个数值表示电子
%   碰撞频率，第二个至最后一个表示离子碰撞频率。
% 1. 修改了nu的数据长度；
% 2. 修改了if语句mode(2)中的赋值方式。
%%----------------------------------------------------------------------%%
%更新内容[Washy 2019/10/08]
%%INFO: 将ISRKM_analysis_factors.m合并进了该函数。并对函数的整体结构进行
%   了重新规划。
% 1. 去除输入变量: ud, nu, B, mode
% 2. 增加输入变量: factors
% 3. 增加输出变量: mode, gordmode
% 4. 修改变量factors.ud: 区分电子、离子输入，修改后的格式为[ude, udi]，
%   udi只有一个值，即默认所有成分具有相同的值。
%%----------------------------------------------------------------------%%
%更新内容[Washy 2019/10/10]
% 1. 增加输出变量: theomode
%%----------------------------------------------------------------------%%
%更新内容[Washy 2020/05/01]
% 1. 去除输入: ion, ne, ni, Te, Ti
% 2. 去除输出: gordmode
%%----------------------------------------------------------------------%%

function [parameters, mode, theomode] = ISR_updateFactors(parameters,factors)

leni = length(parameters.plasmas.namei);

% 缺省值
B = [3.5e-5, 0];

% 模式控制(关: 0; 开: 1)
if isfield(factors,'mode')
    validateattributes(factors.mode,{'double'},{'row','binary'})
    mode = factors.mode;
    
    % 漂移速度
    if isfield(factors,'ud') && mode(1)
        parameters.factors.ude    = factors.ud(1);
        parameters.factors.udi    = factors.ud(2);
    else
        parameters.factors.ude    = 0;
        parameters.factors.udi    = 0;
    end
    
    % 碰撞频率
    if isfield(factors,'nu') && mode(2)
        parameters.factors.nuen = factors.nu(1);
        parameters.factors.nuin = factors.nu(2:end);
    else
        parameters.factors.nuen = 0;
        parameters.factors.nuin = zeros(1,leni);
    end
    
    % 磁场
    if isfield(factors,'B') && mode(3)
        parameters.factors.B     = factors.B(1);
        parameters.factors.alpha = factors.B(2);
    else
        parameters.factors.B     = B(1);
        parameters.factors.alpha = B(2);
    end
    
    % 库仑碰撞频率
    if mode(4)
        [parameters.factors.nuii, parameters.factors.nuei, ...
            parameters.factors.nuee] = ...
            CollisionCoulomb(parameters.plasmas.namei, ...
            parameters.plasmas.Ti, parameters.plasmas.Te, ...
            parameters.plasmas.ni, parameters.plasmas.ne);
    else
        parameters.factors.nuii = zeros(1,leni);
        parameters.factors.nuei = 0;
        parameters.factors.nuee = 0;
    end
    
else
    mode = [0,0,0,0];
    
    parameters.factors.ude   = 0;
    parameters.factors.udi   = 0;
    
    parameters.factors.nuen  = 0;
    parameters.factors.nuin  = zeros(1,leni);
    
    parameters.factors.B     = B(1);
    parameters.factors.alpha = B(2);
    
    parameters.factors.nuii  = zeros(1,leni);
    parameters.factors.nuei  = 0;
    parameters.factors.nuee  = 0;
    
end

parameters.factors.mode = mode;

% 理论谱计算公式选择
if isfield(factors,'theomode')
    validateattributes(factors.theomode,{'double'},{'scalar','integer','>=',1,'<=',2})
    theomode = factors.theomode;
else
    theomode = 1;
end

end