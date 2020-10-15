%%INFO: 本函数用来计算库仑碰撞频率。
%%----------------------------------------------------------------------%%
% Inputs:
%   densties    - 数密度 [ne, ni] [m^-3]
%   temperatures- 温度 [Te, Ti] [K]
%   ions        - 离子信息 cell
%   mode        - 计算模式, 默认为 'e_i'
% Outputs:
%   varargout   - 输出, 根据 mode 自动调整
%       计算电子: 输出只有一个参数 nue 包含 [nuee, nuei]
%       计算离子: 输出只有一个参数 nuii
%       全部计算: 输出包含两个参数 [nue, nuii]
%       nuii    - ion-ion库仑碰撞频率 [Hz]
%       nuei    - ele-ion库仑碰撞频率 [Hz]
%       nuee    - ele-ele库仑碰撞频率 [Hz]
%%----------------------------------------------------------------------%%
%参考文献
% [1] Rober Schunk and Andrew Nagy (2009). Ionospheres: Physics, Plasma
%   Physics, and Chemistry Second edition. [Book]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/01/02
%%----------------------------------------------------------------------%%
%更新内容[2019/08/20]
% 1. 修改ni输入单位为m^-3
% 2. 修改Ti的输入不随ion的元素个数改变
%%----------------------------------------------------------------------%%
%更新内容[2019/12/12]
% 1. 修改函数名为: CollisionCoulomb (原函数名为: Coulomb_C)
%%----------------------------------------------------------------------%%
%更新内容[2020/04/30]
% 1. 取消了csv文件读取, 添加Ion元胞组和Bst数组
%%----------------------------------------------------------------------%%
%更新内容[2020/08/03][重要]
%INFO: 对函数整体进行优化, 现在可以只计算电子[nuee, nuei]或离子[nuii]库仑
%   碰撞频率, 更为重要的是现在可以同时计算多个密度、温度对应的碰撞频率。
% 1. 修改输入为: densties, temperatures, ions, mode 
%   (原为: ion, Ti, Te, ni, ne)
% 2. 修改输出为: varargout 根据 mode 自动调整
%   (原为: nuii, nuei, nuee)
%%----------------------------------------------------------------------%%
%更新内容[2020/10/15]
%INFO: 修复了nuii计算错误的问题。第123行nust计算过程中，需要对先对ni(:,ia)
%   进行reshape，再计算。
%%----------------------------------------------------------------------%%
%输入示例[2020/08/03]
% ne = 1e10;
% ni = [0.5 0.2 0.3]*ne;
% Ti = 1000;
% Te = 1000;
% denstiy = [ne, ni; 0.1*ne, 0.1*ni; 0.2*ne, 0.2*ni];
% temperature = [Te, Ti; 0.95*Te, 0.95*Ti; 0.9*Te, 0.95*Ti];
% ions = {'O+','H+','NO+'};
% mode = 'e_i';
%%----------------------------------------------------------------------%%

function varargout = CollisionCoulomb(densties, temperatures, ions, mode)

% 输入参量个数判断
if nargin == 2
    ions = {};
    mode = 'e';
elseif nargin == 3
    mode = 'e_i';
elseif nargin < 2 || nargin > 4
    error('输入至少2个参量，至多4个参量！')
end

mode_e = {'e', 'ele', 'electron', 'E', 'Ele', 'Electron'};
mode_i = {'i', 'ion', 'ions', 'I', 'Ion', 'Ions'};
mode_e_i = {'e_i', 'ele_ion', 'ele_ions', 'E_I', 'Ele_Ion', 'Ele_Ions'};

% 计算电子碰撞频率
if ismember(mode, [mode_e, mode_e_i])
    % 提取电子参量
    ne = densties(:,1)*1e-6;
    Te = temperatures(:,1);
    
    % ele_ele碰撞频率
    nuee = 54.5/sqrt(2)*ne./Te.^1.5; %[1] eq(4.145)
    
    % ele_ion碰撞频率
    % nuei = 54.5*ni*Z'.^2/power(Te, 1.5); %[1] eq(4.144)
    nuei = 54.5*ne./Te.^1.5;

    % 输出
    varargout{1} = [nuee, nuei];
end

% 计算离子碰撞频率
if ismember(mode, [mode_i, mode_e_i])
    % 提取离子参量
    ni = densties(:, size(temperatures,2):end)*1e-6;
    Ti = temperatures(:,end);
    
    % 可计算离子种类
    defaultIons = {'H+','He+','C+','N+','O+','CO+','N2+','NO+','O2+','CO2+'};
    % 碰撞频率系数: 行 - s; 列 - t
    Bst = [0.900, 1.140, 1.220, 1.230, 1.230, 1.250, 1.250, 1.250, 1.250, 1.260;
        0.280, 0.450, 0.550, 0.560, 0.570, 0.590, 0.590, 0.600, 0.600, 0.610;
        0.102, 0.180, 0.260, 0.270, 0.280, 0.310, 0.310, 0.310, 0.310, 0.320;
        0.088, 0.160, 0.230, 0.240, 0.250, 0.280, 0.280, 0.280, 0.280, 0.300;
        0.077, 0.140, 0.210, 0.220, 0.220, 0.250, 0.250, 0.260, 0.260, 0.270;
        0.045, 0.085, 0.130, 0.140, 0.150, 0.170, 0.170, 0.170, 0.180, 0.190;
        0.045, 0.085, 0.130, 0.140, 0.150, 0.170, 0.170, 0.170, 0.180, 0.190;
        0.042, 0.080, 0.120, 0.130, 0.140, 0.160, 0.160, 0.160, 0.170, 0.180;
        0.039, 0.075, 0.120, 0.120, 0.130, 0.150, 0.150, 0.160, 0.160, 0.170;
        0.029, 0.055, 0.090, 0.090, 0.100, 0.120, 0.120, 0.120, 0.120, 0.140];
    
    % ion_ion碰撞频率
    [~, ia, ib] = intersect(ions, defaultIons, 'stable');
    lenTi = length(Ti);
    lenia = length(ia);
    lenib = length(ib);
    nuii = zeros(lenTi, length(ions));
%     for iTi = 1:lenTi
%         nust = Bst(ib, ib).*repmat(ni(iTi, ia), length(ia), 1)/Ti(iTi)^1.5; %[1] eq(4.143)
%         nuii(iTi, ia) = sum(nust, 2);
%     end
    nust = repmat(reshape(Bst(ib, ib),1,lenib,lenib), lenTi, 1, 1)...
        .*repmat(reshape(ni(:,ia), lenTi, 1, lenia), 1, lenia, 1)./repmat(Ti.^1.5, 1, lenia, lenia);
    nuii(:, ia) = sum(nust, 3);

    % 输出
    if ismember(mode, mode_i)
        varargout{1} = nuii;
    else
        varargout{2} = nuii;
    end
end

end