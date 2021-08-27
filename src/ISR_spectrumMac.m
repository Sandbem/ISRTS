%%INFO: ISR_mian.m的子函数, 计算ISR理论谱。(Mac)
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - 计算过程产生的中间变量 struct
% Outputs:
%   sigma       - 微分散射截面 [s/m]
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
%参考文献
% [1] Swartz, W. E. & Farley, D. T. (1979). Theory of Incoherent-Scattering
%   of Radio-Waves by a Plasma .5. Use of the Nyquist Theorem in General 
%   Quasi-Equilibrium Situations. Journal of Geophysical Research-Space 
%   Physics, 1979, 84, 1930-932. doi:10.1029/JA084iA05p01930
% [2] Dougherty, J. P. & Farley, D. T. (1963). A Theory of Incoherent 
%   Scattering of Radio Waves by a Plasma .3. Scattering in a Partly 
%   Ionized Gas. Journal of Geophysical Research, 1963, 68, 5473-5486.
%   doi:10.1029/jz068i019p05473
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/10
%%----------------------------------------------------------------------%%
%%更新内容[Washy 2020/05/03]
% 1. 更改输入: parameters
% 2. 增加输出: parameters
%%----------------------------------------------------------------------%%

function [sigma,parameters] = ISR_spectrumMac(parameters)

k    = parameters.radar.k;
qe   = parameters.plasmas.qe;
ne   = parameters.plasmas.ne;
Te   = parameters.plasmas.Te;
vTe  = parameters.plasmas.vTe;
he   = parameters.plasmas.he;
omge = parameters.plasmas.omegae;
qi   = parameters.plasmas.qi;
ni   = parameters.plasmas.ni;
Ti   = parameters.plasmas.Ti;
vTi  = parameters.plasmas.vTi;
% hi   = parameters.plasmas.hi;
omgi = parameters.plasmas.omegai;
thee = parameters.dimensionless.thetae;
psien= parameters.dimensionless.psien;
thei = parameters.dimensionless.thetai;
psiin= parameters.dimensionless.psiin;
Je   = parameters.spectrum.Je;
Ji   = parameters.spectrum.Ji;

eta  = ni.*qi.^2/(ne*qe^2);
mu   = eta*Te/Ti;

% 经典电子半径re的平方re2
re2 = 7.940478951426530e-30; % e^2/(4*pi*eps0*me*c^2)

% 电子归一化导纳
ye = (1i+(thee-1i*psien).*Je)./(1-psien*Je); % [2] eq(5)
% 热运动引起的电子项散射谱
specTe = computeSpecT(ye,omge,psien,k,vTe,1);

% 热运动引起的离子项散射谱
row    = length(eta);  % 行数
column = length(omgi); % 列数
specTi = zeros(1, column);
muyi   = complex(zeros(1, column));

for irow = 1:row
    % 离子归一化导纳
    yi = (1i+(thei(irow,:)-1i*psiin(irow)).*Ji(irow,:))./(1-psiin(irow)*Ji(irow,:)); % [2] eq(5)
    
    specTi = specTi + computeSpecT(yi,omgi,psiin(irow),k,vTi(irow),eta(irow));
    muyi   = muyi + mu(irow)*yi;
end

den   = abs(ye+muyi+1i*k^2*he^2).^2;
parameters.spectrum.spece = abs(muyi+1i*k^2*he^2).^2./den.*specTe;
parameters.spectrum.speci = abs(ye).^2./den.*specTi;
% sigma = ne*re2*(spece + speci); % [1] eq(20)
sigma = 2*pi*(parameters.spectrum.spece + parameters.spectrum.speci); % [1] eq(20)

end

%%----------------------------------------------------------------------%%
%%INFO: ISR_spectrumMac.m的内含子函数，用来计算热运动引起的散射谱，并消除
%   中心点位置为NaN的问题。
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   y           - 无量纲导纳
%   omg         - 多普勒频移 [rad*Hz]
%   psin        - 与中性成分碰撞频率 [Hz]
%   k           - 散射差矢 [m^-1]
%   vT          - 热速度 [m/s]
%   eta         - 与离子带电量相关系数
% Outputs:
%   specT       - 热运动引起的项散射谱 [s]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/10
%%----------------------------------------------------------------------%%

function specT = computeSpecT(y,omg,psin,k,vT,eta)

specT = eta/pi*real(y)./omg; % [1] eq(20)

% 重新计算NaN位置的数值
if ~all(omg)
    imagG = sqrt(pi)*fadf(1i*psin);
    specT(omg==0) = eta/(pi*k*vT)*imagG/(1-psin*imagG);
end

end