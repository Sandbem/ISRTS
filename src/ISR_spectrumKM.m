%INFO: ISR_mian.m的子函数, 用来计算带碰撞(中性)、多种离子成分、带
%   漂移的微观法理论谱。(K&M)
%%----------------------------------------------------------------------%%
% Inputs:
%   Je          - 电子归一化Gordeyev积分
%   Ji          - 离子归一化Gordeyev积分
%   thee        - 电子自变量
%   thei        - 离子自变量
%   k           - 散射差矢 [m^-1]
%   ne          - 电子数密度 [m^-3]
%   ni          - 离子数密度 [m^-3]
%   vTe         - 电子热速度 [m/s]
%   vTi         - 离子热速度 [m/s]
%   he          - 电子德拜半径 [m]
%   hi          - 离子德拜半径 [m]
% Output:
%   spec        - 理论谱 [s]
%%----------------------------------------------------------------------%%
%参考文献
% [1] Kudeki, E., & Milla, M. A. (2011). Incoherent Scatter Spectral
%   Theories-Part I: A General Framework and Results for Small Magnetic
%   Aspect Angles. Ieee Transactions on Geoscience and Remote Sensing,
%   49(1), 315-328. doi:10.1109/Tgrs.2010.2057252
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/08/15,20
%%----------------------------------------------------------------------%%

function spec = ISR_spectrumKM(Je,Ji,thee,thei,k,ne,ni,vTe,vTi,he,hi)

% 热运动引起的电子项散射谱
specTe = 2*ne*real(Je)/k/vTe; % [1] eq(39)
sigma_e = (1 - 1i*thee.*Je)/k^2/he^2; % [1] eq(40)

% 热运动引起的离子项散射谱
specTi = zeros(size(thee));
sigma_i = complex(zeros(size(thee)));
for i=1:length(ni)
    specTi = specTi + 2*ni(i)*real(Ji(i,:))/k/vTi(i); % [1] eq(39)
    sigma_i = sigma_i + (1 - 1i*thei(i,:).*Ji(i,:))/k^2/hi(i)^2; % [1] eq(40)
end

den = abs(1+sigma_i+sigma_e);%分母
spece = power(abs(1+sigma_i)./den,2).*specTe;
speci = power(abs(sigma_e)./den,2).*specTi;
spec = spece + speci; % [1] eq(41)

end