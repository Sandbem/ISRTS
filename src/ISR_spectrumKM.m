%INFO: ISR_mian.m���Ӻ���, �����������ײ(����)���������ӳɷ֡���
%   Ư�Ƶ�΢�۷������ס�(K&M)
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - ������̲������м���� struct
% Output:
%   spec        - ������ [s]
%   parameters  - ������̲������м���� struct
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Kudeki, E., & Milla, M. A. (2011). Incoherent Scatter Spectral
%   Theories-Part I: A General Framework and Results for Small Magnetic
%   Aspect Angles. Ieee Transactions on Geoscience and Remote Sensing,
%   49(1), 315-328. doi:10.1109/Tgrs.2010.2057252
% [2] Milla, M. and Kudeki, E. (2009). Particle dynamics description of
%   "BGK collisions" as a Poisson process. Journal of Geophysical
%   Research-Space Physics, 114. doi: 10.1029/2009ja014200
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/08/15,20
%%----------------------------------------------------------------------%%
%%��������[Washy 2020/05/03]
% 1. ��������: parameters
% 2. �������: parameters
%%----------------------------------------------------------------------%%

function [spec,parameters] = ISR_spectrumKM(parameters)

k    = parameters.radar.k;
ne   = parameters.plasmas.ne;
vTe  = parameters.plasmas.vTe;
he   = parameters.plasmas.he;
ni   = parameters.plasmas.ni;
vTi  = parameters.plasmas.vTi;
hi   = parameters.plasmas.hi;
thee = parameters.dimensionless.thetae;
psien= parameters.dimensionless.psien;
thei = parameters.dimensionless.thetai;
psiin= parameters.dimensionless.psiin;
Ge   = parameters.spectrum.Je;
Gi   = parameters.spectrum.Ji;

% ���˶�����ĵ�����ɢ����
Je = Ge./(1 - psien*Ge); %[2] eq(20)
specTe = 2*ne*real(Je)/k/vTe; %[1] eq(39)
sigma_e = (1 - 1i*thee.*Je)/k^2/he^2; %[1] eq(40)

% ���˶������������ɢ����
specTi = zeros(size(thee));
sigma_i = complex(zeros(size(thee)));
for i=1:length(ni)
    Ji = Gi(i,:)./(1 - psiin(i)*Gi(i,:)); %[2] eq(20)
    specTi = specTi + 2*ni(i)*real(Ji)/k/vTi(i); %[1] eq(39)
    sigma_i = sigma_i + (1 - 1i*thei(i,:).*Ji)/k^2/hi(i)^2; %[1] eq(40)
end

den = abs(1+sigma_i+sigma_e); % ��ĸ
parameters.spectrum.spece = power(abs(1+sigma_i)./den,2).*specTe;
parameters.spectrum.speci = power(abs(sigma_e)./den,2).*specTi;
spec = parameters.spectrum.spece + parameters.spectrum.speci; %[1] eq(41)

end