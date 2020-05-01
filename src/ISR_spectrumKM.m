%INFO: ISR_mian.m���Ӻ���, �����������ײ(����)���������ӳɷ֡���
%   Ư�Ƶ�΢�۷������ס�(K&M)
%%----------------------------------------------------------------------%%
% Inputs:
%   Je          - ���ӹ�һ��Gordeyev����
%   Ji          - ���ӹ�һ��Gordeyev����
%   thee        - �����Ա���
%   thei        - �����Ա���
%   k           - ɢ���ʸ [m^-1]
%   ne          - �������ܶ� [m^-3]
%   ni          - �������ܶ� [m^-3]
%   vTe         - �������ٶ� [m/s]
%   vTi         - �������ٶ� [m/s]
%   he          - ���ӵ°ݰ뾶 [m]
%   hi          - ���ӵ°ݰ뾶 [m]
% Output:
%   spec        - ������ [s]
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Kudeki, E., & Milla, M. A. (2011). Incoherent Scatter Spectral
%   Theories-Part I: A General Framework and Results for Small Magnetic
%   Aspect Angles. Ieee Transactions on Geoscience and Remote Sensing,
%   49(1), 315-328. doi:10.1109/Tgrs.2010.2057252
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/08/15,20
%%----------------------------------------------------------------------%%

function spec = ISR_spectrumKM(Je,Ji,thee,thei,k,ne,ni,vTe,vTi,he,hi)

% ���˶�����ĵ�����ɢ����
specTe = 2*ne*real(Je)/k/vTe; % [1] eq(39)
sigma_e = (1 - 1i*thee.*Je)/k^2/he^2; % [1] eq(40)

% ���˶������������ɢ����
specTi = zeros(size(thee));
sigma_i = complex(zeros(size(thee)));
for i=1:length(ni)
    specTi = specTi + 2*ni(i)*real(Ji(i,:))/k/vTi(i); % [1] eq(39)
    sigma_i = sigma_i + (1 - 1i*thei(i,:).*Ji(i,:))/k^2/hi(i)^2; % [1] eq(40)
end

den = abs(1+sigma_i+sigma_e);%��ĸ
spece = power(abs(1+sigma_i)./den,2).*specTe;
speci = power(abs(sigma_e)./den,2).*specTi;
spec = spece + speci; % [1] eq(41)

end