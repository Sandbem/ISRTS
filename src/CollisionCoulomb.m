%%INFO: �������������������ײƵ�ʡ�
%%----------------------------------------------------------------------%%
% Inputs:
%   ion         - ������Ϣ cell
%   Ti          - �����¶� [K]
%   Te          - �����¶� [K]
%   ni          - �������ܶ� [m^-3]
%   ne          - �������ܶ� [m^-3]
% Outputs:
%   vii         - ion-ion������ײƵ�� [Hz]
%   vei         - ele-ion������ײƵ�� [Hz]
%   vee         - ele-ele������ײƵ�� [Hz]
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Rober Schunk and Andrew Nagy (2009). Ionospheres: Physics, Plasma
%   Physics, and Chemistry Second edition.
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/01/02
%%----------------------------------------------------------------------%%
%��������[2019/08/20]
% 1. �޸�ni���뵥λΪm^-3
% 2. �޸�Ti�����벻��ion��Ԫ�ظ����ı�
%%----------------------------------------------------------------------%%
%��������[2019/12/12]
% 1. �޸ĺ�����Ϊ: CollisionCoulomb (ԭ������Ϊ: Coulomb_C)
%%----------------------------------------------------------------------%%
%��������[2020/04/30]
% 1. ȡ����csv�ļ���ȡ, ����IonԪ�����Bst����
%%----------------------------------------------------------------------%%
%����ʾ��[2020/04/30]
% ion = {'O+','H+','NO+'};
% Ti = 1000;
% Te = 1000;
% ne = 1e10;
% ni = [0.5 0.2 0.3]*ne;
%%----------------------------------------------------------------------%%

function [nuii, nuei, nuee] = CollisionCoulomb(ion, Ti, Te, ni, ne)

% ��������
Ion = {'H+','He+','C+','N+','O+','CO+','N2+','NO+','O2+','CO2+'};
% ��ײƵ��ϵ��: �� - s; �� - t
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

% ���㵥λ
ni = ni*1e-6;
ne = ne*1e-6;

% ion_ion��ײƵ��
[~, ia, ib] = intersect(ion, Ion, 'stable');
nust = Bst(ib, ib).*repmat(ni(ia), length(ia), 1)/Ti^1.5; %[1] eq(4.143)
nuii = zeros(1, length(ion));
nuii(ia) = sum(nust, 2);

% ele_ion��ײƵ��
% vei = 54.5*ni*Z'.^2/power(Te, 1.5); %[1] eq(4.144)
nuei = 54.5*ne/power(Te, 1.5);

% ele_ele��ײƵ��
nuee = 54.5/sqrt(2)*ne/power(Te, 1.5); %[1] eq(4.145)

end