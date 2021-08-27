%%INFO: ISR_main.m���Ӻ���, ��������parameters��dimensionless���֡�
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - ������̲������м���� struct
% Output:
%   parameters  - ������̲������м���� struct
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Dougherty, J. P. & Farley, D. T. (1963). A Theory of Incoherent 
%   Scattering of Radio Waves by a Plasma .3. Scattering in a Partly 
%   Ionized Gas. Journal of Geophysical Research, 1963, 68, 5473-5486.
%   doi:10.1029/jz068i019p05473
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%

function parameters = ISR_updateDimensionless(parameters)

k = parameters.radar.k;

%% ����

vTe  = parameters.plasmas.vTe;
rece = 1/k/vTe;

% ���Ӷ�����Ƶ��
parameters.dimensionless.thetae     = rece*parameters.plasmas.omegae; %[1] eq(7)
% ����-������ײƵ��
parameters.dimensionless.psien      = rece*parameters.factors.nuen;
% ���ӹ�һ������Ƶ��
parameters.dimensionless.phie       = rece*parameters.plasmas.Oe;
% ���ӹ�һ��������ײƵ��(ƽ��)
parameters.dimensionless.psiec_par  = rece*parameters.factors.nuei;
% ���ӹ�һ��������ײƵ��(��ֱ)
parameters.dimensionless.psiec_perp = rece*(parameters.factors.nuei+parameters.factors.nuee);

%% ����

vTi  = parameters.plasmas.vTi;
reci = 1/k./vTi;

% ���Ӷ�����Ƶ��
parameters.dimensionless.thetai     = reci'*parameters.plasmas.omegai; %[1] eq(7)
% ����-������ײƵ��
parameters.dimensionless.psiin      = reci.*parameters.factors.nuin;
% ���ӻ���Ƶ��
parameters.dimensionless.phii       = reci.*parameters.plasmas.Oi;
% ���ӿ�����ײƵ��
parameters.dimensionless.psiic      = reci.*parameters.factors.nuii;

end