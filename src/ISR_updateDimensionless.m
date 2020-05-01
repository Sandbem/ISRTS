%%INFO: ISR_main.m���Ӻ���, ��������parameters��dimensionless���֡�
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - ������̲������м���� struct
% Output:
%   parameters  - ������̲������м���� struct
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
parameters.dimensionless.thetae     = rece*parameters.plasmas.omegae;
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
parameters.dimensionless.thetai     = reci'*parameters.plasmas.omegai;
% ����-������ײƵ��
parameters.dimensionless.psiin      = reci.*parameters.factors.nuin;
% ���ӻ���Ƶ��
parameters.dimensionless.phii       = reci.*parameters.plasmas.Oi;
% ���ӿ�����ײƵ��
parameters.dimensionless.psiic      = reci.*parameters.factors.nuii;

end