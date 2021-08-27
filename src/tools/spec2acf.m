%%INFO: �������������渵��Ҷ�任, �������״�Ƶ��תΪʱ��.
%%----------------------------------------------------------------------%%
% Inputs:
%   frequency   - Ƶ�ʷ�Χ [Hz]
%   spectrum    - ������ [s]
% Outputs:
%   tau         - ʱ�� s
%   acf         - ����غ���
%%----------------------------------------------------------------------%%
%�ο�����
% ������ѡ�� ismodel_matlab �е� f_spec2acf.m
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/04/02
%%----------------------------------------------------------------------%%

function [tau, acf] = spec2acf(frequency, spectrum)

Nspec = length(spectrum);
zsize = Nspec*3;
if mod(Nspec,2.0) == 0.0
	zsize = Nspec/2.0;
end

spec = [zeros(1,zsize),spectrum,zeros(1,zsize)];

NFFT = length(spec);
df = frequency(2) - frequency(1);
tau = linspace(-1.0/2/df,1.0/2/df,NFFT);
dtau = tau(2) - tau(1);
acf = fftshift(ifft(ifftshift(spec)))/dtau;

%%
% Nspec = length(spectrum);
% df = frequency(2) - frequency(1);
% tau = linspace(-1.0/2/df,1.0/2/df,Nspec);
% dtau = tau(2) - tau(1);
% acf = fftshift(ifft(ifftshift(spectrum)))/dtau;

end