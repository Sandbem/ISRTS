%%INFO: 本函数用来做逆傅里叶变换, 将理论谱从频域转为时域.
%%----------------------------------------------------------------------%%
% Inputs:
%   frequency   - 频率范围 [Hz]
%   spectrum    - 理论谱 [s]
% Outputs:
%   tau         - 时间 s
%   acf         - 自相关函数
%%----------------------------------------------------------------------%%
%参考文献
% 本函数选自 ismodel_matlab 中的 f_spec2acf.m
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