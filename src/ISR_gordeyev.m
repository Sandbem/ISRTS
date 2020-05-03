%%INFO: ISR_main.m���Ӻ�������������Gordeyev���֡�
%%----------------------------------------------------------------------%%
% Needs: fadf.m; ISR_SingleParticleACF.m; SommerfeldIntegral.m; 
%%----------------------------------------------------------------------%%
% Input:
%   parameters  - ������̲������м���� struct
%   mode        - ģʽ����(��: 0; ��: 1)
%   theomode    - �����׼���ģʽ(Mac: 1; KM: 2)
% Outputs:
%   parameters  - ������̲������м���� struct
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Dougherty, J. P. and Farley, D. T. (1963). A Theory of Incoherent
%   Scattering of Radio Waves by a Plasma .3. Scattering in a Partly
%   Ionized Gas. Journal of Geophysical Research, 68, 5473-5486.
%   doi: 10.1029/jz068i019p05473
% [2] Milla, M. and Kudeki, E. (2009). Particle dynamics description of
%   "BGK collisions" as a Poisson process. Journal of Geophysical
%   Research-Space Physics, 114. doi: 10.1029/2009ja014200
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/30
%%----------------------------------------------------------------------%%

function parameters = ISR_gordeyev(parameters, mode, theomode)

thetae = parameters.dimensionless.thetae;
psien  = parameters.dimensionless.psien;

thetai = parameters.dimensionless.thetai;
psiin  = parameters.dimensionless.psiin;

%% �޴ų�
if mode(3) == 0 && theomode == 2
    Je = sqrt(pi)*fadf(-(thetae - 1i*psien));
    
    [row, column] = size(thetai);
    Ji = complex(zeros(row, column));
    for irow = 1:row
        Ji(irow,:) = sqrt(pi)*fadf(-(thetai(irow,:) - 1i*psiin(irow)));
    end
    
    parameters.spectrum.Je = Je;
    parameters.spectrum.Ji = Ji;
    
	return
end

%% �дų�

func = @ISR_SingleParticleACF;
a = 0;
b = 1;
N = 32;
restriction = [1e-9,512,1e-6,32768];

alpha = parameters.factors.alpha*pi/180;
phie  = parameters.dimensionless.phie;
phii  = parameters.dimensionless.phii;
    
if theomode == 1 % K&M
    psiec_para = parameters.dimensionless.psiec_par;
    psiec_perp = parameters.dimensionless.psiec_perp;
    psiic      = parameters.dimensionless.psiic;
elseif theomode == 2 % Mac
    psiec_para = 0;
    psiec_perp = 0;
    psiic      = 0;
else
    error('���� theomode ���ô���!')
end

funPare = [psien, phie, alpha, psiec_para, psiec_perp];
[Je,~,~] = SommerfeldIntegral(func,a,b,N,thetae,funPare,restriction); %[2] eq(18)

Ji = complex(zeros(size(thetai)));
for ini=1:size(thetai,1)
    funPari = [psiin(ini), phii, alpha, psiic(ini), psiic(ini)];
    Ji(ini,:) = SommerfeldIntegral(func,a,b,N,thetai(ini,:),funPari,restriction); %[2] eq(18)
end

parameters.spectrum.Je = Je;
parameters.spectrum.Ji = Ji;

end