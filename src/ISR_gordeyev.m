%%INFO: ISR_main.m���Ӻ�������������Gordeyev���֡�
%%----------------------------------------------------------------------%%
% Needs: fadf.m; gordeyevSommerfeld.m; gordeyevBessel.m; 
%%----------------------------------------------------------------------%%
% Input:
%   parameters  - ������̲������м���� struct
%   mode        - ģʽ����(��: 0; ��: 1)
%   gordmode    - Gordeyev����ģʽ(Sommerfeld: 1; Bessel: 2)
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
%��������[Washy 2020/05/08]
%INFO: �ش����, �Ժ�������������ع���
%%----------------------------------------------------------------------%%

function parameters = ISR_gordeyev(parameters,mode,gordmode)

thetae = parameters.dimensionless.thetae;
psien  = parameters.dimensionless.psien;

thetai = parameters.dimensionless.thetai;
psiin  = parameters.dimensionless.psiin;

%% �޴ų�
if mode(3) == 0
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

alpha      = parameters.factors.alpha*pi/180;
phie       = parameters.dimensionless.phie;
phii       = parameters.dimensionless.phii;
    
if gordmode == 1 % SommerfeldIntegral
    
    psiec_para = parameters.dimensionless.psiec_par;
    psiec_perp = parameters.dimensionless.psiec_perp;
    psiic      = parameters.dimensionless.psiic;
    
    Je = gordeyevSommerfeld(thetae,psien,alpha,phie,psiec_para,psiec_perp);
    
    Ji = complex(zeros(size(thetai)));
    for ini=1:size(thetai,1)
        Ji(ini,:) = gordeyevSommerfeld(thetai(ini,:),psiin(ini),alpha,phii(ini),psiic(ini), psiic(ini));
    end
    
elseif gordmode == 2 % BesselFunction
    
    Je = gordeyevBessel(thetae,psien,alpha,phie,128);
    
    Ji = complex(zeros(size(thetai)));
    for ini=1:size(thetai,1)
        Ji(ini,:) = gordeyevBessel(thetai(ini,:),psiin(ini),alpha,phii(ini),128);
    end
    
else
    error('���� gordmode ���ô���!')
end

parameters.spectrum.Je = Je;
parameters.spectrum.Ji = Ji;

end