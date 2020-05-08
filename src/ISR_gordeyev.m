%%INFO: ISR_main.m的子函数，用来计算Gordeyev积分。
%%----------------------------------------------------------------------%%
% Needs: fadf.m; gordeyevSommerfeld.m; gordeyevBessel.m; 
%%----------------------------------------------------------------------%%
% Input:
%   parameters  - 计算过程产生的中间变量 struct
%   mode        - 模式控制(关: 0; 开: 1)
%   gordmode    - Gordeyev积分模式(Sommerfeld: 1; Bessel: 2)
%   theomode    - 理论谱计算模式(Mac: 1; KM: 2)
% Outputs:
%   parameters  - 计算过程产生的中间变量 struct
%%----------------------------------------------------------------------%%
%参考文献
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
%更新内容[Washy 2020/05/08]
%INFO: 重大更新, 对函数整体进行了重构。
%%----------------------------------------------------------------------%%

function parameters = ISR_gordeyev(parameters,mode,gordmode)

thetae = parameters.dimensionless.thetae;
psien  = parameters.dimensionless.psien;

thetai = parameters.dimensionless.thetai;
psiin  = parameters.dimensionless.psiin;

%% 无磁场
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

%% 有磁场

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
    error('参数 gordmode 设置错误!')
end

parameters.spectrum.Je = Je;
parameters.spectrum.Ji = Ji;

end