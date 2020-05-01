%INFO: ���� exp(-j*omega*t)*f(t)dt ��(a,b)�ϵĻ��֡�
%%----------------------------------------------------------------------%%
% Inputs:
%   func        - �������� f(t)
%   a           - ��������
%   b           - �������� [����ֵ]
%   N           - �ֶθ��� [����ֵ]
%   omega       - Ƶ��Χ
%   funcPar     - ���������Ա���
%   restriction - Լ������
%       eps1    - f(b)����С�ڸ�ֵ [С��]
%       eps2    - f'(b)����С�ڸ�ֵ [С��]
%       bm      - ����������
%       del     - ���λ��ֽ���Ĳ�ֵƽ������С�ڸ�ֵ [С��]
%       Nm      - ���ֶθ��������ڷ�ֹռ���ڴ����
%       dif     - f(t)������ֵ��ֵ����С�ڸ�ֵ�������жϷֶθ����Ƿ��㹻
% Output:
%   I           - ���ֽ��
%   b           - �������� [���ֵ]
%   N           - �ֶθ��� [���ֵ]
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Ooi, B. L. (2007). USEFUL INTEGRATION QUADRATURE FOR ELECTROMAGNETICS
%   ��THE ERF TRANSFORM. Microwave and Optical Technology Letters, Wiley,
%   2007, 49, 789-791. doi:10.1002/mop.22270
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/26,29
%%----------------------------------------------------------------------%%
%��������[Washy 2019/12/02]
%%INFO: �Ժ����ṹ�����˲��ֵ��Ż�.
% 1. ��������b��ȷ��λ��
%%----------------------------------------------------------------------%%

function [I,b,N] = SommerfeldIntegral(func,a,b,N,omega,funPar,restriction)

% restriction = [1e-9,1e-6,1e3,1e-6,512,0.2];

epsFuncb = restriction(1); % f(b) < eps1
bm   = restriction(2);
del  = restriction(3);
Nm   = restriction(4); % N <= Nm

% �жϻ��������Ƿ�����Ҫ��
while 1
    if abs(func(b, funPar)) < epsFuncb || b >= bm
        break;
    end
    b = b*2;
end

% ��ֹ��������Nֵ������ɼ���������
if N > Nm
    N = Nm;
end

% ���ֽ��ȱʡֵ��ӦԶ���ڻ��ֽ��
Ip = 1e10;

while 1
    
    n  = -N:N;
    h  = 1/N*log(1.05*sqrt(2)*N);
    
    An = cosh(n*h).*exp(-sinh(n*h).^2);
    t  = 0.5*(b+a) + 0.5*(b-a)*erf(sinh(n*h));
    
    ft = func(t, funPar);
    
    In = h*(b-a)/sqrt(pi)*An.*ft;
    I  = In*exp(-1i*t'*omega);
    
    % �ж��Ƿ�I�Ƿ�ϸ�
    dI = max(abs((I-Ip)./Ip).^2);
    if dI < del || N >= Nm 
        break
    end
    
    Ip = I;
    N  = 2*N;
end

end