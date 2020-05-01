%INFO: 计算 exp(-j*omega*t)*f(t)dt 在(a,b)上的积分。
%%----------------------------------------------------------------------%%
% Inputs:
%   func        - 被积函数 f(t)
%   a           - 积分下限
%   b           - 积分上限 [输入值]
%   N           - 分段个数 [输入值]
%   omega       - 频域范围
%   funcPar     - 被积函数自变量
%   restriction - 约束条件
%       eps1    - f(b)必须小于该值 [小量]
%       eps2    - f'(b)必须小于该值 [小量]
%       bm      - 最大积分上限
%       del     - 两次积分结果的差值平方必须小于该值 [小量]
%       Nm      - 最大分段个数，用于防止占用内存过大
%       dif     - f(t)相邻数值差值必须小于该值，用于判断分段个数是否足够
% Output:
%   I           - 积分结果
%   b           - 积分上限 [输出值]
%   N           - 分段个数 [输出值]
%%----------------------------------------------------------------------%%
%参考文献
% [1] Ooi, B. L. (2007). USEFUL INTEGRATION QUADRATURE FOR ELECTROMAGNETICS
%   —THE ERF TRANSFORM. Microwave and Optical Technology Letters, Wiley,
%   2007, 49, 789-791. doi:10.1002/mop.22270
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/26,29
%%----------------------------------------------------------------------%%
%更新内容[Washy 2019/12/02]
%%INFO: 对函数结构进行了部分的优化.
% 1. 调整上限b的确定位置
%%----------------------------------------------------------------------%%

function [I,b,N] = SommerfeldIntegral(func,a,b,N,omega,funPar,restriction)

% restriction = [1e-9,1e-6,1e3,1e-6,512,0.2];

epsFuncb = restriction(1); % f(b) < eps1
bm   = restriction(2);
del  = restriction(3);
Nm   = restriction(4); % N <= Nm

% 判断积分上限是否满足要求
while 1
    if abs(func(b, funPar)) < epsFuncb || b >= bm
        break;
    end
    b = b*2;
end

% 防止出现输入N值过大造成计算量过大
if N > Nm
    N = Nm;
end

% 积分结果缺省值，应远大于积分结果
Ip = 1e10;

while 1
    
    n  = -N:N;
    h  = 1/N*log(1.05*sqrt(2)*N);
    
    An = cosh(n*h).*exp(-sinh(n*h).^2);
    t  = 0.5*(b+a) + 0.5*(b-a)*erf(sinh(n*h));
    
    ft = func(t, funPar);
    
    In = h*(b-a)/sqrt(pi)*An.*ft;
    I  = In*exp(-1i*t'*omega);
    
    % 判断是否I是否合格
    dI = max(abs((I-Ip)./Ip).^2);
    if dI < del || N >= Nm 
        break
    end
    
    Ip = I;
    N  = 2*N;
end

end