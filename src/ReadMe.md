# ReadMe

Washy 2020-04-30

## ISR函数介绍

### ISR_main.m

- INFO: 理论谱计算主函数
- Needs: `ISR_init.m`; `ISR_updateFactors.m`; 

### ISR_init.m

- INFO: `ISR_main`的子函数，用来初始化参量parameters。
- Needs: `constants.m`; `analysisIon.m`; 
- Use: `parameters = ISR_init(ion,ne,Ti,Tr,percent,frequency,fradar,theta_s);` 

### ISR_updateFactors

- INFO: `ISR_main`的子函数，用来更新parameters中factors。
- Needs: `CollisionCoulomb.m`; 
- Use: `[parameters, mode, theomode] = ISR_updateFactors(parameters,factors);` 

### ISR_updatePlasmas

- INFO: `ISR_main.m`的子函数，用来完善parameters中的plasmas部分。
- Needs: `constants.m;` 
- Use: `parameters = ISR_updatePlasmas(parameters);` 

### ISR_updateDimensionless

- INFO: `ISR_main.m`的子函数，用来更新parameters中dimensionless部分。
- Use: `parameters = ISR_updateDimensionless(parameters);` 

### ISR_gordeyev

- INFO: `ISR_main.m`的子函数，用来计算Gordeyev积分。
- Needs: `fadf.m;` `ISR_SingleParticleACF.m;` `SommerfeldIntegral.m;` 
- Use: `[Je,Ji] = ISR_gordeyev(parameters, mode, theomode);` 

### ISR_SingleParticleACF.m

- INFO: ISR_main.m的子函数, 单粒子ACFs。

### ISR_spectrumKM

- INFO: ISR_mian.m的子函数, 用来计算带碰撞(中性)、多种离子成分、带漂移的微观法理论谱。(K&M)
- Use: `spec = ISR_spectrumKM(Je,Ji,thee,thei,k,ne,ni,vTe,vTi,he,hi);` 

### ISR_spectrumMac

- INFO: ISR_mian.m的子函数, 计算ISR理论谱。(Mac)
- Needs: `fadf.m;` 
- Use: `sigma = ISR_spectrumMac(Je,Ji,omge,omgi,k,ne,he,eta,mu,vTe,vTi,psien,psiin);` 

## 辅助函数介绍

### constants.m

- INFO: 提供各种常数。
- Use: `name = {'c', 'kB'}; res = constants(name);` 

### analysisIon.m

- INFO: 获取指定离子成分的质量和带电量
- Needs: `constants.m`; 
- Use: `ion = {'H+','O+'}; [mi, qi] = analysisIon(ion);` 

### CollisionCoulomb.m

- INFO: 计算库仑碰撞频率
- Use: `[nuii, nuei, nuee] = CollisionCoulomb(ion, Ti, Te, ni, ne);` 

### fadf.m

- [The Voigt/complex error function (second version)](https://ww2.mathworks.cn/matlabcentral/fileexchange/47801-the-voigt-complex-error-function-second-version?focused=9c4306e9-9303-b421-9eeb-10569c46ac8b&tab=function) 

### SommerfeldIntegral.m

- INFO: 计算 exp(-j\*omega\*t)\*f(t)dt 在(a,b)上的积分。
- Use: `[I,b,N] = SommerfeldIntegral(func,a,b,N,omega,funPar,restriction);` 