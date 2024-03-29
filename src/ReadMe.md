# 函数简介

Washy 2020-04-30

## 更新内容

- `2021-04-02` 主函数电子离子温度比`Tr`改为电子温度`Te`，并调换了`Te`和`Ti`的顺序：（1）在设定参数时，使用电子温度作为输入更加的灵活；（2）为了统一，同类型参量，一概按照先电子后离子的顺序。

- `2020-10-15` 修复了`CollisionCoulomb.m`函数中`nuii`计算错误的问题。第123行`nust`计算过程中，需要对先对`ni(:,ia)`进行`reshape`再计算。

- `2020-08-03` （1）对`CollisionCoulomb.m`做了重要更新：a. 可以单独计算电子和离子的库仑碰撞频率；b. 可以一次计算多个密度、温度。（2）多个函数进行了细微调整。
- `2020-05-08` 增加了`gordeyevBessel.m`，并对`ISR_gordeyev.m`的结构进行了调整。
- `2020-04-30` 首次提交。

## ISR函数介绍

### `ISR_main.m` 

- INFO: 理论谱计算主函数
- Needs: `ISR_init.m;` `ISR_updateFactors.m;` `ISR_updatePlasmas.m;` `ISR_updateDimensionless.m;` `ISR_gordeyev.m;` `ISR_spectrumKM.m;` `ISR_spectrumMac.m;` 

### `ISR_init.m` 

- INFO: `ISR_main`的子函数，用来初始化参量parameters。
- Needs: `constants.m`; `analysisIon.m`; 

### `ISR_updateFactors.m` 

- INFO: `ISR_main`的子函数，用来更新parameters中factors。
- Needs: `CollisionCoulomb.m`; 

### `ISR_updatePlasmas.m` 

- INFO: `ISR_main.m`的子函数，用来完善parameters中的plasmas部分。
- Needs: `constants.m;` 

### `ISR_updateDimensionless.m` 

- INFO: `ISR_main.m`的子函数，用来更新parameters中dimensionless部分。

### `ISR_gordeyev.m` 

- INFO: `ISR_main.m`的子函数，用来计算Gordeyev积分。
- Needs: `fadf.m;` `gordeyevSommerfeld.m;` `gordeyevBessel.m;` 

### `ISR_spectrumKM.m` 

- INFO: `ISR_mian.m`的子函数, 用来计算带碰撞(中性)、多种离子成分、带漂移的微观法理论谱。(K&M)

### `ISR_spectrumMac.m` 

- INFO: `ISR_mian.m`的子函数, 计算ISR理论谱。(Mac)
- Needs: `fadf.m;` 

## 辅助函数介绍

### `constants.m` 

- INFO: 提供各种常数。

### `analysisIon.m` 

- INFO: 获取指定离子成分的质量和带电量
- Needs: `constants.m`; 

### `CollisionCoulomb.m` 

- INFO: 计算库仑碰撞频率

### `fadf.m` 

- [The Voigt/complex error function (second version)](https://ww2.mathworks.cn/matlabcentral/fileexchange/47801-the-voigt-complex-error-function-second-version?focused=9c4306e9-9303-b421-9eeb-10569c46ac8b&tab=function) 

### `gordeyevSommerfeld.m` 

- INFO: 使用Sommerfeld积分求解Gordeyev积分。
- Needs: `SingleParticleACF.m;` `SommerfeldIntegral.m;` 

### `SingleParticleACF.m` 

- INFO: 单粒子ACFs。

### `SommerfeldIntegral.m` 

- INFO: 计算 exp(-j\*omega\*t)\*f(t)dt 在(a,b)上的积分。

### `gordeyevBessel.m` 

- INFO: 使用Bessel函数求解Gordeyev积分。
- Needs: `fadf.m;` 