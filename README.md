# ISRTS

Washy 2020-05-01

## 功能

- 计算非相干散射理论谱：多种离子成分、漂移、碰撞、磁场

## 如何使用

共9个输入参量，分别为：

- `ion`：离子成分，元胞数组。例：`ion = {'O+'};`或`ion = {'O+','H+'};` 
- `ne`：电子数密度，单位$m^{-3}$。例：`ne = 5e10;` 
- `Ti`：离子温度，单位$K$。例：`K = 1000;` 
- `Tr`：电子离子温度比，$Tr = T_e/T_i$。例：`Tr = 1;` 
- `percent`：离子各成分比例，行向量，需要和`ion`的长度相同且各项之和为1。例：`percent = 1;`或`percent = [0.8, 0.2];` 
- `frequency`：多普勒频率，单位Hz。例：`frequency = linspace(-100e3, 100e3, 1000);` 
- `fradar`：雷达频率，单位Hz。例：`fradar = 440e6;` 
- `theta`：散射角，单位度。例：`theta = 180;` 
- `factors`：其他补充因素参数，结构体，包含以下各参量
  - `factors.ud`：漂移速度，单位m/s，一维数组，`[ude, udi]`，第1个参数为电子漂移速度，第2-end参数为离子漂移速度。例：`factors.ud = [0, 0];`或`factos.ud = [0, 0, 0];` 
  - `factors.nu`：与中性成分的碰撞频率，单位Hz，一维数组，`[nuen, nuin]`，第1个参数为电子碰撞频率，第2-end参数为离子碰撞频率。例：`factors.nu = [0, 0];` 或 `factors.nu = [0, 0, 0];` 
  - `factors.B`：地磁场参数，一维数组，单位T，`[B0, alpha]`，第1个参数为地磁场磁感应强度，第二个参数为散射差矢与磁场的夹角，`alpha = 90`表示垂直于地磁场。例：`factors.B = [3.5e-5, 60];` 
  - `factors.mode`：模式控制，1×4数组，0：关，1：开，从前至后依次控制`ud`，`nu`，`B`以及库仑碰撞的开关。例：`factors.mode = [0, 0, 0, 0];` 
  - `factors.gordmode`：Gordeyev积分模式选择，1：Sommerfeld积分求解，2：Bessel函数展开式求解。例：`factors.gordmode = 1;` 
  - `factors.theomode`：理论公式选择，1：使用Kudeki&Milla（2006,2011）公式；2：使用Farley等公式。

主函数调用方法：

设定好以上各参数之后，使用语句：

`[spec,parameters]=ISR_main(ion,ne,Ti,Tr,percent,frequency,fradar,theta);` 

或

`[spec,parameters]=ISR_main(ion,ne,Ti,Tr,percent,frequency,fradar,theta,factors);` 

参数`factors`不输入时，使用默认设置，即不考虑factors中的所有参数影响。