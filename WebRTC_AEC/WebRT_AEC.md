## WebRTC_AEC回声消除

### 1.算法整理框图

使用该AEC算法要注意两点：

1）延时要小，因为算法默认滤波器长度是分为12块，每块64点，按照8000采样率，也就是12*8ms=96ms的数据，而且超过这个长度是处理不了的。

2）延时抖动要小，因为算法是默认10块也计算一次参考数据的位置（即滤波器能量最大的那一块），所以如果抖动很大的话找参考数据时不准确的，这样回声就消除不掉了。

对fullaec.m文件里的每一行都做了注释，红色部分是第一次看的疑问，绿色部分是看完几遍之后的解答，最后还有几个红色问题依旧不太清楚

<img src="C:\Users\Admin\Desktop\简历\项目笔记\笔记图\aec.jpg" alt="aec" style="zoom:200%;" />



### 2. 延迟估计

为了估计近端数据相对于远端数据的延迟，算法需要在缓存中缓存一定时长的远端历史数据，尽可能覆盖到出现的最大延迟，然后用近端数据和远端缓存中的每一帧历史数据进行匹配，找到相似度最高的那一帧远端数据，两者时间索引的差值就是估计出来的延迟。

相较于传统互相关算法具有计算复杂度低、准确度高的优点。

**实现流程**：

- 将参与运算的128点远端和近端信号变换到频域，获得65点的频域信息。
- 取第12个频点到43个频点之间的32个子带数据，更新出每个频点对应的阈值。
- 将32个子带数据与阈值做比较，大于阈值记为1，小于记为0，这样可以用“二元谱”来表示。
- 将远端信号二元谱集体后移一位，把最新的二元谱放在缓存中最新的位置。
- 将近端信号和远端信号的二元谱按位异或，统计1的个数，其越小表示越相似，二者差值则为延迟。

苹果手机延迟基本120ms，Andriod手机延迟基本200ms-300ms

### 3.自适应滤波（PBFDAF算法）

```matlab
% 初始化变量，读取远端数据xk，近端数据dk
% 拼凑远端信号，前一帧与当前帧拼接
xx = [xo;xk];
tmp = fft(xx); XX=tmp(1:N+1);
% 拼凑近端信号，前一帧与当前帧拼接
dd = [do;dk];
tmp = fft(dd); DD=tmp(1:N+1);
% 计算远端功率谱，用于更新公式归一化参考信号（远端数据）
pn0 = (1-alp)*pn0+alp*real(XX.*conj(XX));
% 是否开启舒适降噪
if(CNon) end
% 当前帧赋值，取16帧数据，为65×16矩阵
XFm(:,1) = XX;
% 进行滤波输出
YFb(:,m) = XFm(:,m).*WFm(:,m);
% 将16块信息按行累加，得1块65×1的矩阵
yfk = sum(YFb,2);
% 取滤波后信号的反变换，再取后N个点作为当前一帧输出
ykt = real(ifft(tmp));ykfb = ykt(end-N+1:end);
% 计算当前帧误差信号
ekfb = dk -ykfb;


% 将时域的误差信号前面补64个0，做FFT变换，得到频域误差Ek
tmp = fft([zm;ekfb]);Ek = tmp(1;N+1);
% 做误差信号的归一化，如同NLMS归一化公式
EK2 = Ek./(pn + 0.001);
% 对Ek2做限幅处理，比门限大取门限，比门限小保持不变
absEf = max(abs(Ek2),threshold); absEf = ones(N+1,1)*threshold./absEf; Ek2 = Ek2.*absEf;
% 做收敛速度的选取，mufb做步长取0.5
mEk = mufb.*Ek2;
% PP是16块远端信号频谱共轭乘以误差信号频谱，用于更新公式的预调整量，做处理后可直接叠加到WFb上，
PP = conj(XFm).*(ones(M,1)*mEk')';
% 需要将PP变换到时域上，将后面64点置零，再做FFT，得到用于叠加的调整量，目的是为了避免线性卷积变成循环卷积
% 最后叠加更新
WFb = WFb + FPH;  
```



### 4.非线性处理（NLP）

```matlab
% 每过几帧之后，计算WFbD权重最大的那个块dIdx的索引
if  mod(kk,10*mult)==0 end

% 进行NLP处理

% 拼凑误差信号，前一帧与当前帧拼接
ee = [eo;ekfb];
% 把xx、dd、ee加窗后做FFT变换
fft();
% 取最近16块中dIdx所对应的那一块传给xf用做后续计算
xfwm(:,1) = xf; xf = xfwm(:,dIdx);
% 计算ef、df、xf的功率谱
Se = gamma*Se + (1-gamma)*real(ef.*conj(ef));
% 计算互功率谱，计算远端信号和近端信号互功率谱Sxd，误差信号和近端信号的互功率谱Sed
Sxd = gamma*Sxd + (1-gamma)*real(xf.*conj(df));
% 计算误差信号和近端信号的相关性cohed，远端信号和近端信号相关性cohxd
cohed = real(Sed.*conj(Sed))./(Se.*Sd+1e-10);
cohxd = real(Sxd.*conj(Sxd))./(Sx.*Sd+1e-10);
% 取1-cohxd和cohed的最小值，表示该变量越大，回声越小
hnled = min(1-cohxd,cohed);
% 计算cohed和1-cohxd在echoBandRange范围的平均相干性
cohedMean = mean(cohed(echoBandRange));
hnlSortQ = mean(1-cohxd(echoBandRange));
% 对1-cohxd进行升序排序，对hnled进行升序排序
[hnlSort,hnlSortIdx] = sort(1-cohxd(echoBandRange));
[hnlSort2,hnlSortIdx2] = sort(hnled(echoBandRange));
% 计算两个hnled的值，hnlPrefAvg在排序后的3/4处；hnlPrefAvgLow在排序后的1/2处

% 判断抑制等级
if cohedMean > 0.98 & hnlSortQ > 0.9 suppState = 0;
elseif cohedMean < 0.95 | hnlSortQ < 0.8 suppState = 1;
end
% 如果suppState = 0不需要进行回声抑制，设置
hnled = cohed; hnlPrefAvg = cohedMean; hnlPrefAvgLow = cohedMean;
% 判断hnlSortQ是否小于前一次的不相关性，且小于0.75，若一次都不满足，则说明回声很小，要使用较小的ovrd值和较大的hnled值，避免发生过抑制
ovrd = 3;hnled = 1-cohxd; hnlPrefAvg = hnlSort; hnlPrefAvgLow = hnlSort;
% 判hnlPrefAvgLow是否小于之前的最小值，且小于0.6，连续两次循环满足先 更新hnlLocalMin，设置抑制等级ovrd
ovrd =  ;
% 对hnlLocalMin做处理，防止大于1
% 平滑更新ovrdSm，使ovrdSm是一个较大的值，为了尽量抑制回声
ovrdSm = 0.9*ovrdSm + 0.1*ovrd;

% 做发散处理，误差能量大于近端能量，则用近端频谱更新误差频谱，设置发散状态1；误差能量的1.05倍小于近端能量，设置发散状态0
% 误差能量大于近端能量的19.95倍，滤波器发散，重新置零


% 计算权重曲线wCurve
wCurve = [0;0.3*sqrt(linspace(0,1,N))'+0.1];
% 平滑hnled，wCurve分布让频率高的本次hnled占比大，频率低的上一次hnled平滑结果占比大;min(hnlPrefAvg,hnled)为hnled不超过3/4最大值
hnled = wCurve.*min(hnlPrefAvg,hnled)+(1-weight).*hnled;
% 用产生用于更新的幂指数 
od = ovrdSm*（sqrt(linspace(0,1,N+1))+1）;
hnled =hnled.^(od.*sshift);
% 用hnled系数预误差频谱频域相乘，对ef中的残余回声进一步抑制
ef = ef.*(hnled);

% 是否开启舒适降噪
if(CNon) end
% 使用重叠相加法获得时域平滑信号
mixw = real(ifft(ef));
mola = mbuf(end-N+1:end) + mixw(1:N);

```



## Speex_AEC回声消除算法

### 1.变步长计算

Speex的AEC是以NLMS为基础，采用时变的步长因子。

根据算法原理，在双降时，Speex会计算出一个很小的步长因子，这使得滤波器系数在双降阶段能够保持稳定，避免发散。因为双降时，回声中混杂着近端人说话的声音，这时近端信号不能作为自适应滤波器的期望输出，此时为了保持滤波器的稳定，停止系数更新。而在远端单讲时，算法会将步长因子动态调大，使得滤波器能够快速收敛。

算法用MDF频域实现，最终推导出最优步长估计：**残余回声与误差之比。最优步长等于残余回声方差与误差信号方差之比**

### 2.MDF滤波器

（1）将输入信号分块处理，使用Overlap-and-Save方法计算卷积的分块结果。

（2）分块卷积转入频域计算，使用FFT将计算复杂度从O(N^2)降到O(Nlog_2N)。

（3）进一步将FIR滤波器系数进行分段，一组卷积分解为多组卷积之和。输入和输出的分块长度变短，大大减小滤波器时延.

频域LMS自适应滤波：[](https://www.jianshu.com/p/e4ee7b6496e1)

### 3.双线性滤波器结构

Speex 采用了双线性滤波器的结构，分为迭代更新的背景滤波器和一个非自适应的前景滤波器。

背景滤波器按照正常的频域分块NLMS算法做自适应更新，其误差信号不作为系统输出。对于前景滤波器，其滤波器系数是在一定条件下复制背景滤波器所得，前景滤波器不做自适应，但误差信号会作为系统输出。

算法在每次循环时都会比较更新后的背景滤波器和当前的前景滤波器抑制回声的能力，如果背景滤波器的抑制能力强，就会将背景滤波器系数复制给前景滤波器。此外，如果背景滤波器发散，也会将前景滤波器的值再赋值回背景滤波器。双降滤波器的结构设计使得算法鲁棒性更强，在处理双讲时，滤波器更稳定。

![img](https://upload-images.jianshu.io/upload_images/9127311-e6586b687b9f7796.png?imageMogr2/auto-orient/strip|imageView2/2/format/webp)

## Speex与WebRtc回声消除小节

回声消除AEC包含：  延时估计对齐+线性自适应滤波器+NLP(双讲检测、处理)+舒适噪声CNG

### 1.speex AEC

（1）没有NLP

（2）只考虑实时DSP系统，即是没有延时对齐等

（3）自适应滤波（MDF）使用双滤波器结构，自适应滤波器因子自动更新

### 2.WebRtc AEC

（1）双讲检测没有，双讲时远端的声音会消没了

（2）自适应滤波（PBFDAF），固定自适应因子 0.6

（3）抑制是使用相关性技术，近端误差，近端远端，由低频段相关性参数求出gain值

### 3.speex 与WebRtc 效果对比

对于aec，webrtc主要依赖NLP，speex主要是自适应滤波器（双滤波器）

实际效果对比：如果样本非线性不严重，两者的效果都不错；对于非线性speex效果就很差了，webrtc的效果好；双讲时，webrtc出来的音质就很差，有吃音现象。至于webrtc的aecm音质差，单讲会有吱吱声。

优化点：可对webrtc的aec加入双讲检测，双讲处理。
