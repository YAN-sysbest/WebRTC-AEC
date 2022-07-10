clc;
clear all;

fid=fopen('near_NLP.pcm', 'rb'); % 麦克风信号——回声+近端
ssin=fread(fid,inf,'int16');
fclose(fid);
%far is speaker played music
fid=fopen('far_NLP.pcm', 'rb'); % 扬声器信号——参考
rrin=fread(fid,inf,'int16');
fclose(fid);
rand('state',13);
fs=16000;
if fs == 8000
cohRange = 2:3;
elseif fs==16000
cohRange = 2;
end

M = 16; % Number of partitions 块数p
N = 64; % Partition length 长度k
L = M*N; % Filter length 滤波器长度

mufb = 0.1;


alp = 0.15; % Power estimation factor 功率估计遗忘因子
alc = 0.1; % Coherence estimation factor 相关估计遗忘因子

%%
threshold=0.5e-3;

len=length(ssin);
w=zeros(L,1); % Sample-based TD(time domain) NLMS
WFb=zeros(N+1,M); % Block-based FD(frequency domain) NLMS
WFbOld=zeros(N+1,M); % Block-based FD NLMS
YFb=zeros(N+1,M);
erfb=zeros(len,1);
edk=zeros(len,1);
eyk=zeros(len,1);

zm=zeros(N,1);
XFm=zeros(N+1,M);
YFm=zeros(N+1,M);
pn0=10*ones(N+1,1);
pn=zeros(N+1,1);
NN=len;
Nb=floor(NN/N)-M;
start=1;
xo=zeros(N,1);
do=xo;
eo=xo;


xfwm=zeros(N+1,M);
dfm=zeros(N+1,M);
WFbD=ones(N+1,1);

cohxdMax = 0;
for kk=1:Nb
    pos = N * (kk-1) + start;  
    %far is speaker played music
    xk = rrin(pos:pos+N-1);
    %near is micphone captured signal
    dk = ssin(pos:pos+N-1);
    
    %far end signal process
    xx = [xo;xk];
    xo = xk;
    tmp = fft(xx);
    XX = tmp(1:N+1);%this is overlap save need, end half needed(because fftshift not used)  
    %near end signal process
    dd = [do;dk]; % Overlap
    do = dk;
    tmp = fft(dd); % Frequency domain
    DD = tmp(1:N+1);
    
    %far end Power estimation
    pn0 = (1 - alp) * pn0 + alp * real(XX.* conj(XX));
    pn = pn0;    
    
    %-------------Filtering
    XFm(:,1) = XX;
    for mm=0:(M-1)
        m=mm+1;
        YFb(:,m) = XFm(:,m) .* WFb(:,m);
    end
    yfk = sum(YFb,2);
    tmp = [yfk ; flipud(conj(yfk(2:N)))];
    ykt = real(ifft(tmp));
    ykfb = ykt(end-N+1:end);
    
    % --------Error estimation
    ekfb = dk - ykfb;
    %For robustness
%     if sum(abs(ekfb)) < sum(abs(dk))
%         ekfb = dk - ykfb;
%     erfb(pos:pos+N-1) = ekfb;
%     else
%         ekfb = dk;
%     erfb(pos:pos+N-1) = dk;
%     end

    erfb(pos:pos+N-1) = ekfb;
    edk(pos:pos+N-1) = dk;
    eyk(pos:pos+N-1) = ykfb;
    tmp = fft([zm;ekfb]); % FD version for cancelling part (overlap-save)
    Ek = tmp(1:N+1);

    % ------------------------ Adaptation
%   Ek2 = Ek ./(M*pn + 0.001); % Normalized error
    Ek2 = Ek ./(pn + 0.001); % Normalized error
    
    absEf = max(abs(Ek2), threshold);
    absEf = ones(N+1,1)*threshold./absEf;
    Ek2 = Ek2.*absEf;

    mEk = mufb.*Ek2;
    PP = conj(XFm).*(ones(M,1) * mEk')';
    tmp = [PP ; flipud(conj(PP(2:N,:)))];
    IFPP = real(ifft(tmp));
    PH = IFPP(1:N,:);
    tmp = fft([PH;zeros(N,M)]);
    FPH = tmp(1:N+1,:);
    WFb = WFb + FPH;

    % Filter update
%    Ek2 = Ek ./(100*pn + 0.001); % Normalized error


% Shift old FFTs
    XFm(:,2:end) = XFm(:,1:end-1);
    YFm(:,2:end) = YFm(:,1:end-1);
    xfwm(:,2:end) = xfwm(:,1:end-1);
    dfm(:,2:end) = dfm(:,1:end-1);


end

out = erfb;
figure(9);
    subplot(3,1,1);plot(ssin);
    subplot(3,1,2);plot(rrin);
    subplot(3,1,3);plot(out);
    
fidout=fopen('E:\Matlab_workspace\myWebRTC-audio-processing\out.pcm','wb');
fwrite(fidout,out,'int16');
fclose(fidout);
hold off
