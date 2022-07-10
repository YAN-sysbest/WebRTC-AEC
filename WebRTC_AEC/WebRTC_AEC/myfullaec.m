clc;
clear all;
% Partitioned block frequency domain adaptive filtering NLMS and
% standard time-domain sample-based NLMS
%near is micphone captured signal
fid=fopen('near_NLP.pcm', 'rb'); % Load far end 麦克风信号——回声+近端
ssin=fread(fid,inf,'int16');
fclose(fid);
%far is speaker played music
fid=fopen('far_NLP.pcm', 'rb'); % Load fnear end 扬声器信号——参考
rrin=fread(fid,inf,'int16');
fclose(fid);
rand('state',13);
fs=16000;
mult=fs/8000;

% Flags
NLPon=1; % NLP on
CNon=1; % Comfort noise on

M = 16; % Number of partitions
N = 64; % Partition length
L = M*N; % Filter length

mufb = 0.01;

alp = 0.15; % Power estimation factor 
% alc = 0.1; % Coherence estimation factor
%% Changed a little %%
step = 0.1875;%0.1875; % Downward step size
%%
threshold=0.5e-3;
echoBandRange = ceil(50/fs*N):floor(3000*2/fs*N);

suppState = 1;
ramp = 1.0003; % Upward ramp
rampd = 0.999; % Downward ramp

len=length(ssin);
xk=zeros(N,1);
dk=zeros(N,1);
% w=zeros(L,1); % Sample-based TD(time domain) NLMS
WFb=zeros(N+1,M); % Block-based FD(frequency domain) NLMS
% WFbOld=zeros(N+1,M); % Block-based FD NLMS
YFb=zeros(N+1,M);
erfb=zeros(len,1);
% erfb3=zeros(len,1);

ercn=zeros(len,1);
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

hnledAvg=zeros(Nb+1,1);
hnlxdAvg=zeros(Nb+1,1);
ovrdV=zeros(Nb+1,1);
dIdxV=zeros(Nb+1,1);
hnlSortQV=zeros(Nb+1,1);
hnlPrefAvgV=zeros(Nb+1,1);
hnled = zeros(N+1, 1);
weight=zeros(N+1,1);

xfwm=zeros(N+1,M);
dfm=zeros(N+1,M);
WFbD=ones(N+1,1);

hnlLocalMin = 1;
cohxdLocalMin = 1;
hnlLocalMinV=zeros(Nb+1,1);
cohxdLocalMinV=zeros(Nb+1,1);
hnlMinV=zeros(Nb+1,1);
dkEnV=zeros(Nb+1,1);
ekEnV=zeros(Nb+1,1);
ovrd = 2;
ovrdSm = 2;
hnlMin = 1;
minCtr = 0;
dIdx = 1;
hnlMinCtr = 0;
hnlNewMin = 0;
divergeState = 0;

Sy=ones(N+1,1);
Sym=1e7*ones(N+1,1);
wins=[0;sqrt(hanning(2*N-1))];
mbuf=zeros(2*N,1);

cohxd = zeros(N+1,1);
Se = zeros(N+1,1);
Sd = zeros(N+1,1);
Sx = zeros(N+1,1);
Sed = zeros(N+1,1);
Sxd = zeros(N+1,1);

delay = 0;
for kk=1:Nb - 1
    pos = N * (kk-1) + start; 
    %far is speaker played music
    xk = rrin(pos:pos+N-1);
    %near is micphone captured signal
    dk = ssin(pos:pos+N-1);
    %延迟位置设置（模拟dk+延迟xk）
%     pos_delay = N * (kk-1) - delay + 1;
%     if pos_delay >= 1
%         dk = ssin(pos_delay:pos_delay+N-1);
%         xk = rrin(pos_delay:pos_delay+N-1);
%     end
    
    
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
    % 开启舒适降噪
    if(CNon)
        Yp = real(conj(DD).*DD);
        Sy = (1 - alp)*Sy + alp*Yp;
        mm = min(Sy,Sym);
        diff = Sym - mm;
%         if(kk > 50)
%             Sym = (mm + step*diff)*ramp;
%         end
        Sym = (mm + step*diff)*ramp;
    end
    %-------------Filtering
    XFm(:,1) = XX;
    for mm=1:M
        YFb(:,mm) = XFm(:,mm) .* WFb(:,mm);
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

    if mod(kk, 10*mult) == 0
        WFbEn = sum(real(WFb.*conj(WFb)));
        %WFbEn = sum(abs(WFb));
        [tmp, dIdx] = max(WFbEn);
        
        WFbD = sum(abs(WFb(:, dIdx)),2);
        WFbD = min(max(WFbD, 0.5), 4);
    end
    dIdxV(kk) = dIdx;
    
%     dIdx = 2;
    
    % NLP
    if (NLPon)  
        ee = [eo;ekfb];
        eo = ekfb;
        window = wins;
        if fs == 8000
            gamma = 0.9;
        else 
            gamma = 0.93;
        end
        tmp = fft(xx.*window);
        xf = tmp(1:N+1);
        tmp = fft(dd.*window);
        df = tmp(1:N+1);
        tmp = fft(ee.*window);
        ef = tmp(1:N+1);

        xfwm(:,1) = xf;
        xf = xfwm(:,dIdx);
        dfm(:,1) = df;
        
        SxOld = Sx;

        Se = gamma*Se + (1-gamma)*real(ef.*conj(ef));
        Sd = gamma*Sd + (1-gamma)*real(df.*conj(df));
        Sx = gamma*Sx + (1 - gamma)*real(xf.*conj(xf));

        Sxd = gamma*Sxd + (1 - gamma)*xf.*conj(df);
        Sed = gamma*Sed + (1-gamma)*ef.*conj(df);

        cohed = real(Sed.*conj(Sed))./(Se.*Sd + 1e-10);
        cohxd = real(Sxd.*conj(Sxd))./(Sx.*Sd + 1e-10);
        
%         freqSm = 0.55;
%         cohxd(2:end) = filter(freqSm, [1 -(1-freqSm)], cohxd(2:end));
%         cohxd(end:2) = filter(freqSm, [1 -(1-freqSm)], cohxd(end:2));

        hnled = min(1 - cohxd, cohed);

        cohedMean = mean(cohed(echoBandRange));
        [hnlSort, 	hnlSortIdx] = sort(1-cohxd(echoBandRange));
        [xSort, xSortIdx] = sort(Sx);

        hnlSortQ = mean(1 - cohxd(echoBandRange));
        %hnlSortQ = mean(1 - cohxd);
        [hnlSort2, hnlSortIdx2] = sort(hnled(echoBandRange));
        %[hnlSort2, hnlSortIdx2] = sort(hnled);
        
        hnlQuant = 0.75;
        hnlQuantLow = 0.5;
        qIdx = floor(hnlQuant*length(hnlSort2));
        qIdxLow = floor(hnlQuantLow*length(hnlSort2));
        hnlPrefAvg = hnlSort2(qIdx);
        hnlPrefAvgLow = hnlSort2(qIdxLow);

%         suppState = 1;
        if cohedMean > 0.98 & hnlSortQ > 0.9
            suppState = 0;
        elseif cohedMean < 0.95 | hnlSortQ < 0.8
            suppState = 1;
        end
        if hnlSortQ < cohxdLocalMin & hnlSortQ < 0.75
            cohxdLocalMin = hnlSortQ;
        end
        
        if cohxdLocalMin == 1 %此时回声小，抑制等级低
            ovrd = 3;
            hnled = 1 - cohxd;
            hnlPrefAvg = hnlSortQ;
            hnlPrefAvgLow = hnlSortQ;
        end

        if suppState == 0 %此时几乎没有回声，不需要抑制
            hnled = cohed;
            hnlPrefAvg = cohedMean;
            hnlPrefAvgLow = cohedMean;
        end

        %if hnlPrefAvg < hnlLocalMin & hnlPrefAvg < 0.6
        if hnlPrefAvgLow < hnlLocalMin & hnlPrefAvgLow < 0.6
            %hnlLocalMin = hnlPrefAvg;
            %hnlMin = hnlPrefAvg;
            hnlLocalMin = hnlPrefAvgLow;
            hnlMin = hnlPrefAvgLow;
            hnlNewMin = 1;
            hnlMinCtr = 0;
        end
        if hnlNewMin == 1
            hnlMinCtr = hnlMinCtr + 1;
        end
        if hnlMinCtr == 2
            hnlNewMin = 0;
            hnlMinCtr = 0;
%             ovrd = max(log(0.0001)/log(hnlMin), 2);
%             ovrd = max(log(0.0001)/(log(hnlMin + 1e-10) + 1e-10), 5);
%             ovrd = max(log(0.00001)/(log(hnlMin + 1e-10) + 1e-10), 5);
%             ovrd = max(log(0.0001)/log(hnlPrefAvg), 2);
            ovrd = max(log(0.001)/log(hnlMin), 8);
%             ovrd = max(log(0.00000001)/(log(hnlMin + 1e-10) + 1e-10), 3);
        end
        hnlLocalMin = min(hnlLocalMin + 0.0008/mult, 1);
        cohxdLocalMin = min(cohxdLocalMin + 0.0004/mult, 1);

        if ovrd < ovrdSm
            ovrdSm = 0.99*ovrdSm + 0.01*ovrd;
        else
            ovrdSm = 0.9*ovrdSm + 0.1*ovrd;
        end

%         ovrdSm = 2;
        ekEn = sum(Se);
        dkEn = sum(Sd);

        %发散处理
        if divergeState == 0
            if ekEn > dkEn
                ef = df;
                divergeState = 1;
            end
        else
            if ekEn*1.05 < dkEn
                divergeState = 0;
            else
                ef = df;
            end
        end
        if ekEn > dkEn*19.95 %此时滤波器发散，置0重新收敛
            WFb=zeros(N+1,M); % Block-based FD NLMS
        end

        ekEnV(kk) = ekEn;
        dkEnV(kk) = dkEn;

        hnlLocalMinV(kk) = hnlLocalMin;
        cohxdLocalMinV(kk) = cohxdLocalMin;
        hnlMinV(kk) = hnlMin;
        
        aggrFact = 0.3;
        wCurve = [0; aggrFact*sqrt(linspace(0,1,N))' + 0.1];
        weight = wCurve;

        hnled = weight.*min(hnlPrefAvg, hnled) + (1 - weight).*hnled;

%         od = ovrdSm*fliplr(sqrt(linspace(0,1,N+1))' + 1);
        od = ovrdSm*(sqrt(linspace(0,1,N+1))' + 1);
%         sshift = 1.5*ones(N+1,1);
        sshift = ones(N+1,1);

        hnled = hnled.^(od.*sshift);

        hnl = hnled;
        ef = ef.*(hnl);
%         ef = ef.*(min(1 - cohxd, cohed).^2);
        
        ovrdV(kk) = ovrdSm;
        hnledAvg(kk) = 1-mean(1-cohed(echoBandRange));
        hnlxdAvg(kk) = 1-mean(cohxd(echoBandRange));
        hnlSortQV(kk) = hnlPrefAvgLow;
        hnlPrefAvgV(kk) = hnlPrefAvg;

%         Fmix = ef;
        if (CNon)
            snn = sqrt(Sym);
            snn(1) = 0;% Reject LF noise
            Un = 10.0*snn.*exp(1j*2*pi.*[0;rand(N-1,1);0]);%增加了幅度
            % Weight comfort noise by suppression
            Un = sqrt(1 - hnled.^2).*Un;
            Fmix = ef + Un;
        else
            Fmix = ef;
        end
        % Overlap and add in time domain for smoothness
        tmp = [Fmix ; flipud(conj(Fmix(2:N)))];
        mixw = wins.*real(ifft(tmp));
        mola = mbuf(end-N+1:end) + mixw(1:N);
        mbuf = mixw;
        ercn(pos:pos+N-1) = mola;%%%%%----you can hear the effect by sound(100*ercn,16000),add by Shichaog
    end % NLPon

% Shift old FFTs
    XFm(:,2:end) = XFm(:,1:end-1);
    YFm(:,2:end) = YFm(:,1:end-1);
    xfwm(:,2:end) = xfwm(:,1:end-1);
    dfm(:,2:end) = dfm(:,1:end-1);


end

out = floor(erfb);
aecout = floor(ercn);
figure(9);
    subplot(4,1,1);plot(rrin);
    subplot(4,1,2);plot(ssin(delay+1:end));
    subplot(4,1,3);plot(erfb);
    subplot(4,1,4);plot(ercn);
    
fidout=fopen('E:\Matlab_workspace\myWebRTC-audio-processing\out.pcm','wb');
fwrite(fidout,out,'int16');
fclose(fidout);    
fidaecout=fopen('E:\Matlab_workspace\myWebRTC-audio-processing\aecout.pcm','wb');
fwrite(fidaecout,aecout,'int16');
fclose(fidaecout);
hold off
