roomDimensions = [10 8 4];%定义房间尺寸，以米为单位（分别为宽度、长度和高度）。
sourceCoord = [2 2 2];    %将发射器视为 Room 空间内的一个点。假设接收器是半径为 8.75 cm 的球体。
receiverCoord = [5 5 1.8];
r = 0.0875;
h = figure;
plotRoom(roomDimensions,receiverCoord,sourceCoord,h);%绘制房间空间以及接收器（红色圆圈）和发射器（蓝色 x）
N = 5000;                         %设置光线数。
rng(0)
rays = RandSampleSphere(N);       %使用辅助函数 RandSampleSphere 生成光线。rays 是一个 N×3 矩阵。每行光线都包含三维光线矢量方向。
size(rays);
FVect = [125 250 500 1000 2000 4000 8000];%相关的频率

A = [.02, .02, .03, .03, .04, .05, .05;  %频率吸收系数
    .02, .02, .03, .03, .03, .04, .04;
    .02, .02, .03, .03, .04, .07, .07;
    .02, .02, .03, .03, .04, .05, .05;
    .02, .02, .03, .03, .04, .05, .05;
    .02, .02, .03, .03, .04, .05, .05];
R = 1-A;%反射系数
D = [.13, .56, .95, .95, .95, .95, 0.95;%散射系数定义为 1 减去镜面反射声能与总反射声能之间的比率。
    .5, .9, .95, .95, .95, .95, 0.95;
    .5, .9, .95, .95, .95, .95, 0.95;
    .13, .56, .95, .95, .95, .95, 0.95;
    .13, .56, .95, .95, .95, .95, 0.95;
    .13, .56, .95, .95, .95, .95, 0.95];
histTimeStep = 0.0040;%设置直方图时间分辨率（以秒为单位）
impResTime = 0.5;    %脉冲响应时间
nTBins = round(impResTime/histTimeStep);
nFBins = length(FVect);
TFHist = zeros(nTBins,nFBins);
%% 通过跟踪每个频带上的射线来计算接收能量直方图。当光线照射到表面时，
% 根据漫射降雨模型 [2] 记录在接收器处看到的漫射光线能量。
% 撞击表面时的新光线方向是镜面反射和随机反射的组合。继续跟踪光线，直到其传播时间超过脉冲响应持续时间。
for iBand = 1:nFBins
    % Perform ray tracing independently for each frequency band.
    for iRay = 1:size(rays,1)
        % Select ray direction
        ray = rays(iRay,:);
        % All rays start at the source/transmitter
        ray_xyz = sourceCoord;
        % Set initial ray direction. This direction changes as the ray is
        % reflected off surfaces.
        ray_dxyz = ray;
        % Initialize ray travel time. Ray tracing is terminated when the
        % travel time exceeds the impulse response length.
        ray_time = 0;
        % Initialize the ray energy. Energy decreases when the ray hits 
        % a surface.
        ray_energy = 1/size(rays,1);

        while (ray_time <= impResTime)

            % Determine the surface that the ray encounters
            [surfaceofimpact,displacement] = getImpactWall(ray_xyz,...
                                             ray_dxyz,roomDimensions);
            
            % Determine the distance traveled by the ray
            distance = sqrt(sum(displacement.^2));

            % Determine the coordinates of the impact point
            impactCoord = ray_xyz+displacement;

            % Update ray location/source
            ray_xyz = impactCoord;

            % Update cumulative ray travel time
            c = 343; % speed of light (m/s)
            ray_time = ray_time+distance/c;

            % Apply surface reflection to ray's energy
            % This is the amount of energy that is not lost through
            % absorption.
            ray_energy = ray_energy*R(surfaceofimpact,iBand);

            % Apply diffuse reflection to ray energy
            % This is the fraction of energy used to determine what is
            % detected at the receiver
            rayrecv_energy = ray_energy*D(surfaceofimpact,iBand);

            % Adjust the amount of reflected energy
            ray_energy = ray_energy * (1-D(surfaceofimpact,iBand));

            % Determine impact point-to-receiver direction.
            rayrecvvector = receiverCoord-impactCoord;

            % Determine the ray's time of arrival at receiver.
            distance = sqrt(sum(rayrecvvector.*rayrecvvector));
            recv_timeofarrival = ray_time+distance/c;

            if recv_timeofarrival>impResTime
                break
            end

            % Determine amount of diffuse energy that reaches the receiver.
            % See (5.20) in [2].

            % Compute received energy
            N = getWallNormalVector(surfaceofimpact);
            cosTheta = sum(rayrecvvector.*N)/(sqrt(sum(rayrecvvector.^2)));
            cosAlpha = sqrt(sum(rayrecvvector.^2)-r^2)/sqrt(sum(rayrecvvector.^2));
            E = (1-cosAlpha)*2*cosTheta*rayrecv_energy;

            % Scale the enrgy by the travelled distance
            E = E/(recv_timeofarrival*c)^2;

            % Update energy histogram
            tbin = floor(recv_timeofarrival/histTimeStep + 0.5);
            TFHist(tbin,iBand) = TFHist(tbin,iBand) + E;

            % Compute a new direction for the ray.
            % Pick a random direction that is in the hemisphere of the
            % normal to the impact surface.
            d = rand(1,3);
            d = d/norm(d);
            if sum(d.*N)<0
                d = -d;
            end

            % Derive the specular reflection with respect to the incident
            % wall
            ref = ray_dxyz-2*(sum(ray_dxyz.*N))*N;

            % Combine the specular and random components
            d = d/norm(d);
            ref = ref/norm(ref);
            ray_dxyz = D(surfaceofimpact,iBand)*d+(1-D(surfaceofimpact,iBand))*ref;
            ray_dxyz = ray_dxyz/norm(ray_dxyz);
        end
    end
end

figure;
bar(histTimeStep*(0:size(TFHist,1)-1),TFHist);
grid on;
xlabel("Time (s)");
legend(["125 Hz","250 Hz","500 Hz","1000 Hz","2000 Hz","4000 Hz"]);
fs = 44100;
V = prod(roomDimensions);
t0 = ((2*V*log(2))/(4*pi*c^3))^(1/3); % eq 5.45 in [2]
poissonProcess = [];
timeValues = [];
t = t0;
while (t<impResTime)
    timeValues = [timeValues t]; %#ok
    % Determine polarity.
    if (round(t*fs)-t*fs) < 0 
        poissonProcess = [poissonProcess 1]; %#ok
    else
        poissonProcess = [poissonProcess -1];%#ok
    end
    % Determine the mean event occurrence (eq 5.44 in [2])
    mu = min(1e4,4*pi*c^3*t^2/V); 
    % Determine the interval size (eq. 5.44 in [2])
    deltaTA = (1/mu)*log(1/rand); % eq. 5.43 in [2])
    t = t+deltaTA;
end
randSeq = zeros(ceil(impResTime*fs),1);
for index=1:length(timeValues)
    randSeq(round(timeValues(index)*fs)) = poissonProcess(index);
end
flow =  FVect/2;
fhigh = FVect*2;
NFFT = 8192;
sfft = dsp.STFT(FFTLength=NFFT,FrequencyRange="onesided");
F = sfft.getFrequencyVector(fs);
RCF = zeros(length(FVect),length(F));
for index0 = 1:length(FVect)
    for index=1:length(F)
        f = F(index);
        if f<FVect(index0) && f>=flow(index0)
            RCF(index0,index) = .5*(1+cos(2*pi*f/FVect(index0)));
        end
        if f<fhigh(index0) && f>=FVect(index0)
            RCF(index0,index) = .5*(1-cos(2*pi*f/(2*FVect(index0))));
        end
    end
end
figure
semilogx(F,RCF(1,:))
hold on
semilogx(F,RCF(2,:))
semilogx(F,RCF(3,:))
semilogx(F,RCF(4,:))
semilogx(F,RCF(5,:))
semilogx(F,RCF(6,:))
semilogx(F,RCF(6,:))
xlabel("Frequency (Hz)")
ylabel("Response")
grid on;
RCF = fftshift(ifft(RCF,NFFT,2,'symmetric'),2);
y = zeros(length(randSeq),6);
for index=1:length(FVect)
    y(:,index) = conv(randSeq,RCF(index,:),'same');
end
impTimes = (1/fs)*(0:size(y,1)-1);
hisTimes = histTimeStep/2 + histTimeStep*(0:nTBins);
W = zeros(size(impTimes,2),numel(FVect));
EdgeF = getBandedgeFrequencies(FVect,fs);
BW = diff(EdgeF);
for k=1:size(TFHist,1)
    gk0 = floor((k-1)*fs*histTimeStep)+1;
    gk1 = floor(k*fs*histTimeStep);
    yy = y(gk0:gk1,:).^2;
    val = sqrt(TFHist(k,:)./sum(yy,1)).*sqrt(BW/(fs/2));
    for iRay=gk0:gk1
        W(iRay,:)= val;
    end
end
y = y.*W;
ip = sum(y,2);
ip = ip./max(abs(ip));\
figure
plot((1/fs)*(0:numel(ip)-1),ip);
grid on;
xlabel("Time (s)");
ylabel("Impulse Response");
%将脉冲响应应用于音频信号。
[audioIn,fs] = audioread("FunkyDrums-44p1-stereo-25secs.mp3");
audioIn = audioIn(:,1);
%通过使用脉冲响应进行过滤来模拟接收到的音频
audioOut = filter(ip,1,audioIn);
audioOut = audioOut/max(audioOut);
%请听几秒钟的原始音频。
T = 10;
sound(audioIn(1:T*fs),fs)
pause(T)
sound(audioOut(1:T*fs),fs)%收听几秒钟的接收音频。