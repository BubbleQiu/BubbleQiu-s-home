%% mmwcas-rf-EVM signal process（range,doppler,angle）
clc;clear all;close all;
tic;
%% 读数据
path =  'C:\Users\11738\Desktop\毕设matlab代码\12_18.06'; %更改目录文件

[data_all,frames] = read_mmwcas(path,1);
%一帧的回波数据和总帧数，outData维度为[samples,chirps,Rx,Tx]

dimension = size(data_all);%[256,32,16,12]
samples = dimension(1);
chirps = dimension(2);
Rx_number = dimension(3);
Tx_number = dimension(4);

Fs = 8e6; %采样率 
B = 160.28e6;  
K = 5e12;     

T_chirp= B/K;
c = 3e8;
fc = 77e9;
lamda = c/(fc+B/2);
%

idle = 5e-6;
ramp_time = 40e-6;
Tc = idle+ramp_time;
% 实际时间 = 闲置or抖动时间 + 斜率时间

d_res = c/2/B;  %距离分辨率
v_res = c/2/(fc+B/2)/chirps/Tx_number/Tc;  %速度分辨率

D_max = d_res * samples;  %最大距离
V_max = v_res * chirps/2;   %最大速度
 
%% 回波处理
micro_doppler = zeros(chirps,frames); %向速度维投影
micro_range = zeros(samples,frames); %向距离维投影
micro_doppler_inte = zeros(chirps,frames); %波束赋形后向速度维投影
micro_range_inte = zeros(samples,frames); %波束赋形后向距离维投影

Rx = [1:16]; %接收天线选择
Tx = [1:12]; %发送天线选择
RxForMIMOProcess = [13:16,1:4,9:12,5:8];
% [13:16,1:4,9:12,5:8]
all_data = zeros(chirps,samples,length(Rx),length(Tx));

H = zeros(1,size(frames,2)-1);
V = zeros(1,size(frames,2)-1);
for a = 2:frames%frames
    data_all = read_mmwcas(path,a);
    
    %% 2D-FFT处理
    for j = 1:length(Rx)
        for k = 1:length(Tx)
            
            r = zeros(chirps,samples); %每一帧出一张RDM图，一行一个快时间，一列一个慢时间
            r = data_all(:,:,Rx(j),Tx(k)).';%转置 32*256
            r = r-repmat(mean(r),chirps,1); %去直流分量[1,256]%减去列平均
            
            %range-fft %快时间
%             max_range = zeros(chirps,2);% [ 最大值坐标 , 最大值 ]
            range = zeros(chirps,samples);
            for i=1:chirps
                b = r(i,:);%一行一个快时间
                range_win = hamming(samples).'; %加窗
                din_win = b.* range_win;
                datafft = fft(din_win,samples);
                range(i,:) = datafft;
%                 max_range(i,1) = find(abs(datafft)==max(abs(datafft)));%最大值坐标
%                 max_range(i,2) = max(abs(datafft));%最大值
%调试 距离fft 每一个chirp，256个点
%                 i
%                 plot([1:samples]*d_res,abs(datafft))
%                 hold on;
%                 plot([1:samples], datafft );
            end
%             relation = zeros(chirps,chirps);
%             range = abs(r);
%             for a = 1:chirps
%                     for b=1:chirps
%                             ave_a = mean(range(a,:));
%                             ave_b = mean(range(b,:));
%                             relation(a,b)=(sum((range(a,:)-ave_a).*(range(b,:)-ave_b)))...
%                                 /(sqrt(sum((range(a,:)-ave_a).^2))*sqrt(sum((range(b,:)-ave_b).^2)));
%                      end
%                 end
%                 a
            %doppler-fft %慢时间
            data = zeros(chirps,samples);
            for i = 1:samples
                range_win = hamming(chirps).';
                din_win = (range(:,i).').*range_win;
                data(:,i) = fft(din_win).';
            end
            data = fftshift(data,1);  %上半和下半
            %将速度维的零频分量搬移至频谱中心,前半行速度为负，后半行速度为正
            %[0 2*pi]  -->  [-pi pi]
            all_data(:,:,j,k) = data;
        end
    end
    
    micro_doppler(:,a) = sum(data,2); %仅使用了最后一个循环的数据（未波束赋形）（最后一个tx和rx）
    micro_range(:,a) = sum(data,1).';%range对应samples256个点，所以对列求和
%%
    %波束赋形
    all_data = all_data(:,:,RxForMIMOProcess,:); %天线重新排序
    data_for_beamforming = antenna_reshape(all_data,lamda,Tc); %选取无重复的86根虚拟天线
    data_inte = sum(data_for_beamforming,3)/86; %正前方波束赋形后数据
    micro_doppler_inte(:,a) = sum(data_inte,2);
    micro_range_inte(:,a) = sum(data_inte,1).';

%     figure(); %正前方波束赋形幅度分布图
%     aaa = abs(reshape(data_inte,1,chirps*samples)); %查看幅度分布情况
%     histfit(aaa,50,'rayleigh') %'normal'、'rician''weibull'
%     xlabel('幅值','Fontsize',16);ylabel('数量','Fontsize',16);
%     title(strcat('第',num2str(a),'帧，幅度分布图'),'Fontsize',16);
%     figure(); %Quantile-Quantile Plot,qq图
%     pd = makedist('rayleigh');
%     qqplot(aaa,pd);

%     test_cdf = makedist('rayleigh','B',135.472);
%     [h ,p]= kstest(aaa,'CDF',test_cdf)

%%
%波束赋形 指定角度的幅度分布
    linear_theta_scan = pi/2;
    W = exp(-1i*(1:86).'*pi*cos(linear_theta_scan)).'/86;%[ 1*86 ]    
    M = reshape(repmat(W,chirps*samples,1),chirps,samples,86);
    data_inte = sum(data_for_beamforming.*M ,3);
    
%直方图
    aaa = abs(reshape(data_inte,1,chirps*samples)); %查看幅度分布情况
    figure();
    grid on;hold on;box on;
    histogram(aaa,50,'FaceColor',[0.3010 0.5450 0.9330],'BinWidth',29.2)%50bin的各个bin宽度29.2
    
%瑞利分布
    pd = fitdist(aaa.','rayleigh');
    N = 2;
    x = 0:1:1000;
    M_2 = mean(aaa.^2);
    sigma = sqrt(M_2/2^(N/2)/gamma(1+N/2));
    y = x/sigma/sigma.*exp(-x.^2/2/sigma^2)*chirps*samples*29.2;%29.2对于50bins的宽度
    plot(x,y,'--','LineWidth',3,'Color', [0 0.52 0.10])

%K分布
    M_2 = mean(aaa.^2);
    M_4 = mean(aaa.^4);
    nu = (M_4/2/M_2^2-1)^(-1);
    Ka = 1/2*sqrt(M_2/nu);
    y=(2/Ka/gamma(nu))*(x/2/Ka).^nu.* besselk(nu-1,x/Ka)*chirps*samples*29.2;%29.2对于50bins的宽度
    plot(x,y,'-r','LineWidth',3)

%韦布尔分布
    pd = fitdist(aaa.','weibull');
    M = log(aaa);
    N_Xn = chirps*samples;
%     p = (6/pi^2*N_Xn/(N_Xn-1)*(mean(aaa.^2)-mean(aaa).^2)).^(-1/2);
%     q = exp(mean(aaa)+0.5772/p);
    p= 1.77511;
    q=185.94;
    y = p/q*(x/q).^(p-1).*exp(-(x/q).^p)*chirps*samples*29.2;%29.2对于50bins的宽度
    plot(x,y,'-.g','LineWidth',3,'Color',[0.5 0 0.8])

    legend('实测数据直方图','瑞利分布','K分布','韦布尔分布')
    xlim([0,800])
    set(gcf,'WindowState','maximized','color','white');%白色背景+全屏
    xlabel('幅值','Fontsize',24);ylabel('数量','Fontsize',24);
    set(gcf,'color','white');%白色背景
    set(gca,'FontSize',24);%坐标数值大小
    
%     print(figure(1),'../Ray_K_Weibull_60.png','-r600','-dpng');%-r600:分辨率

%%
    %功率谱
    S_w = 1/(2*pi*Fs/256)*0.5*abs(micro_range_inte(:,a)).^2;%高频地波雷达复杂背景杂波信息提取p43

%功率谱
%     figure();
%     grid on;box on;hold on;
%     plot(2*pi/256*[1:256],S_w,'LineWidth',2)
%     xlim([0,2*pi])
%     xticklabels({'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi',})
%     set(gcf,'WindowState','maximized','color','white');%白色背景+全屏
%     xlabel('数字角频率{\itw}','Fontsize',24);ylabel('功率谱S({\itw})','Fontsize',24);
%     set(gcf,'color','white');%白色背景
%     set(gca,'FontSize',24);%坐标数值大小
%     hold off;

    m_0 = sum(S_w);
    m_0 = var(S_w);
    H_1_3 = 4* sqrt(m_0);%有效浪高
    H(1,a-1) = H_1_3 ; %每一帧的有效浪高 进行保存

%样本一阶矩 近似= 实测数据的样本一阶矩 估计风速
    mu_1 = mean(S_w); %样本一阶矩 
    alpha = 8.1*10^3;
    beta = 0.74;
    g = 9.8;
    %矩估计
    syms w
    delt=zeros(1,41);
    i=1;
%数值求解 一阶矩不完全相等 做差求最小
    for u = 0.8:0.001:1.2
        fun =  @(w) (alpha*g^2/w^5)*exp(-beta*(g/(u*w))^4);
        A_1 = int(fun,w,0,2*pi)/2/pi;
        delt(1,i) = A_1-mu_1;
        i=i+1;
    end
    V(1,a-1) = find(abs(delt)== min(abs(delt)))*0.001+0.979;%每一帧的速度 进行保存

   
% %%
%     %range-angle图
%     chirp_num = 20; %选取0.2752m/s速度分量处的数据
%     %v=c*chirp_num/2/(fc+B/2)/Tc/chirps
%     linear_theta_scan = 0:(pi/180):pi;
%     W = exp(-1i*(0:85).'*pi*cos(linear_theta_scan))/86;%[ 86 * 181 ]
%     y_Output = abs(W'* squeeze(data_for_beamforming(chirp_num,:,:)).');
%     %[181 * 86] [  86 * 256 ] 181种theta，86天线的阵列响应对每一点samples求和
%  
% %     检验接收端波束赋形
% %     linear_theta_scan = pi/6;
% %     W = exp(-1i*(0:85).'*pi*cos(linear_theta_scan))/86; 
% %     i = 0:(pi/180):pi
% %     x = exp( 1i*(0:85).'*pi*cos(i));
% %     Y = x.*W;
% %     Y = abs(sum(Y,1));
% %     plot(Y)
%       
%     sine_theta = sin(linear_theta_scan);
%     cos_theta = cos(linear_theta_scan);
%     indices_1D = (1:256);
%     [R_mat, sine_theta_mat] = meshgrid(indices_1D*c/2/B,sine_theta);
%     [~, cos_theta_mat] = meshgrid(indices_1D,cos_theta);
%     x_axis = R_mat.*cos_theta_mat; %R_mat :[256 181].*[1 181]
%     y_axis = R_mat.*sine_theta_mat;
% %     figure();
% %     surf(y_axis, x_axis, y_Output,'EdgeColor','none');
%     xlabel('距离 (m)','Fontsize',16);ylabel('距离 (m)','Fontsize',16);
%     title(strcat('第',num2str(a),'帧，range-angle图'),'Fontsize',16);
%     
% % %%    
%     %doppler-angle图
%     samples_num = 128; 
%     linear_theta_scan = 0:(pi/180):pi;
%     W = exp(-1i*(0:85).'*pi*cos(linear_theta_scan))/86;%[ 86 * 181 ]
%     y_Output = abs(W'* squeeze(data_for_beamforming(:,samples_num,:)).');
% %     surf(y_Output)
% %     sine_theta = sin(linear_theta_scan);
% %     cos_theta = cos(linear_theta_scan);
% %     indices_1D = (1-16.5:32-16.5);
% %     [V_mat, sine_theta_mat] = meshgrid(indices_1D*v_res,sine_theta);
% %     [~, cos_theta_mat] = meshgrid(indices_1D,cos_theta);
% %     x_axis = V_mat.*cos_theta_mat; 
% %     y_axis = V_mat.*sine_theta_mat;
% % 	figure();
% %     surf(y_axis, x_axis, y_Output,'EdgeColor','none');
%     x_label = (1-16.5:30/5:32-16.5) *v_res;
%     x_label = round(x_label*100)/100; 
%     set(gca,'xticklabel',x_label);
% %     y_label = [(32-16.5:-32/32:1-16.5),(1-16.5:32/32:32-16.5)] *v_res;
% %     y_label = round(y_label*100)/100; 
% %     set(gca,'yticklabel',y_label);
%     xlabel('速度 (m/s)','Fontsize',16);ylabel('角度 (°)','Fontsize',16);
%     title(strcat('第',num2str(a),'帧，doppler-angle图'),'Fontsize',16);
% %%    
%     %range-doppler图
%     sig_integrate = (sum(sum(abs(all_data),3),4)/192 + 1);
%     %192天线的32x256的数据对位求和
% %     %为啥要+1
% %     figure();
% %     surf(abs(sig_integrate)) %range-doppler图
% %%%%     hold on;
% %%%%     plot(detect(2,:),detect(1,:),'ro')
%     x_label_range = (0:50:300) * (c/2/B);
%     x_label_range = round(x_label_range*100)/100;set(gca,'xticklabel',x_label_range);
%     y_label = (0-17.5:5:35-17.5) * 0.1101; %(lamda/2/((45e-6)*32*12));用帧时间？
%     y_label = round(y_label*100)/100;set(gca,'yticklabel',y_label);
%     set(gca,'yticklabel',y_label);
%     xlabel('距离 (m)','Fontsize',16);ylabel('速度（m/s）','Fontsize',16);title(strcat('第',num2str(a),'帧，range-doppler图'),'Fontsize',16);
%         
    close all;
end
close all;
%% 绘图
% %速度维投影坐标
% micro_doppler_filted = medfilt2(abs(micro_doppler),[3,3]);
% figure();
% % colormap(jet);
% imagesc((10*log10(abs(micro_doppler_filted))));colorbar
% y_label = (5-17.5:30/6:30-17.5) *(0.1101);
% y_label = round(y_label*100)/100; 
% set(gca,'yticklabel',y_label);
% xlabel('帧号','Fontsize',16);ylabel('速度 (m/s)','Fontsize',16);title('time-doppler图','Fontsize',16);
% 
% %距离维投影坐标
% figure();
% micro_range_filted = medfilt2(abs(micro_range),[3,3]);
% % colormap(jet);
% imagesc((10*log10(abs(micro_range_filted))));colorbar
% y_label_range = (50:50:250) * c/2/B;
% y_label_range = round(y_label_range*100)/100;%小数点第三位四舍五入
% set(gca,'yticklabel',y_label_range);
% xlabel('帧号','Fontsize',16);ylabel('距离 (m)','Fontsize',16);title('time-range','Fontsize',16);
% 
% %(正前方接收波束赋形)
% %速度维投影坐标
% micro_doppler_inte_filted = medfilt2(abs(micro_doppler_inte),[3,3]);
% figure();
% % colormap(jet);
% imagesc((10*log10(abs(micro_doppler_inte))));colorbar
% y_label = (5-17.5:30/6:30-17.5) *(0.1101);
% y_label = round(y_label*100)/100;set(gca,'yticklabel',y_label);
% set(gca,'yticklabel',y_label);
% xlabel('帧号','Fontsize',16);ylabel('速度 (m/s)','Fontsize',16);title('time-doppler图(波束赋形)','Fontsize',16);

toc;
