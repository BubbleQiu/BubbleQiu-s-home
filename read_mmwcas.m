function [outData,frame] = read_mmwcas(path,a)
%输入
%路径和帧号

%输出
%每一帧的回波数据和总帧数，outData维度为[samples,chirps,Rx,Tx]

%% 读取数据
%studio配置的参数
samples = 256; %一个chirp256个采样点
numchirpperloop = 12; %一个Loop12个chirp,因为有12个发射天线
numLoops = 32; %一帧32个Loop
numchirpperframe = numchirpperloop * numLoops;
numTx=numchirpperloop;
numRx=4;


Fs =8e6; %采样率
B = 2.56e9;%带宽
K = 80e12;%斜率
Tc= B/K;%chirp的时长

% 每Tc就是一个chirp，Tc*Fs= 256 =samples

c = 3e8;%光速
fc = 77e9;%载波频率
lamda = c/fc;%波长

% %校正参数
fs_calib = 8e6;
Slope_calib = 80e12;

filepath_master= strcat(path,'\master_0000_data.bin');%读取雷达数据，retVal的维度是[rxnums,numChirps*numADCSamples]
filepath_slave1= strcat(path,'\slave1_0000_data.bin');
filepath_slave2= strcat(path,'\slave2_0000_data.bin');
filepath_slave3= strcat(path,'\slave3_0000_data.bin');
Expected_Num_SamplesPerFrame = samples*numchirpperloop*numLoops*numRx*2;
%实部+虚部
fp_master = fopen(filepath_master, 'r');
fp_slave1 = fopen(filepath_slave1, 'r');
fp_slave2 = fopen(filepath_slave2, 'r');
fp_slave3 = fopen(filepath_slave3, 'r');
range = [];
radar_data_Rxchain = [];

%% 获取帧数
frame_master = fread(fp_master,'uint16');
frame = length(frame_master)/Expected_Num_SamplesPerFrame;

%% 向距离维投影
    %?(a-1)*Expected_Num_SamplesPerFrame*2  为什么x2
    fseek(fp_master,(a-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData1 = fread(fp_master,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData1, 16));%最高位是1,逻辑索引
    adcData1(neg)    = adcData1(neg) - 2^16;
    adcData1 = adcData1(1:2:end) + 1j*adcData1(2:2:end);%复数
    adcData1Complex = reshape(adcData1, numRx, samples, numchirpperloop, numLoops);
    adcData1Complex = permute(adcData1Complex, [2 4 1 3]);
    
    fseek(fp_slave1,(a-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData2 = fread(fp_slave1,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData2, 16));
    adcData2(neg)    = adcData2(neg) - 2^16;
    adcData2 = adcData2(1:2:end) + 1j*adcData2(2:2:end);
    adcData2Complex = reshape(adcData2, numRx, samples, numchirpperloop, numLoops);
    adcData2Complex = permute(adcData2Complex, [2 4 1 3]);
    
    fseek(fp_slave2,(a-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData3 = fread(fp_slave2,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData3, 16));
    adcData3(neg)    = adcData3(neg) - 2^16;
    adcData3 = adcData3(1:2:end) + 1j*adcData3(2:2:end);
    adcData3Complex = reshape(adcData3, numRx, samples, numchirpperloop, numLoops);
    adcData3Complex = permute(adcData3Complex, [2 4 1 3]);
    
    fseek(fp_slave3,(a-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData4 = fread(fp_slave3,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData4, 16));
    adcData4(neg)    = adcData4(neg) - 2^16;
    adcData4 = adcData4(1:2:end) + 1j*adcData4(2:2:end);
    adcData4Complex = reshape(adcData4, numRx, samples, numchirpperloop, numLoops);% 4 256 12 32
    adcData4Complex = permute(adcData4Complex, [2 4 1 3]);% 256 32 4 12
    
    radar_data_Rxchain(:,:,1:4,:) = adcData1Complex; %master数据，
    radar_data_Rxchain(:,:,5:8,:) = adcData2Complex; %slave1
    radar_data_Rxchain(:,:,9:12,:) =adcData3Complex; %slave2
    radar_data_Rxchain(:,:,13:16,:) = adcData4Complex;%slave3，全部维度[256 32 16 12]
    
    
%% 取TX1/RX1的数据进行FFT，将这组数据当作reference
load('calibrateResults_high.mat')
%TI 提供 校正
    
    for iTX = 1: numTx
        % 校正矩阵与数据相乘得到校正后的数据
        %use first enabled TX1/RX1 as reference for calibration
        TXind = iTX;
        %construct the frequency compensation matrix
        freq_calib = (calibResult.RangeMat(TXind,:)-calibResult.RangeMat(1,1))*fs_calib/Fs *K/Slope_calib;
        freq_calib = 2*pi*(freq_calib)/(samples * 5);
        correction_vec = (exp(1i*((0:samples-1)'*freq_calib))');
        
        
        freq_correction_mat = repmat(correction_vec, 1, 1,numLoops);
        freq_correction_mat = permute(freq_correction_mat, [2 3 1]);
        outData1TX = radar_data_Rxchain(:,:,:,iTX).*freq_correction_mat;
        
        
        %construct the phase compensation matrix
        phase_calib = calibResult.PeakValMat(1,1)./calibResult.PeakValMat(TXind,:);
        %remove amplitude calibration
        %         if phaseCalibOnly == 1
        phase_calib = phase_calib./abs(phase_calib);
        %         end
        phase_correction_mat = repmat(phase_calib.', 1,samples, numLoops);
        phase_correction_mat = permute(phase_correction_mat, [2 3 1]);
        outData(:,:,:,iTX) = outData1TX.*phase_correction_mat; %每一帧校正后的数据
    end
    fclose all;
end
