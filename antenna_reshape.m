function Out = antenna_reshape(data,lamda,Tc)
%Input:
    %data----12发16收所有的回波数据,维度为[chirps samples Rnum Tnum]，16收已按实际天线布局调整
%Output:
    %Out----维度为[86 1]的虚拟天线的回波数据，用于后续方位角FFT
    chirps = size(data,1);
    samples = size(data,2);
    Rnum = size(data,3);
    Tnum = size(data,4);
    
    %% Doppler 相位补偿
%     v = repmat(  (([1:chirps]-(1+chirps)/2) * (lamda/2/(Tc*chirps*Tnum))).',1,samples  ); %[32,256]
%     compensation = 4*pi*v*Tc/lamda;
%     for i = 1:Tnum
%         for j = 1:Rnum
%             data(:,:,j,i) = squeeze(data(:,:,j,i)).*exp(-1i*compensation*(i-1));
%         end
%     end
%     
    %% 天线重新排列
    load('indU.mat')
    data_reshape = zeros(chirps,samples,144);
    ant = 0;
    for T = 12:-1:4
        for R = 1:Rnum
            ant = ant+1; %16*9 =144
            data_reshape(:,:,ant) = squeeze(data(:,:,R,T));
        end
    end
    
%     virtual_num = size(indU,1);
    Out = data_reshape(:,:,indU);

    
end