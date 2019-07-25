clear
tablesDVBT
DVBT_symbol_interleaver

channelSelection = [1:4096 4864+1:8960];
channelMapping = mod((0:6816)-(kmin+kmax)/2,8960)+1;

dataset = 'L:\sdr_data\20180906DVBT_CAR\06-Sep-2018 074935.195 522.000MHz 000.wav';

file1 = util.io.SdrConsoleRiff(dataset);
offset = 254027200;
%offset = 540000000;
%offset = 574222222+254027200;
%offset = 1;
fileReader = file1.getReader(11200,offset);
nSamples = fileReader.nSamples - offset;

of = util.comm.OFDMReceiver('guardIntervalDuration',2240,'symbolDuration',8960,'trackingWindow',3,'samplePoint',0.9);
oc = util.comm.OFDMChannelEstimation('pilotIndices',pilotIndex,'pilotValues',pilotValue);

nSymbols = floor(nSamples/11200)-1;
%nSymbols = 9000;
posA = zeros(nSymbols,1);
posA2 = zeros(nSymbols,1);
noiseSpectralDensity = zeros(nSymbols,1);
syncA = zeros(nSymbols,6817);
%syncB = zeros(nSymbols,6817);
syncC = zeros(nSymbols,6817);
tsA =  zeros(nSymbols,1);
tpsA =  zeros(nSymbols,68);
siA =  zeros(nSymbols,1);
sicA =  zeros(nSymbols,1);
vide = zeros(2268,nSymbols,'uint8');
%channelDataA =  zeros(nSymbols,36288);
tic
constellation =[1;33;41;9;11;43;35;3;17;49;57;25;27;59;51;19;21;53;61;29;31;63;55;23;5;37;45;13;15;47;39;7;4;36;44;12;14;46;38;6;20;52;60;28;30;62;54;22;16;48;56;24;26;58;50;18;0;32;40;8;10;42;34;2];
constellation = [32 33 37 36 52 53 49 48 34 35 39 38 54 55 51 50 42 43 47 46 62 63 59 58 40 41 45 44 60 61 57 56 8 9 13 12 28 29 25 24 10 11 15 14 30 31 27 26 2 3 7 6 22 23 19 18 0 1 5 4 20 21 17 16];
qamd = comm.RectangularQAMDemodulator(64,'SymbolMapping','custom',...
    'CustomSymbolMapping',constellation, ...
    'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');
vitd = comm.gpu.ViterbiDecoder('TracebackDepth',136);
deint = comm.ConvolutionalDeinterleaver('NumRegisters',12,'RegisterLengthStep',17);
N = 255;
K = 239;
S = 188;
gp = rsgenpoly(N,K,[],0);
rsd = comm.RSDecoder(N,K,gp,S);
pnseq = step(comm.PNSequence('Polynomial', '1+z^1+z^15', 'InitialConditions', ([1 0 0 1 0 1 0 1 0 0 0 0 0 0 0]), 'Mask', [0 0 0 0 0 0 0 0 0 0 0 0 0 1 1], 'SamplesPerFrame', 1503*8));
pnseq =[255;bi2de(reshape(pnseq,8,[]).','left-msb')];
pnseq(189:188:end) = 0;
for ind = 1:length(posA)
    dataIn = double(fileReader.step());
    [syncspec, pos, pos2] = of.step(dataIn);
        
    if ~isempty(syncspec) && ind>3
        noiseSpectralDensity(ind-1) = mean(abs(syncspec(3410:5552)).^2);
        posA(ind-1)=pos;
        posA2(ind-1)=pos2;
        ofdmSymbol = syncspec(channelMapping).';
        syncA(ind-1,:) = ofdmSymbol;
        [~, pilotState,syncC(ind-1,:)] = oc.step(ofdmSymbol);
        %[syncB(ind-1,:), pilotState,syncC(ind-1,:)] = oc.step(ofdmSymbol);
%         pilotState = mod(pilotState-2,4)+1;
%         dataChannels = channelIndex{pilotState}+1;
%         channelData = syncA(ind-2,dataChannels)./syncB(ind-1,dataChannels)*sqrt(42);
%         channelTransferFun = syncB(ind-1,dataChannels);
%         if mod(pilotState,2) == 1
%             channelData = channelData(H8k+1);
%             channelTransferFun = channelTransferFun(H8k+1);
%         else
%             channelData(H8k+1) = channelData;
%             channelTransferFun(H8k+1) = channelTransferFun;
%         end
%         qamdd = qamd(channelData.');
%         qamdd(HB8k) = reshape(reshape(qamdd,6,[]).*abs(channelTransferFun).^2,[],1);
% 
%         %channelDataA(ind-1,:) = qamdd;
%         vide(:,ind-1) = uint8(bi2de(reshape(vitd.step(qamdd),8,[]).','left-msb'));
%         %plot(20*log10(abs(ifft(syncspec(channelSelection).*scorrf))));
        
       % plot(20*log10(abs(syncspec(channelSelection))));
       % ylim([40,120])
       % drawnow;
        
    end
    if mod(ind,893)==0
        disp("At " + num2str(ind) + "/" + num2str(nSymbols))
        toc
        tic
    end
end
%clockEstimate = polyfit(0:length(posA)-2,posA(1:end-1).',1);
%clock = polyval(clockEstimate,0:length(posA)-2);
%clockError = posA(1:end-1).'-clock-(0.1*2240); %Correct for samplePoint
%imagesc(angle(syncA(1:6000,:).*pilotValues.'.*exp(-2i*pi/8960*clockError(1:6000).'*((0:6816)-(kmin+kmax)/2))));
%saa = syncA(1:end-1,:).*exp(-2i*pi/8960*clockError.'*((0:6816)-(kmin+kmax)/2));
%clear syncA
%%


syncC(:,(sum(abs(syncC)>0) <= size(syncC,1)/4)) = syncC(:,(sum(abs(syncC)>0) <= size(syncC,1)/4))*4;
kcc = fftshift(fft2(syncC));
%kcc(1:3/8*end,:)=0;
%kcc(5/8*end:end,:)=0;
%kcc(:,2/3*end:end)=0;
%kcc(:,1:1/3*end)=0;

kcc(:,2/3*end:end)=0;
kcc(:,1:1/3*end)=0;
kcc(1:30000,:)=0;
kcc(70000:end,:)=0;
kcc(1:46000,1:3300)=0;
kcc(1:46000,3500:end)=0;
kcc(55000:end,1:3300)=0;
kcc(55000:end,3500:end)=0;

syncC=ifft2(ifftshift(kcc))*3;

%%
fileReader = file1.getReader(11200,offset);
of = util.comm.OFDMReceiver('guardIntervalDuration',2240,'symbolDuration',8960,'trackingWindow',3,'samplePoint',0.9);
oc = util.comm.OFDMChannelEstimation('pilotIndices',pilotIndex,'pilotValues',pilotValue);
for ind = 1:length(posA)
    dataIn = double(fileReader.step());
    [syncspec, pos, ts] = of.step(dataIn);
        
    if ~isempty(syncspec) && ind>3
        noiseSpectralDensity(ind-1) = mean(abs(syncspec(3410:5552)).^2);
        posA(ind-1)=pos;
        ofdmSymbol = syncspec(channelMapping).';
        syncA(ind-1,:) = ofdmSymbol;
        [~, pilotState] = oc.step(ofdmSymbol);
        %pilotState = mod(pilotState-2,4)+1;
        dataChannels = channelIndex{pilotState}+1;
        channelData = syncA(ind-1,dataChannels)./syncC(ind-1,dataChannels)*sqrt(42);
        tpsA(ind-1,:) = syncA(ind-1,tpsCarrier+1)./syncC(ind-1,tpsCarrier+1);
        channelTransferFun = syncC(ind-1,dataChannels);
        if mod(pilotState,2) == 1
            channelData = channelData(H8k+1);
            channelTransferFun = channelTransferFun(H8k+1);
        else
            channelData(H8k+1) = channelData;
            channelTransferFun(H8k+1) = channelTransferFun;
        end
        qamdd = qamd(channelData.');
        qamdd(HB8k) = reshape(reshape(qamdd,6,[]).*abs(channelTransferFun).^2,[],1);

        %channelDataA(ind-1,:) = qamdd;
        vide(:,ind-1) = uint8(bi2de(reshape(vitd.step(qamdd),8,[]).','left-msb'));
        %plot(20*log10(abs(ifft(syncspec(channelSelection).*scorrf))));
        
       % plot(20*log10(abs(syncspec(channelSelection))));
       % ylim([40,120])
       % drawnow;
        
    end
    if mod(ind,893)==0
        disp("At " + num2str(ind) + "/" + num2str(nSymbols))
        toc
        tic
    end
end
%bchd = comm.BCHDecoder(67,53,bchgenpoly(127,113));
%bche = comm.BCHEncoder(67,53,bchgenpoly(127,113));
%%
deint.release
rsd.release
deint.reset
rsd.reset
ofs1  = 1;
ofs2=1;
[rsdata, rserr] = rsd(deint(vide((1:floor((numel(vide)-17-2268*ofs1-204*ofs2)/204)*204)+17+2268*ofs1+ofs2*204).'));
pnlong = repmat(uint8(pnseq),ceil(length(rsdata)/188/8),1);
mpegts = bitxor(rsdata,pnlong(1:length(rsdata)));
imagesc(reshape(mpegts,188,[]))
%%
mpegts(1:188:end) = 71;
mpegts(2:188:end)  = bitset(mpegts(2:188:end),8,rserr == -1);
% fid = fopen('dvbt_temp.ts','w');
% fwrite(fid,mpegts);
% fclose(fid);
save('dvbt.h5','syncA', 'syncC', 'posA', 'posA2','tpsA','noiseSpectralDensity','-v7.3')
%fid = fopen('dvbt_temp.viterb','w');
%fwrite(fid,vide((1:floor((numel(vide)-17-2268*ofs1-204*ofs2)/204)*204)+17+2268*ofs1+ofs2*204))
%fclose(fid);