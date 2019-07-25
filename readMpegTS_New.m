clear;
tic
file = 'L:\sdr_data\20180906DVBT_CAR\522_H_8000_cut.ts';
file = 'L:\MATLAB\dvbt_temp_noerr.ts';

fid = fopen(file,'r');
header = swapbytes(fread(fid,Inf,'*uint32',184));
fclose(fid);

%%
tablesDVBT;
DVBT_symbol_interleaver
sync_byte = uint8(bitshift(header,-24));
transport_error_indicator = uint8(bitget(header,24));
payload_unit_start_indicator = uint8(bitget(header,23));
transport_priority = uint8(bitget(header,22));
pid = uint16(bitand(bitshift(header,-8),8191));
transport_scrambling_control = uint8(bitand(bitshift(header,-6),3));
adaption_field_control = uint8(bitand(bitshift(header,-4),3));
continuity_counter = uint8(bitand(header,15));

programAssocitationTableIndex = find(pid == 0 & transport_error_indicator == 0)-1;
timeTableIndex = find(pid == 20 & transport_error_indicator == 0)-1;
networkSyncIndex = find(pid == 21 & transport_error_indicator == 0)-1;

networkInformationTableIndex = find(pid == 16 & transport_error_indicator == 0)-1;

megaFrameIndex = mod((0:length(sync_byte)-1).'-networkSyncIndex(1)-1,networkSyncIndex(2)-networkSyncIndex(1));

npindex = find(pid == hex2dec('1fff'));
ttable = readTDT(readMpegTSPackets(file,timeTableIndex));
netSync = readMIP(readMpegTSPackets(file,networkSyncIndex));
pat = readPAT(readMpegTSPackets(file,programAssocitationTableIndex));
nit = readNIT(readMpegTSPackets(file,networkInformationTableIndex));


nullp = readMpegTSPackets(file,npindex-1);

packetTimestamps = datetime(ttable.utc_time_mjd,'ConvertFrom','mjd')+duration(ttable.utc_time_bcd);
mfTransmitTimestamps = 100e-9*double(cumsum([int64(netSync.synchronization_time_stamp(1));mod(diff(int64(netSync.synchronization_time_stamp)),1e7)]));
mfTSec = interp1(timeTableIndex,seconds(packetTimestamps-packetTimestamps(1)),networkSyncIndex);

currentNetworkDescriptor = nit.network_descriptors{find(nit.table_id==64,1)};
currentCellId = currentNetworkDescriptor([currentNetworkDescriptor.descriptor_tag] == 108).descriptor.cell_id;
currentTpsData = de2bi(netSync.tps_mip(1),32,'left-msb');
tps.Constellation = currentTpsData(1:2);
tps.Hierarchy = currentTpsData(3:5);
tps.CodeRateHP = currentTpsData(6:8);
tps.CodeRateLP = currentTpsData(6:8);
tps.GuardInterval = currentTpsData(9:10);
tps.TransmissionMode = currentTpsData(11:12);
tps.CellID = de2bi(currentCellId,16,'left-msb');
tps.DVBH = currentTpsData(16:17);

tpsData = zeros(53,4);
tpsData(1:16,1:2:end) = repmat([0 0 1 1 0 1 0 1 1 1 1 0 1 1 1 0].',1,2);
tpsData(1:16,2:2:end) = repmat([1 1 0 0 1 0 1 0 0 0 0 1 0 0 0 1].',1,2);
tpsData(17:22,:) = repmat([0 1 1 1 1 1].',1,4);
tpsData(23:24,:) = de2bi((0:3),2,'left-msb').';
tpsData(25:26,:) = repmat(tps.Constellation.',1,4);
tpsData(27:29,:) = repmat(tps.Hierarchy.',1,4);
tpsData(30:32,:) = repmat(tps.CodeRateHP.',1,4);
tpsData(33:35,:) = repmat(tps.CodeRateLP.',1,4);
tpsData(36:37,:) = repmat(tps.GuardInterval.',1,4);
tpsData(38:39,:) = repmat(tps.TransmissionMode.',1,4);
tpsData(40:47,1:2:end) = repmat(tps.CellID(1:8).',1,2);
tpsData(40:47,2:2:end) = repmat(tps.CellID(9:16).',1,2);
bch = comm.BCHEncoder(67,53);
tpsData = reshape(bch.step(tpsData(:)),67,4);

tpsSignal = zeros(68,68,4);
tpsSignal(:,1,:) = repmat(tpsInitValue.',1,1,4);
for ind = 2:68
tpsSignal(:,ind,:) = reshape((-1).^tpsData(ind-1,:),1,1,4).*tpsSignal(:,ind-1,:);
end
tpsSignal = reshape(tpsSignal,68,68*4);

pnseq = step(comm.PNSequence('Polynomial', '1+z^1+z^15', 'InitialConditions', ([1 0 0 1 0 1 0 1 0 0 0 0 0 0 0]), 'Mask', [0 0 0 0 0 0 0 0 0 0 0 0 0 1 1], 'SamplesPerFrame', 1503*8));
pnseq =uint8([255;bi2de(reshape(pnseq,8,[]).','left-msb')]);
pnseq(189:188:end) = 0;
N = 255;
K = 239;
S = 188;
gp = rsgenpoly(N,K,[],0);
rse = comm.RSEncoder(N,K,gp,S);
convint = comm.ConvolutionalInterleaver('NumRegisters',12,'RegisterLengthStep',17);




fid = fopen(file,'r');
tsdata = fread(fid,[188,inf],'*uint8');
fclose(fid);

for ind=1:8
    pncut = pnseq((1:188)+188*mod(megaFrameIndex(ind),8));
    tsdata(:,ind:8:end) = bitxor(tsdata(:,ind:8:end),pncut);
end

tsdatars = rse.step(tsdata(:));
tsdatars = convint.step(tsdatars);
tsdatarsbin = reshape(de2bi(tsdatars(:),'left-msb',8).',[],1);
%%
convcode = comm.ConvolutionalEncoder();
constellation = [32 33 37 36 52 53 49 48 34 35 39 38 54 55 51 50 42 43 47 46 62 63 59 58 40 41 45 44 60 61 57 56 8 9 13 12 28 29 25 24 10 11 15 14 30 31 27 26 2 3 7 6 22 23 19 18 0 1 5 4 20 21 17 16];
qamm = comm.RectangularQAMModulator(64,'SymbolMapping','custom',...
    'CustomSymbolMapping',constellation, ...
    'BitInput',true);
mfstart = (find(megaFrameIndex==0,1)-1)*1632+1;
tsdatarsbincut = tsdatarsbin(mfstart:end);
numSymbols = floor(length(tsdatarsbincut)/18144);
syncZ = zeros(numSymbols,6817); 
tpsChannels = tpsCarrier+1;
for ind=1:numSymbols
    convdata = convcode.step(tsdatarsbincut(1+(ind-1)*18144:18144+(ind-1)*18144));
    convdata = convdata(HB8k);
    channelData = qamm.step(convdata)/sqrt(42);
    pilotState = mod(ind-1,4)+1;
    if mod(pilotState,2) == 0
        channelData = channelData(H8k+1);
    else
        channelData(H8k+1) = channelData;
    end
    
    dataChannels = channelIndex{pilotState}+1;
    pilotChannels = pilotIndex{pilotState}+1;
    pilotValues = pilotValue{pilotState};
    syncZ(ind,dataChannels) = channelData;
    syncZ(ind,tpsChannels) = tpsSignal(:,mod(ind-1,272)+1);
    syncZ(ind,pilotChannels) = pilotValues;
    ind
end
toc
load('dvbt.h5','-mat','syncA','posA','posA2')
syncA = syncA(10848:10847+93008,:)./syncZ;
posA = posA(10848:10847+93008);
posA2 = posA2(10848:10847+93008);
freqa=[0:ceil(6817/2)-1,-floor(6817/2):-1]';
syncA = ifftshift(syncA,2);
posAA = (posA-posA2).';


symbolDuration = 10240*7e-6/64;
rangeSampleInterval = 8192/6817*7e-6/64; 
timeToFirstSampleOffset = packetTimestamps(1);
timeToFirstSampleOffsetGPS = round(util.math.unix2gps(seconds(timeToFirstSampleOffset-datetime(1970,01,01))));
timeToFirstPulse = (mod(mfTransmitTimestamps(1)+100e-9*double(netSync.maximum_delay(1)),1));
timeToFirstSample = timeToFirstPulse+(0:size(syncA,1)-1)*symbolDuration;

timeErr = (posAA-(0:length(posAA)-1).'*11200)/(10e6*rangeSampleInterval);
timeErr = timeErr-timeErr(1);
%%

util.io.readRTKlibSolution('L:\sdr_data\20180906DVBT_CAR\FDR\pos\20180906T073850_G0.pos')
pos0 = util.io.readRTKlibSolution('L:\sdr_data\20180906DVBT_CAR\FDR\pos\20180906T073850_G0.pos');
positionTimestamps = (seconds(pos0.time-datetime(1980,01,06))-timeToFirstSampleOffsetGPS);
position = interp1(positionTimestamps,pos0.position.',timeToFirstSample).';


transmitterPosition1 = utm2ecef(32,432492,5637462,354+103+48);
transmitterPosition2 = utm2ecef(32,413038,5667028,669+150+48);
transmitterPosition3 = utm2ecef(32,456494,5672918,793+173+48);

range1 = sqrt(sum((position-transmitterPosition1).^2));
range2 = sqrt(sum((position-transmitterPosition2).^2));
range3 = sqrt(sum((position-transmitterPosition3).^2));
delay1s = range1/Constants.c0/rangeSampleInterval;
delay2s = range2/Constants.c0/rangeSampleInterval;
delay3s = range3/Constants.c0/rangeSampleInterval;
%%


function data=readMpegTSPackets(file,indices)
fid = fopen(file,'r');
data = zeros(length(indices),188,'uint8');
for ind = 1:length(indices)
    fseek(fid,indices(ind)*188,-1);
    data(ind,:) = fread(fid,188,'*uint8');
end
fclose(fid);
end

function packet = readTDT(data)
packet.pointer = data(:,5);
packet.table_id = data(:,6);
packet.section_syntax_indicator = bitget(data(:,7),8);
packet.section_length = bitshift(uint16(bitand(data(:,7),15)),8)+uint16(data(:,8));
packet.utc_time_mjd = bitshift(uint16(data(:,9)),8)+uint16(data(:,10));
packet.utc_time_bcd = [bitshift(data(:,11),-4)*10+bitand(data(:,11),15), ...
    bitshift(data(:,12),-4)*10+bitand(data(:,12),15), ...
    bitshift(data(:,13),-4)*10+bitand(data(:,13),15)];

end

function packet = readMIP(data)
packet.syncronization_id = data(:,5);
packet.section_length = data(:,6);
packet.pointer = bitshift(uint16(data(:,7)),8)+uint16(data(:,8));
packet.periodic_flag = bitget(data(:,9),8);
packet.synchronization_time_stamp = bitshift(uint32(data(:,11)),16)+bitshift(uint32(data(:,12)),8)+uint32(data(:,13));
packet.maximum_delay = bitshift(uint32(data(:,14)),16)+bitshift(uint32(data(:,15)),8)+uint32(data(:,16));
packet.tps_mip = bitshift(uint32(data(:,17)),24)+bitshift(uint32(data(:,18)),16)+bitshift(uint32(data(:,19)),8)+uint32(data(:,20));
packet.individual_addressing_length = data(:,21);
if any(packet.individual_addressing_length~=0)
    error('Individual Addressing in MIP not implemented');
end
packet.crc_32 = bitshift(uint32(data(:,22)),24)+bitshift(uint32(data(:,23)),16)+bitshift(uint32(data(:,24)),8)+uint32(data(:,25));
end

function packet = readPAT(data)
packet.pointer = data(:,5);
packet.table_id = data(:,6);
packet.section_syntax_indicator = bitget(data(:,7),8);
packet.section_length = bitshift(uint16(bitand(data(:,7),15)),8)+uint16(data(:,8));
packet.transport_stream_id = bitshift(uint16(data(:,9)),8)+uint16(data(:,10));
packet.version_number = bitand(bitshift(data(:,11),-1),31);
packet.current_next_indicator = bitget(data(:,11),1);
packet.section_number = data(:,12);
packet.last_section_number = data(:,13);
for ind = 14:4:max(packet.section_length)
    packet.program_number(:,(ind-14)/4+1) = bitshift(uint16(bitand(data(:,ind),15)),8)+uint16(data(:,ind+1));
    packet.pid(:,(ind-14)/4+1) = bitshift(uint16(bitand(data(:,ind+2),15)),8)+uint16(data(:,ind+3));
end
%packet.crc_32 = bitshift(uint32(data(:,packet.section_length)),24)+bitshift(uint32(data(:,packet.section_length)),16)+bitshift(uint32(data(:,packet.section_length+1)),8)+uint32(data(:,packet.section_length+4));
end

function packet = readNIT(data)
packet.pointer = data(:,5);
packet.table_id = data(:,6);
packet.section_syntax_indicator = bitget(data(:,7),8);
packet.section_length = bitshift(uint16(bitand(data(:,7),15)),8)+uint16(data(:,8));
packet.network_id = bitshift(uint16(data(:,9)),8)+uint16(data(:,10));
packet.version_number = bitand(bitshift(data(:,11),-1),31);
packet.current_next_indicator = bitget(data(:,11),1);
packet.section_number = data(:,12);
packet.last_section_number = data(:,13);
packet.network_descriptors_length = bitshift(uint16(bitand(data(:,14),15)),8)+uint16(data(:,15));
bytes_nd = packet.network_descriptors_length;
for indp = 1:size(data,1)
    bytes_ndd = bytes_nd(indp);
    dpoint = 16;
    indd = 1;
    clear descr;
    while bytes_ndd>0
        descr(indd).descriptor_tag = data(indp,dpoint);
        dpoint=dpoint+1;
        descr(indd).descriptor_length = data(indp,dpoint);
        dpoint=dpoint+1;
        switch(descr(indd).descriptor_tag)
            case 64
                descr(indd).descriptor = readNetworkNameDescriptor(data(indp,dpoint:dpoint-1+descr(indd).descriptor_length));
            case 108
                descr(indd).descriptor = readCellListDescriptor(data(indp,dpoint:dpoint-1+descr(indd).descriptor_length));
        end
        dpoint=dpoint+double(descr(indd).descriptor_length);
        bytes_ndd = int32(bytes_ndd)-int32(descr(indd).descriptor_length)-2;
        indd = indd+1;
    end
    descriptors{indp,:} = descr;
    packet.transport_stream_loop_length(indp,:) = bitshift(uint16(bitand(data(indp,dpoint),15)),8)+uint16(data(indp,dpoint+1));
    dpoint = dpoint+2;
    bytes_ndd = packet.transport_stream_loop_length(indp,:);
    clear descr;
    while bytes_ndd>0
        descr.transport_stream_id = bitshift(uint16(data(indp,dpoint)),8)+uint16(data(indp,dpoint+1));
        dpoint = dpoint+2;
        descr.original_network_id = bitshift(uint16(data(indp,dpoint)),8)+uint16(data(indp,dpoint+1));
        dpoint = dpoint+2;
        descr.transport_descriptors_length = bitshift(uint16(bitand(data(indp,dpoint),15)),8)+uint16(data(indp,dpoint+1));
        dpoint = dpoint+2;
        bytes_nd2 = descr.transport_descriptors_length;
        clear descr2;
        indd2 = 1;
        while bytes_nd2>0
            descr2(indd2).descriptor_tag = data(indp,dpoint);
            dpoint=dpoint+1;
            descr2(indd2).descriptor_length = data(indp,dpoint);
            dpoint=dpoint+1;
            switch(descr2(indd2).descriptor_tag)
                case 90
                    descr2(indd2).descriptor = readTerrestrialDeliverySystemDescripter(data(indp,dpoint:dpoint-1+descr2(indd2).descriptor_length));
                %case 65
                %    descr2(indd2).descriptor = readNetworkNameDescriptor(data(indp,dpoint:dpoint-1+descr2(indd2).descriptor_length));
                %case 109
                %    descr2(indd2).descriptor = readNetworkNameDescriptor(data(indp,dpoint:dpoint-1+descr2(indd2).descriptor_length));
            end
            dpoint=dpoint+double(descr2(indd2).descriptor_length);
            bytes_nd2 = int32(bytes_nd2)-int32(descr2(indd2).descriptor_length)-2;
            indd2 = indd2+1;
        end
        bytes_ndd = int32(bytes_ndd)-int32(descr.transport_descriptors_length)-6;
        descr.descriptors = descr2;
    end
    descriptors2{indp,:} = descr;
end
packet.network_descriptors = descriptors;
packet.transport_streams = descriptors2;
%packet.crc_32 = bitshift(uint32(data(:,packet.section_length)),24)+bitshift(uint32(data(:,packet.section_length)),16)+bitshift(uint32(data(:,packet.section_length+1)),8)+uint32(data(:,packet.section_length+4));
end

function desc = readNetworkNameDescriptor(data)
desc = char(data);
end

function desc = readCellListDescriptor(data)
desc.cell_id = bitshift(uint16(data(1)),8)+uint16(data(2));
desc.cell_latitude = double(typecast(bitshift(uint16(data(3)),8)+uint16(data(4)),'int16'))/32768*90;
desc.cell_longitude = double(typecast(bitshift(uint16(data(5)),8)+uint16(data(6)),'int16'))/32768*180;
desc.cell_extent_of_latitude = double(bitshift(typecast(bitshift(uint16(data(7)),8)+uint16(data(8)),'int16'),-4))/32768*90;
desc.cell_extent_of_longitude = double(bitshift(typecast(bitshift(uint16(data(8)),12)+bitshift(uint16(data(9)),4),'int16'),-4))/32768*180;
desc.subcell_info_loop_length = data(10);
if desc.subcell_info_loop_length>0
    warning("subcell_info_loop not implemented");
end
end

function desc = readTerrestrialDeliverySystemDescripter(data)
bw = [8e6,7e6,6e6,5e6,NaN];
desc.centre_frequency = 10*double(bitshift(uint32(data(1)),24)+bitshift(uint32(data(2)),16)+bitshift(uint32(data(3)),8)+uint32(data(4)));
desc.bandwidth = bw(bitshift(data(5),-5)+1);
desc.priority = bitget(data(5),5);
desc.Time_Slicing_indicator = bitget(data(5),4);
desc.MPE_FEC_indicator  = bitget(data(5),3);
desc.constellation = bitshift(data(6),-6);
desc.hierarchy_information = bitand(bitshift(data(6),-3),7);
desc.code_rate_HP_stream = bitand(data(6),7);
desc.code_rate_LP_stream  = bitshift(data(7),-5);
desc.guard_interval = bitand(bitshift(data(7),-3),3);
desc.transmission_mode = bitand(bitshift(data(7),-1),3);
desc.other_frequency_flag = bitget(data(7),1);
end

% function desc = readCellListDescriptor(data)
% desc.cell_id = bitshift(uint16(data(1)),8)+uint16(data(2));
% desc.cell_latitude = double(typecast(bitshift(uint16(data(3)),8)+uint16(data(4)),'int16'))/32768*90;
% desc.cell_longitude = double(typecast(bitshift(uint16(data(5)),8)+uint16(data(6)),'int16'))/32768*180;
% desc.cell_extent_of_latitude = double(bitshift(typecast(bitshift(uint16(data(7)),8)+uint16(data(8)),'int16'),-4))/32768*90;
% desc.cell_extent_of_longitude = double(bitshift(typecast(bitshift(uint16(data(8)),12)+bitshift(uint16(data(9)),4),'int16'),-4))/32768*180;
% desc.subcell_info_loop_length = data(10);
% if desc.subcell_info_loop_length>0
%     warning("subcell_info_loop not implemented");
% end
% end