%
% demo file on how to use Nlx2MatSpike_v3. for further details,see the original help in Nlx2MatSpike_v3.m
%
%

%Filename = '/home/urut/MPI/data/testData/TT5.ntt';
Filename = '/Users/ueli/data/turtleA1_feb13/SE3_CSC36.nse';

%     1. Timestamps   
%     2. Sc Numbers
%     3. Cell Numbers
%     4. Params
%     5. Data Points
FieldSelection(1) = 1;
FieldSelection(2) = 1;
FieldSelection(3) = 1;
FieldSelection(4) = 1;
FieldSelection(5) = 1;

ExtractHeader = 1;
ExtractMode = 1;

ModeArray=[]; %all.

[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints,header] = Nlx2MatSpike_v3( Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
which Nlx2MatSpike_v3

spikes=squeeze(DataPoints(:,1,1:1000));
figure;
plot(spikes);