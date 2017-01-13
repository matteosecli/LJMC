function [ allData ] = ImportTooBigData( largedatafile, numLines, colToImport )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(largedatafile, 'rt');
%formatString = repmat('%f,', 1, 1200);
%formatString(end) = [];
formatString = repmat('%*f64', 1, colToImport-1);
formatString = [formatString, '%f64%*[^\n]'];
allData = zeros(numLines, 1);
iLine = 1;
LinesPerTime = 100000;
while iLine <= numLines-LinesPerTime
  allData(iLine:iLine+LinesPerTime-1,1) = cell2mat(textscan(fid, formatString, LinesPerTime));
  %allData = [allData; [data{:}]];
  iLine = iLine + LinesPerTime;
  fprintf('Reading lines from %d to %d. Complete: %.2f%%.\r',iLine,iLine+LinesPerTime-1,(iLine+LinesPerTime-1)/numLines*100);
end
fclose(fid);

fprintf('Import finished.\n')

end