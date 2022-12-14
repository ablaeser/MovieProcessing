% Path to the SBX file to be z-projected
input_sbx = 'D:\2photon\Simone\Simone_Neutrophils\NaVAi6-1\211008\009\NaVAi6-1_211008_009.sbxfix';

% which planes to project?
zProj = {11:15, 19-23};

% max or mean projection? (default = 'mean')
projType = 'mean'; 

% which channels? ('both' will try to project both separatley if 2 channels exist)
projChan = 'both';



[projData, binLims, projTifPath] = WriteSbxZproj(input_sbx, [], 'z',zProj, 'chan',projChan, 'Nscan',100);

%{
, 'dir',runsDir, 'name',tempName, 'sbxType','raw', 'projType',projParam.type, 'monochrome',true,...
    'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite
%}