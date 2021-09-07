@echo off
echo Copying MSGFPlus.jar to C:\DMS_Programs\ and \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\
pause

@echo on
xcopy ..\target\msgfplus.jar C:\DMS_Programs\MSGFPlus\ /Y
xcopy ..\target\msgfplus.jar \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\MSGFPlus\ /Y

@echo off
pause
