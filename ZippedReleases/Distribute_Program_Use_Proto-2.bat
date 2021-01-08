@echo off
echo Copying MSGFPlus.jar to C:\DMS_Programs\ and \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\
pause

@echo on
xcopy \\proto-2\CI_Publish\MSGFPlus\msgfplus.jar C:\DMS_Programs\MSGFDB /Y
xcopy \\proto-2\CI_Publish\MSGFPlus\msgfplus.jar \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\MSGFDB\ /y

@echo off
pause
