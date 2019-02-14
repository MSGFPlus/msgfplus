@echo off

for /f "tokens=1,2,3,4* delims=/ " %%i in ('date /t') do set DAYOFWEEK=%%i&set MonthNum=%%j&set DayOfMonth=%%k&set YearNum=%%l

Set DateStamp=%YearNum%%MonthNum%%DayOfMonth%
Set ZipFileName=MSGFPlus_v%DateStamp%.zip

Set "WorkingDirectory=%cd%"

Set StagingDirectory=%Temp%\MSGFPlusFilesToZip
if not exist %StagingDirectory%                    mkdir %StagingDirectory%
if not exist %StagingDirectory%\Docs               mkdir %StagingDirectory%\Docs
if not exist %StagingDirectory%\Docs\Examples      mkdir %StagingDirectory%\Docs\Examples
if not exist %StagingDirectory%\MzidToTsvConverter mkdir %StagingDirectory%\MzidToTsvConverter

echo Copying files to %StagingDirectory%

xcopy ..\target\MSGFPlus.jar %StagingDirectory% /d /y
xcopy ..\README.md           %StagingDirectory% /d /y
xcopy ..\LICENSE.txt         %StagingDirectory% /d /y
xcopy ReferenceFiles\*       %StagingDirectory% /d /y

xcopy ..\Docs\* %StagingDirectory%\Docs /d /y
xcopy ..\Docs\examples\* %StagingDirectory%\Docs\Examples /d /y

xcopy C:\DMS_Programs\MzidToTsvConverter\MzidToTsvConverter.exe %StagingDirectory%\MzidToTsvConverter /d /y
xcopy C:\DMS_Programs\MzidToTsvConverter\MzidToTsvConverter.pdb %StagingDirectory%\MzidToTsvConverter /d /y
xcopy C:\DMS_Programs\MzidToTsvConverter\PRISM.dll              %StagingDirectory%\MzidToTsvConverter /d /y
xcopy C:\DMS_Programs\MzidToTsvConverter\PSI_Interface.dll      %StagingDirectory%\MzidToTsvConverter /d /y
xcopy C:\DMS_Programs\MzidToTsvConverter\Readme.md              %StagingDirectory%\MzidToTsvConverter /d /y

echo Creating "%WorkingDirectory%\%ZipFileName%"
pause

pushd %StagingDirectory%
"C:\Program Files\7-Zip\7z.exe" a "%WorkingDirectory%\%ZipFileName%" *
popd

pause
