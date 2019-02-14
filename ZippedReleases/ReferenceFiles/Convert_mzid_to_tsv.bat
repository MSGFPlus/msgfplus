rem Option 1:
MzidToTsvConverter.exe C:\Work\DatasetName_msgfplus.mzid -mzid:C:\Work\DatasetName_msgfplus.mzid -tsv:C:\Work\DatasetName_msgfplus.tsv -unroll -showDecoy

rem Option 2:
"C:\Program Files\Java\jre8\bin\java.exe" -Xmx2000M -cp C:\MSGFPlus\MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i C:\Work\DatasetName_msgfplus.mzid -o C:\Work\DatasetName_msgfplus.tsv -showQValue 1 -showDecoy 1 -unroll 1

