rem Option 1
"C:\Program Files\Java\jre8\bin\java.exe" -Xmx4000M -jar C:\MSGFPlus\MSGFPlus.jar -s C:\Work\DatasetName.mgf -o C:\Work\DatasetName_msgfplus.mzid -d C:\SequenceDB\MyDatabase.fasta  -t 10ppm -m 0 -inst 1 -e 1 -ti -1,2 -ntt 1 -tda 1 -minLength 6 -maxLength 50 -n 1 -thread 7 -mod MSGFPlus_Mods.txt

rem Option 2
"C:\Program Files\Java\jre8\bin\java.exe" -Xmx4000M -jar C:\MSGFPlus\MSGFPlus.jar -s C:\Work\DatasetName.mgf -o C:\Work\DatasetName_msgfplus.mzid -d C:\SequenceDB\MyDatabase.fasta  -conf Docs\Examples\MSGFPlus_Params.txt
