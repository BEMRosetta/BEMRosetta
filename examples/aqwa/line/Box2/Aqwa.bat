echo Start: %date% %time% >  time.txt
call "C:\Program Files\ANSYS Inc\v231\aqwa\bin\winx64\Aqwa.exe" Analysis.dat
echo End:   %date% %time% >> time.txt
