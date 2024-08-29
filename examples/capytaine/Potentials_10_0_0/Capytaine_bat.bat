call conda activate eBase
python "Potentials.py"
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error"