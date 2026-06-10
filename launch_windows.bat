@echo off
REM Gene Expression Explorer - Windows launcher
REM Double-click to start the Shiny app. Closes when you close the browser tab
REM AND press any key in this window.

setlocal enableextensions

REM Move to the folder containing this script (works from any shortcut location)
cd /d "%~dp0"

REM Find Rscript: try PATH first, then the standard install location
set "RSCRIPT="
where Rscript >nul 2>&1
if %errorlevel%==0 (
    set "RSCRIPT=Rscript"
) else (
    for /f "delims=" %%D in ('dir /b /ad /o-n "C:\Program Files\R" 2^>nul') do (
        if exist "C:\Program Files\R\%%D\bin\Rscript.exe" (
            set "RSCRIPT=C:\Program Files\R\%%D\bin\Rscript.exe"
            goto :found
        )
    )
    for /f "delims=" %%D in ('dir /b /ad /o-n "C:\Program Files\R" 2^>nul') do (
        if exist "C:\Program Files\R\%%D\bin\x64\Rscript.exe" (
            set "RSCRIPT=C:\Program Files\R\%%D\bin\x64\Rscript.exe"
            goto :found
        )
    )
)
:found

if "%RSCRIPT%"=="" (
    echo.
    echo ERROR: Could not find Rscript.
    echo Install R from https://cran.r-project.org/ and try again.
    echo.
    pause
    exit /b 1
)

echo Starting Gene Expression Explorer...
echo Using: %RSCRIPT%
echo.
echo Leave this window open while you use the app.
echo Close the browser tab AND this window when you are done.
echo.

"%RSCRIPT%" run.R

if errorlevel 1 (
    echo.
    echo The app exited with an error. See messages above.
    pause
)

endlocal
