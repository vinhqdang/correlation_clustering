@echo off
rem Runs all local jobs (KaPoCE root bounds, RAMA on the GPU, PACE bound jobs)
rem from Windows CMD.  The work is done inside WSL (Ubuntu) by
rem experiments/run_all_local.sh, in this checkout of the repository.
rem
rem   experiments\run_all_local.bat                  (conda environment "cc")
rem   experiments\run_all_local.bat PACE_JOBS=12     (settings as KEY=VALUE)
rem   set CONDA_ENV=myenv ^& experiments\run_all_local.bat
rem
rem Needs WSL with Ubuntu, and inside it: conda with an environment that has
rem Python >= 3.10, cmake, g++, make and the CUDA toolkit (nvcc).
setlocal
set "ENV=%CONDA_ENV%"
if "%ENV%"=="" set "ENV=cc"

where wsl >nul 2>nul
if errorlevel 1 (
    echo WSL not found. In an administrator terminal run:  wsl --install -d Ubuntu
    echo then reboot and set up conda inside Ubuntu.
    goto :fail
)

for /f "usebackq delims=" %%p in (`wsl wslpath -a "%~dp0.."`) do set "REPO=%%p"
if "%REPO%"=="" (
    echo could not translate "%~dp0.." into a WSL path
    goto :fail
)
echo repository in WSL: %REPO%
echo conda environment: %ENV%
echo settings: %*

rem bash -i reads ~/.bashrc, where "conda init" sets up conda
wsl bash -ic "cd '%REPO%' && find experiments -name '*.sh' -exec sed -i 's/\r$//' {} + && { conda activate %ENV% || { echo 'conda activate %ENV% failed: create it with conda create -y -n %ENV% python=3.11'; exit 1; }; } && env %* bash experiments/run_all_local.sh"
if errorlevel 1 goto :fail
echo.
echo finished
pause
exit /b 0

:fail
echo.
echo FAILED, see the messages above (logs in logs_local\ of the repository)
pause
exit /b 1
