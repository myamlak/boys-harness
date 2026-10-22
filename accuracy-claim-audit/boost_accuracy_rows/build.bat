@echo off
rem Build the gate with the local MSVC and the Boost headers in this tree.
setlocal
set VSDIR=C:\Program Files\Microsoft Visual Studio\18\Community
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 2)
if not exist "%~dp0build" mkdir "%~dp0build"
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /O2 /W3 /MD /bigobj ^
   /I"%~dp0boost_1_92_0" ^
   /Fo"%~dp0build\\" /Fe"%~dp0build\%~1.exe" "%~dp0src\%~1.cpp"
exit /b %errorlevel%
