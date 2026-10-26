@echo off
rem Build the Boost.Math Airy witness with MSVC Release, /fp:precise.
rem Header-only Boost 1.92.0 is expected at third_party\boost_1_92_0.
setlocal
set "ROOT=%~dp0.."
set "VSDIR=C:\Program Files\Microsoft Visual Studio\18\Community"
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (
  echo vcvars64.bat failed
  exit /b 2
)
echo === compiler ===
cl 2>&1 | findstr /C:"Version"
cl /nologo /std:c++17 /EHsc /O2 /W3 /MD /fp:precise ^
   /I"%ROOT%\third_party\boost_1_92_0" ^
   /Fo"%ROOT%\logs\boost_airy_witness.obj" ^
   /Fe"%ROOT%\raw\boost_airy_witness.exe" ^
   "%ROOT%\src\boost_airy_witness.cpp"
exit /b %errorlevel%
