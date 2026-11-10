@echo off
rem Build the held-out sweep against the exact v5.3.1 tag tree.
rem Run from a Visual Studio developer shell with BOYS_TREE set to that tree.
if not defined BOYS_TREE (
  echo set BOYS_TREE to the v5.3.1 tag checkout 1>&2
  exit /b 2
)
where cl >nul 2>nul
if errorlevel 1 (
  echo cl is not on PATH; run from a Visual Studio developer shell 1>&2
  exit /b 2
)
if not exist build mkdir build
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I "%BOYS_TREE%\include" /Fobuild\sweep_heldout.obj /Febuild\sweep_heldout.exe sweep_heldout.cpp
exit /b %ERRORLEVEL%
