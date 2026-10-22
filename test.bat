@echo off
rem gnclib -- проверки (MSVC VS2026), как make test; библиотека -- build.bat
rem CSV проверки плескания -> build\slosh_out
if not exist build\gnc.lib call build.bat
if not exist build\tests mkdir build\tests
if not exist build\slosh_out mkdir build\slosh_out
set CLT=cl /nologo /std:c++17 /O2 /W3 /EHs-c- /GR- /utf-8 /Iinclude /Fobuild\tests\
%CLT% tests\test_integrate.cpp /Febuild\tests\test_integrate.exe
build\tests\test_integrate.exe
%CLT% tests\test_ini.cpp build\gnc.lib /Febuild\tests\test_ini.exe
build\tests\test_ini.exe
%CLT% tests\test_point.cpp /Febuild\tests\test_point.exe
build\tests\test_point.exe
%CLT% tests\test_atmosphere.cpp /Febuild\tests\test_atmosphere.exe
build\tests\test_atmosphere.exe
%CLT% tests\test_atm_upper.cpp /Febuild\tests\test_atm_upper.exe
build\tests\test_atm_upper.exe
%CLT% tests\test_common.cpp build\gnc.lib /Febuild\tests\test_common.exe
build\tests\test_common.exe
%CLT% tests\test_slosh.cpp build\gnc.lib /Febuild\tests\test_slosh.exe
build\tests\test_slosh.exe build\slosh_out
