@echo off
rem gnclib -- статическая библиотека GNC (MSVC VS2026)
if not exist build mkdir build
cl /std:c++17 /O2 /W3 /EHs-c- /GR- /utf-8 /Iinclude /c src\*.cpp /Fobuild\
lib /OUT:build\gnc.lib build\*.obj
