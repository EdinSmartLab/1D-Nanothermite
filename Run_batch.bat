@echo OFF

set proc=6
set fold1=nano\70nm\1a
set fold2=nano\70nm\2a
set fold3=nano\70nm\3a
set fold4=nano\70nm\4a
set fold5=nano\70nm\5a

REM mpiexec -n %proc% python main.py %fold1%\Input_File_nt.txt %fold1%
REM python Post.py %fold1%
REM TIMEOUT /T 300

REM mpiexec -n %proc% python main.py %fold2%\Input_File_nt.txt %fold2%
REM python Post.py %fold2%
REM TIMEOUT /T 300

REM mpiexec -n %proc% python main.py %fold3%\Input_File_nt.txt %fold3%
REM python Post.py %fold3%
REM TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold4%\Input_File_nt.txt %fold4%
python Post.py %fold4%
TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold5%\Input_File_nt.txt %fold5%
python Post.py %fold5%