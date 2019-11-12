@echo OFF

set proc=6
set fold1=nano\40nm_gasgen\1
set fold2=nano\40nm_gasgen\2
set fold3=nano\40nm_gasgen\3
set fold4=nano\40nm_gasgen\4
set fold5=nano\40nm_gasgen\5
set fold6=nano\40nm_gasgen\6

mpiexec -n %proc% python main.py %fold1%\Start.txt %fold1%
REM python Post.py %fold1%
REM TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold2%\Start.txt %fold2%
REM python Post.py %fold2%
REM TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold3%\Start.txt %fold3%
REM python Post.py %fold3%
REM TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold4%\Start.txt %fold4%
REM python Post.py %fold4%
REM TIMEOUT /T 300

mpiexec -n %proc% python main.py %fold5%\Start.txt %fold5%
REM python Post.py %fold5%

mpiexec -n %proc% python main.py %fold6%\Start.txt %fold6%
REM python Post.py %fold6%