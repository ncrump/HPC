How to Profile C/C++/Fortran Code in Unix
-----------------------------------------
Uses standard gcc/gfortran compilers
Uses standard GNU profiling tool gprof
-----------------------------------------

1. Write your code 'mycode.c'

2. Compile code with profile directive
	-----------------------------------------------------------
	gcc -Wall -pg mycode.c -o mycode_gprof
	-----------------------------------------------------------
		-Wall - show all compile warnings (optional)
		-pg   - compile with profiler
		-o    - name output file

3. Output executable file 'mycode_gprof.exe' created

4. Run your output executable file
	-----------------------------------------------------------
	time ./mycode_gprof.exe
	-----------------------------------------------------------
		'time' - unix tool to time code (optional)

5. Output binary file 'gmon.out' created for profiler

6. Run gprof profiler on executable and binary files
	-----------------------------------------------------------
	gprof -p -b mycode_gprof.exe gmon.out > mycode_profile.txt
	-----------------------------------------------------------
		-p - print only flat profile
		-b - suppress verbose
		 > - save profile output to .txt file

7. Now optimize your code!
