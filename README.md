

A collection of Xiao-Dong's large scale structure codes

# How to compile

gfortran or ifort. ifort maybe more efficient for some high dimensional array staff

	cd src
	make nompi # compile those do not requiring mpi
	make all # compile all programe

# How to use it

	Add this to ~/.bashrc
		export PATH=/home/yourusername/software/LSSCode/bin:${PATH}
		export LSSPATH=/home/yourusername/software/LSSCode


	Then you can call any program from terminal; type LSS_ and Tab to view them


