# make all source files

all : 
	cd CaMa_v405/src && make
	cd share && make
	cd mksrfdata && make
	cd mkinidata && make
	cd main && make
	cd postprocess && make

clean : 
	cd share && make clean
	cd mksrfdata && make clean
	cd mkinidata && make clean
	cd main && make clean
	cd postprocess && make clean
	cd CaMa_v405/src && make clean
