# make all source files

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

ifneq (${CaMa},\#define)

all : 
	cd share && make
	cd mksrfdata && make
	cd mkinidata && make
	cd main && make
	cd postprocess && make

else

all : 
	cd CaMa_v405/src && make
	cd share && make
	cd mksrfdata && make
	cd mkinidata && make
	cd main && make
	cd postprocess && make

endif

clean : 
	cd share && make clean
	cd mksrfdata && make clean
	cd mkinidata && make clean
	cd main && make clean
	cd postprocess && make clean
	cd CaMa_v405/src && make clean
