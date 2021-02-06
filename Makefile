# This is a c++ code. Using g++ compiler.
CC=g++

# Compiler flags.
CFLAGS= -g

vpath %.cc /home/wuyi/usr/work/als_both_single/cudd
vpath %.c /home/wuyi/usr/CUDD/cudd-2.5.0/cudd: 


WHERE	= /home/wuyi/usr/CUDD/cudd-2.5.0
#WHERE_ABC	= ../abc
LIBS	= libabc.a $(WHERE)/dddmp/libdddmp.a $(WHERE)/cudd/libcudd.a \
	$(WHERE)/mtr/libmtr.a $(WHERE)/st/libst.a $(WHERE)/util/libutil.a \
	$(WHERE)/epd/libepd.a 
#LIBS	= libabc.a 


MAINSRC = main.cc 
OTHSRC1 = bnet.cc queue.cc write_func.cc read_file.cc helper.cc comp_real_er.cc call_abc.cc
OTHSRC2 = loc_sim_main_v2.cc sim_new_abc_max.cc sim_new_abc_ave.cc basics.cc exdc_helper.cc simu_ckt.cc
OTHSRC3 = cudd_build_v2.cc ntr.cc stack.cc btree.cc red_same.cc simu_ckt_ss.cc sub_sim_ckt.cc
OTHSRC4 = exdc_factor_new_v3_ave.cc exdc_factor_new_v2_max.cc

SRC = $(MAINSRC) $(OTHSRC1)  $(OTHSRC2) $(OTHSRC3) $(OTHSRC4)
OBJ = $(SRC:.cc=.o)
TARGET = main_both_max_weight

# make all runs.
all: $(TARGET)


%.o: %.cc
	$(CC) $(CFLAGS) -o $@  -c $< 

	

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ)  $(LIBS) -lm -ldl -rdynamic -lpthread -lrt 


# make clean
clean:
	rm -f $(OBJ) $(TARGET)

