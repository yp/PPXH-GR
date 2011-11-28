.PHONY: all reall main remain build
#####################
# Compiler options for production code
#
OPTP=-O2 -march=native
#####################

#####################
# Compiler options for debugging code
#
OPTD=-g -Wall
#####################

#####################
# Pre-processor symbols
#
DFLAGS=-DNDEBUG
#####################

COMPFLAGS=$(OPTP) $(OPTD)

INCLUDE=-I .

LIBS=-lm  -lgsl -lgslcblas

CC=gcc
CXX=g++

CFLAGS= $(COMPFLAGS) $(DFLAGS) $(INCLUDE)
BIN_DIR= bin
SRC_DIR= src
OBJ_DIR= obj
DOC_DIR= doc
ALL_DIR= $(BIN_DIR) $(SRC_DIR) $(OBJ_DIR)



SOURCE= \
        $(SRC_DIR)/abort.c       \
        $(SRC_DIR)/abort.h       \
        $(SRC_DIR)/bit_matrix.c  \
        $(SRC_DIR)/bit_matrix.h  \
        $(SRC_DIR)/bit_vector.c  \
        $(SRC_DIR)/bit_vector.h  \
        $(SRC_DIR)/my_time.c     \
        $(SRC_DIR)/my_time.h     \
        $(SRC_DIR)/rnd.c         \
        $(SRC_DIR)/rnd.h         \
        $(SRC_DIR)/util.c        \
        $(SRC_DIR)/util.h	 \
	$(SRC_DIR)/heuristic.c

OBJ= \
        $(OBJ_DIR)/abort.o      \
        $(OBJ_DIR)/bit_matrix.o \
        $(OBJ_DIR)/bit_vector.o \
        $(OBJ_DIR)/my_time.o    \
        $(OBJ_DIR)/rnd.o        \
        $(OBJ_DIR)/util.o       \
	$(OBJ_DIR)/heuristic.o


PROG=$(BIN_DIR)/ppxh-gr

all: build 
	@echo 'All compiled!'

.makefile: Makefile

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c .makefile
	@echo 'Compiling $<'; \
	mkdir -pv $(dir $@) ; \
	$(CC) -std=c99 $(CFLAGS) -o $@ -c $<

$(SRC_DIR)/%.c: $(SRC_DIR)/%.h
	@touch .makefile

build: heu
	@echo 'All build!'


heu	:  $(OBJ)
	@echo ' * Linking' $(PROG); \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(PROG) $(CFLAGS) $^ $(LIBS)

clean 	:
	@echo 'Cleaning objects and programs' ; \
	rm -f $(OBJ) $(PROG)

