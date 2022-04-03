C := gcc
FLAGS := -std=c99 -lm -O3

BIN_PATH := ./bin
SRC_PATH := ./src
COMP_FNAME := icovirus
SRCS := $(wildcard ${SRC_PATH}/*.c)


comp:
		@clear
		@echo -n "Compiling... "
		@${CC} ${SRCS} ${FLAGS} -o ${BIN_PATH}/${COMP_FNAME}
		@echo "Done!"
					
