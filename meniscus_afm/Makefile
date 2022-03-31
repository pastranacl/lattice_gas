C := gcc
FLAGS := -std=c99 -lm -O2

BIN_PATH := ./bin
SRC_PATH := ./src
COMP_FNAME := afmeniscus
SRCS := $(wildcard ${SRC_PATH}/*.c)


comp:
		@clear
		@echo -n "Compiling... "
		@${CC} ${SRCS} ${FLAGS} -o ${BIN_PATH}/${COMP_FNAME}
		@echo "Done!"
					
