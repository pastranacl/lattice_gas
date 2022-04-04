C := gcc
FLAGS := -std=c99 -lm -O3
FLAGS_DEBUG := -std=c99 -g -lm

BIN_PATH := ./bin
SRC_PATH := ./src
COMP_FNAME := icovirus
SRCS := $(wildcard ${SRC_PATH}/*.c)


comp:
		@clear
		@echo -n "Compiling... "
		@${CC} ${SRCS} ${FLAGS} -o ${BIN_PATH}/${COMP_FNAME}
		@echo "Done!"

debug:
	@clear
	@echo -n "Compiling (debug mode)..."
	@${CC} ${SRCS} ${FLAGS_DEBUG} -o ${BIN_PATH}/${COMP_FNAME}
