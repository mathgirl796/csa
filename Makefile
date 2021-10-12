top=$(CURDIR)
SRC_DIR=$(top)/src
BUILD_DIR=$(SRC_DIR)
src=$(wildcard $(SRC_DIR)/*.c) 
obj=$(patsubst $(SRC_DIR)/%.c, %.o, $(src))
target=./a.out

all:
	cd ./src && $(MAKE) \
	-f $(top)/Makefile a.out CURDIR=$(top) && \
	mv a.out ../a.out  
        
a.out:$(obj)
	gcc $^ -o $@

.c.o:
	gcc -c $<

clean:
	rm -f $(BUILD_DIR)/*.o