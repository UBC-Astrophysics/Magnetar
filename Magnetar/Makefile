C = gcc
AR = ar
CFLAGS = -std=c99 -Wall -Wextra -O3 -fpic
SRC = emissivity.c
OBJECT = $(patsubst %.c, %.o, $(SRC))
SHARED_LIB = $(patsubst %.c, %.so, $(SRC))
STATIC_LIB = $(patsubst %.c, %.a, $(SRC))

.PHONY : clean

all: $(SHARED_LIB) 

$(OBJECT): $(SRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(SHARED_LIB): $(OBJECT)
	$(CC) -shared $< -o $@

$(STATIC_LIB): $(OBJECT)
	$(AR) cr $(STATIC_LIB) $(OBJECT)

clean:
	rm -f $(OBJECT) $(SHARED_LIB) $(STATIC_LIB)
