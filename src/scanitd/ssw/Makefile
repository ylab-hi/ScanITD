CC ?= gcc
CFLAGS = -O3 -fPIC -std=c99 -Wall -Wextra -I. $(ARCHFLAGS)
LDFLAGS = -lm
SRC = ssw.c str_util.c
OBJ = $(SRC:.c=.o)
TARGET = ssw.so

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(ARCHFLAGS) -shared -o $@ $^ $(LDFLAGS)

%.o: %.c ssw.h str_util.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
