.POSIX:

CC = clang
CFLAGS += -Wall -Wextra -Wshadow -Wundef -O3
LD = clang
LDFLAGS = 

OBJ = build/keygen.o build/shake.o build/rng.o build/fpr.o build/fft.o build/vrfy.o

all: build/libkeygen.a

build:
	-mkdir build

clean:
	-rm -f $(OBJ) build/libkeygen.a

build/libkeygen.a: build $(OBJ)
	ar rcs build/libkeygen.a $(OBJ)

build/fft.o: fft.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/fft.o fft.c

build/fpr.o: fpr.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/fpr.o fpr.c

build/keygen.o: keygen.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/keygen.o keygen.c

build/rng.o: rng.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/rng.o rng.c

build/shake.o: shake.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/shake.o shake.c

build/vrfy.o: vrfy.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o build/vrfy.o vrfy.c
