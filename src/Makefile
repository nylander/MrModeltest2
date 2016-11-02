CC = gcc
CFLAGS  = -O3 -Wno-unused-result -Wall 
LIBS = -lm

TARGET = mrmodeltest2

SRC = $(TARGET).c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	$(RM) $(TARGET)
