CC = gcc
CFLAGS  = -O3
LIBS = -lm

TARGET = mrmodeltest2

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LIBS)

clean:
	$(RM) $(TARGET)
