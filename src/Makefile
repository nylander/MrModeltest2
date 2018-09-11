CFLAGS= -Wall -Wextra -Wpedantic

LDLIBS= -lm

TARGET= mrmodeltest2

.PHONY: all win clean

all: $(TARGET)

win: CFLAGS += -DWIN=1

win: all

clean:
	$(RM) $(TARGET)
