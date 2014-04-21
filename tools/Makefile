EXEC = reformat2

all: $(EXEC)

%: %.cpp
	$(CXX) $< -o $@

.PHONY : clean
clean:
	rm -f *.o $(EXEC)

