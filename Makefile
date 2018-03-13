CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

run: main.o
	g++ main.cpp hash.cpp include/MurmurHash3.cpp -o main

mapper: mapper.o
	g++ mapper.cpp hash.cpp include/MurmurHash3.cpp -o mapper

clean:
	rm -f *.o *.out results/*.txt results/*.pdf data/reference.txt data/longreads.txt
