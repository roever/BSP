all: example1 example2 example3

example1: example1.cpp bsptree.hpp
	clang++ -g --std=c++14 example1.cpp -o example1 -W -Wall

example2: example2.cpp bsptree.hpp
	clang++ -g --std=c++14 example2.cpp -o example2 -W -Wall

example3: example3.cpp bsptree.hpp
	clang++ -g --std=c++14 example3.cpp -o example3 -W -Wall

