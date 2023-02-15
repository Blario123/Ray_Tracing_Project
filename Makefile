all: Project FastProject ValidateProject
.PHONY : all

Project: main.cpp include
	g++ -std=c++20 main.cpp -o ./build/Project -pthread

FastProject: main.cpp include
	g++ -std=c++20 -O3 main.cpp -o ./build/FastProject -pthread

ValidateProject: main.cpp include
	g++ -std=c++20 main.cpp -DVALIDATE -o ./build/ValidateProject -pthread

clean:
	rm -f all Project FastProject ValidateProject
.PHONY : clean
