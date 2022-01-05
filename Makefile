Project: Project.cc
	g++ -std=c++20 Project.cc -o Project -pthread

TestProject: Project.cc
	g++ -std=c++20 Project.cc -DTEST -o TestProject -pthread