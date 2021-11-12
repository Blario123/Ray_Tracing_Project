Project: Project.cc
	g++ Project.cc -o Project -pthread

TestProject: Project.cc
	g++ Project.cc -DTEST -o TestProject -pthread