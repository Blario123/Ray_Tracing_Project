Project: Project.cc
	g++ Project.cc -o Project

TestProject: Project.cc
	g++ Project.cc -DTEST -o TestProject