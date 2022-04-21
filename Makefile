Project: Project.cc Headers/Ray.h Headers/Image.h Headers/Vec3.h Headers/stb_image_write.h
	g++ -std=c++20 Project.cc -o Project -pthread

FastProject: Project.cc Headers/Ray.h Headers/Image.h Headers/Vec3.h Headers/stb_image_write.h
	g++ -std=c++20 -O3 Project.cc -o FastProject -pthread

ValidateProject: Project.cc Headers/Ray.h Headers/Image.h Headers/Vec3.h Headers/stb_image_write.h
	g++ -std=c++20 Project.cc -DVALIDATE -o ValidateProject -pthread
