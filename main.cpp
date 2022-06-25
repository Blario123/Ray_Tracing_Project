#include <iostream>
#include <memory> // Allow usage of unique pointers
#include <thread> // Allow usage of threads
#include <random> // Allow usage of random number generation
#include <chrono> // Allow usage of timing
#include <fstream> // Allow output to file
#include <algorithm> // Allow usage of max
#include <QApplication>
#include <QStyleFactory>

#include "include/main.h"

// Create a type "Radiance" to be a Vec3 containing RGB components

int main(int argc, char *argv[]) {
	QApplication::setStyle(QStyleFactory::create("fusion"));
	QApplication a(argc, argv);
	GUI g;
	g.show();
	/*
	Observer observer(Vec3(0.0, 0.0, 1.0),
					  Vec3(1.0, 0.0, 0.0),
					  Vec3(0.0, 0.0, 1.0),
					  pi / 2.5,
					  std::vector<unsigned>{1000, 500});
					  
	DOFObserver observer2(Vec3(0.0, 0.0, 0.95),
                        Vec3(1.0, 0.0, 0.0),
                        Vec3(0.0, 0.0, 1.0),
                        pi / 2.5,
                        0.265,
                        0.01,
                        std::vector<unsigned>{1000, 500});
                        
	SceneRender scene(observer);

	double radius = 10000.0;

	Sphere floor(Vec3(1.0, 0.0, 0.9 - radius), radius);
	Sphere ceiling(Vec3(2.5, 0.0, 1.1 + radius), radius);
	Sphere left_wall(Vec3(2.5, radius + 0.1, 2.5), radius);
	Sphere right_wall(Vec3(2.5, -radius - 0.1, 2.5), radius);
	Sphere back_wall(Vec3(radius + 0.5, 0.0, 2.5), radius);
	Sphere sphere_1(Vec3(0.25, 0.05, 0.93), 0.03);
	Sphere sphere_2(Vec3(0.4, -0.035, 0.95), 0.05);

	floor.BRDF_Fct_Pt = white_BRDF;
	ceiling.BRDF_Fct_Pt = white_BRDF;
	left_wall.BRDF_Fct_Pt = red_BRDF;
	right_wall.BRDF_Fct_Pt = blue_BRDF;
	back_wall.BRDF_Fct_Pt = white_BRDF;
	sphere_1.BRDF_Fct_Pt = pink_BRDF;
	sphere_2.BRDF_Fct_Pt = purple_BRDF;

	ceiling.Light_Emitted_Fct_Pt = ceiling_light_emitted;

	scene.Add_Object(std::make_unique<Sphere>(floor));
	scene.Add_Object(std::make_unique<Sphere>(ceiling));
	scene.Add_Object(std::make_unique<Sphere>(left_wall));
	scene.Add_Object(std::make_unique<Sphere>(right_wall));
	scene.Add_Object(std::make_unique<Sphere>(back_wall));
	scene.Add_Object(std::make_unique<Sphere>(sphere_1));
	scene.Add_Object(std::make_unique<Sphere>(sphere_2));

	unsigned nbounces = 10;
	unsigned nsamples = 10;
	double threshold_intersection_distance = 1.0e-8;
	double threshold_no_intersection_distance = 100.0;
	double finite_difference_size = 1.0e-8;
	bool silent = false;

	scene.Render_Image_Multithreaded(nbounces, nsamples, threshold_intersection_distance, threshold_no_intersection_distance, finite_difference_size).Save("Images/somebodyonce.png");
	*/
	return QApplication::exec();
}
