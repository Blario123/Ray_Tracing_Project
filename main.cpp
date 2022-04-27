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


#ifdef VALIDATE
// Validation case
// Create a light emitted function for the validation case (Only the x_value
// shall be non-zero for convenience)
Radiance validation_light_emitted(const Vec3 &position)
{
  // Return this arbitrarily chosen RGB radiance value
  return Radiance(0.125, 0.0, 0.0);
}

// Create a BRDF for the validation case (Only the x_value shall be non-zero for
// convenience)
Radiance validation_BRDF(const Vec3 &position,
                         const Vec3 &incident_light_vector,
                         const Vec3 &outgoing_light_vector)
{
  // Implement the Lambertian BRDF with a reflectivity of (0.25, 0.5, 0.75) in
  // the RGB components
  return Vec3(0.5 * pi_reciprocal, 0.0, 0.0);
}

int main()
{
  // Create an observer for the validation case. This observer is positioned at
  // the origin, faces towards positive x with upwards orientated with the
  // z-axis. The horizontal field of view is 90 degrees, and the resolution of
  // the image is 1x1.
  Observer validation_observer(Vec3(0.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 0.0),
                               Vec3(0.0, 0.0, 1.0),
                               pi,
                               std::vector<unsigned>(2, 1));

  // Create the scene
  SceneRender validation_scene(validation_observer);

  // Create a single sphere for the validation case
  Sphere validation_sphere(Vec3(0.0, 0.0, 0.0), 1.0);

  // Set the light emitted by this validation sphere
  validation_sphere.Light_Emitted_Fct_Pt = validation_light_emitted;

  // Set the BRDF of the validation sphere to a Lambertian BRDF
  validation_sphere.BRDF_Fct_Pt = validation_BRDF;

  // Add the validation sphere to the validation scene
  validation_scene.Add_Object(std::make_unique<Sphere>(validation_sphere));

  // Storage for variables
  double expected_answer = 0.0;
  double actual_answer = 0.0;
  unsigned maximum_number_of_bounces_validation = 100;
  unsigned number_of_samples = 1000;
  unsigned number_of_Monte_Carlo_samples_no_RR = 1;
  unsigned maximum_number_of_Monte_Carlo_samples = 500;
  double rr_variance = 0.0;
  double current_answer = 0.0;
  unsigned fixed_bounces = 15;
  double threshold_intersection_distance = 1.0e-8;
  double threshold_no_intersection_distance = 1.0e3;
  double finite_difference_size = 1.0e-8;


  // Create an output data file
  std::ofstream validation_output_file;
  validation_output_file.open("Output/Answer_Vs_Expected");

  // Output the headings for each column
  validation_output_file.width(10);
  validation_output_file << "No_Bounces";
  validation_output_file.width(17);
  validation_output_file << "Expected_Answer";
  validation_output_file.width(15);
  validation_output_file << "Actual_Answer" << std::endl;

  // Output data for the validation case with a varying number of bounces
  for (unsigned number_of_bounces = 0;
       number_of_bounces < maximum_number_of_bounces_validation;
       number_of_bounces++)
  {
    // Render the single pixel of the validation case multiple times, and find
    // the mean value of the red component.
    actual_answer = 0.0;
    for (unsigned j = 0; j < number_of_samples; j++)
    {
      actual_answer += validation_scene
                         .Render_Image(number_of_bounces + 1,
                                       number_of_Monte_Carlo_samples_no_RR,
                                       threshold_intersection_distance,
                                       threshold_no_intersection_distance,
                                       finite_difference_size,
                                       true)(0, 0)
                         .x;
    }
    actual_answer /= number_of_samples;

    // Calculate the expected value of the red component of the pixel colour via
    // the power series L = L_e * sigma(k = 0 to number_of_bounces) (pi *
    // brdf)^k.
    expected_answer = 0.0;
    for (unsigned j = 0; j < number_of_bounces + 1; j++)
    {
      expected_answer += pow(pi * validation_BRDF(Vec3(), Vec3(), Vec3()).x, j);
    }
    expected_answer *= validation_light_emitted(Vec3()).x;

    // Output the number of bounces
    validation_output_file.width(10);
    validation_output_file << number_of_bounces;

    // Output the expected answer
    validation_output_file.width(17);
    validation_output_file << expected_answer;

    // Output the actual answer
    validation_output_file.width(15);
    validation_output_file << actual_answer << std::endl;
  }

  // Close the output file
  validation_output_file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Validation case for Variance
  //////////////////////////////////////////////////////////////////////////////

  // Open a new output file
  validation_output_file.open("Output/Variance");

  // Output the headings of the output data file
  validation_output_file.width(10);
  validation_output_file << "No_Samples";
  validation_output_file.width(20);
  validation_output_file << "Expected_Answer";
  validation_output_file.width(25);
  validation_output_file << "Actual_Answer";
  validation_output_file.width(25);
  validation_output_file << "Variance" << std::endl;

  // Calculate the expected value of the red component of the pixel colour via
  // the power series L = L_e * sigma(k = 0 to number_of_bounces) (pi *
  // brdf)^k.
  expected_answer = 0.0;
  for (unsigned j = 0; j < fixed_bounces + 1; j++)
  {
    expected_answer += pow(pi * validation_BRDF(Vec3(), Vec3(), Vec3()).x, j);
  }
  expected_answer *= validation_light_emitted(Vec3()).x;

  // Run the validation case with increasing number of Monte-Carlo samples,
  // increasing in steps of 10
  for (unsigned i = 1; i < maximum_number_of_Monte_Carlo_samples / 10; i++)
  {
    // Reset the mean and variance with each different sample size taken.
    actual_answer = 0.0;
    rr_variance = 0.0;

    // Find the mean and the variance of the rendering algorithm
    for (unsigned j = 0; j < number_of_samples; j++)
    {
      current_answer = validation_scene
                         .Render_Image(fixed_bounces + 1,
                                       i * 10,
                                       threshold_intersection_distance,
                                       threshold_no_intersection_distance,
                                       finite_difference_size,
                                       true)(0, 0)
                         .x;

      actual_answer += current_answer;

      rr_variance += pow((current_answer - expected_answer), 2);
    }
    actual_answer /= number_of_samples;
    rr_variance /= number_of_samples;

    // Output the data for each value of the Monte-Carlo method sample size
    validation_output_file.precision(15);
    validation_output_file.width(10);
    validation_output_file << i * 10;
    validation_output_file.width(20);
    validation_output_file << expected_answer;
    validation_output_file.width(25);
    validation_output_file << actual_answer;
    validation_output_file.width(25);
    validation_output_file << rr_variance << std::endl;
  }

  validation_output_file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Validation case for Russian Roulette
  //////////////////////////////////////////////////////////////////////////////


  // Open a new output file
  validation_output_file.open("Output/RR_Variance");

  // Output the headings of the output data file
  validation_output_file.width(10);
  validation_output_file << "No_Samples";
  validation_output_file.width(20);
  validation_output_file << "Expected_Answer";
  validation_output_file.width(25);
  validation_output_file << "Actual_Answer";
  validation_output_file.width(25);
  validation_output_file << "Variance" << std::endl;

  // Calculate the exact answer for the output case using the formula L =
  // L_e/(1.0 - pi * brdf)
  expected_answer = validation_light_emitted(Vec3()).x /
                    (1.0 - pi * validation_BRDF(Vec3(), Vec3(), Vec3()).x);

  // Run the validation case with increasing number of Monte-Carlo samples,
  // increasing in steps of 10
  for (unsigned i = 1; i < maximum_number_of_Monte_Carlo_samples / 10; i++)
  {
    // Reset the mean and variance with each different sample size taken.
    actual_answer = 0.0;
    rr_variance = 0.0;

    // Find the mean and the variance of the rendering algorithm with Russian
    // Roulette implemented
    for (unsigned j = 0; j < number_of_samples; j++)
    {
      current_answer =
        validation_scene
          .Render_Image_Russian(i * 10,
                                threshold_intersection_distance,
                                threshold_no_intersection_distance,
                                finite_difference_size,
                                true)(0, 0)
          .x;

      actual_answer += current_answer;

      rr_variance += pow((current_answer - expected_answer), 2);
    }
    actual_answer /= number_of_samples;
    rr_variance /= number_of_samples;

    // Output the data for each value of the Monte-Carlo method sample size
    validation_output_file.width(10);
    validation_output_file << i * 10;
    validation_output_file.width(20);
    validation_output_file << expected_answer;
    validation_output_file.width(25);
    validation_output_file.precision(15);
    validation_output_file << actual_answer;
    validation_output_file.width(25);
    validation_output_file << rr_variance << std::endl;
  }

  // Close the output file
  validation_output_file.close();
}
#else

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

    SceneRender scene(observer);

    double radius = 1000.0;

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

    scene
            .Render_Image(nbounces,
                          nsamples,
                          threshold_intersection_distance,
                          threshold_no_intersection_distance,
                          finite_difference_size)
            .Save("Images/Test2.png");
    */
    return QApplication::exec();
}

#endif
