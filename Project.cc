// Define basic input/output stream objects
#include <iostream>

// Included for the ability to generate random values from a uniform
// distribution
#include <random>

// Include image write implementation
#include "Image.h"

// Include 3D vectors (Already included in Image.h but kept here for
// readability)
#include "Vec3.h"


// Create a type "Colour" to be a Vec3 containing RGB components
typedef Vec3 Colour;

double pi = M_PI;
double pi_reciprocal = M_1_PI;

// The base class for physical objects present in the scene
class PhysicalObject
{
public:
  // This function returns the components of the light emitted from this
  // object's surface given a position on the surface and the direction of the
  // light.
  virtual Colour Light_Emitted(const Vec3 &position_vector,
                               const Vec3 &output_direction) = 0;

  // This function returns the Bidirection Reflectance Distribution Function
  // given a position, incident light vector and outgoing light vector.
  virtual Colour BRDF(const Vec3 &position,
                      const Vec3 &incident_light_vector,
                      const Vec3 &outgoing_light_vector) = 0;

  // This function returns true if a light ray with a given initial position and
  // direction intersects with this object. The third argument (passed by
  // reference) will return the minimum distance the light ray took to intersect
  // with this object.
  virtual bool Intersection_Check(const Vec3 &initial_position,
                                  const Vec3 &direction_vector,
                                  double &distance) = 0;
}; // End of PhysicalObject

// The derived class of Sphere representing spheres in the scene
class Sphere : public PhysicalObject
{
public:
  // The constructor for a sphere
  Sphere(const Vec3 &centre_,
         const double &radius_,
         const Colour &reflectivity_)
  {
#ifdef TEST
    if (radius_ <= 0.0)
    {
      throw std::invalid_argument(
        "Can not create a sphere with a non-positive radius");
    }

    if (reflectivity_.x < 0.0 || reflectivity_.x > 1.0)
    {
      throw std::invalid_argument(
        "Reflectivity components must be between 0 and 1");
    }

    if (reflectivity_.y < 0.0 || reflectivity_.y > 1.0)
    {
      throw std::invalid_argument(
        "Reflectivity components must be between 0 and 1");
    }

    if (reflectivity_.z < 0.0 || reflectivity_.z > 1.0)
    {
      throw std::invalid_argument(
        "Reflectivity components must be between 0 and 1");
    }
#endif

    // Define the private member data values
    centre = centre_;
    radius = radius_;

    // The reflectivity divided by pi will be used as the BRDF so it is more
    // convenient to simply have that value stored instead of the reflectivity
    reflectivity_over_pi = reflectivity_ * pi_reciprocal;
  }

  Colour Light_Emitted(const Vec3 &position_vector,
                       const Vec3 &output_direction)
  {
    // Temporary
    Vec3 zero_vector;
    return zero_vector;
  }

  Colour BRDF(const Vec3 &position,
              const Vec3 &incident_light_vector,
              const Vec3 &outgoing_light_vector)
  {
    // Return the Lambertian BRDF for now
    return reflectivity_over_pi;
  }

  bool Intersection_Check(const Vec3 &initial_position,
                          const Vec3 &direction_vector,
                          double &distance)
  {
    Vec3 normalised_direction_vector =
      direction_vector / direction_vector.norm();

    // By substituting the equation of a line into the equation of a sphere, we
    // obtain a quadratic equation for the distance a light ray must travel with
    // initial_position initial position and normalised_direction_vector
    // direction to intersect with this sphere.

    // Get the coefficients of this quadratic equation
    double quadratic_coefficient = normalised_direction_vector.norm2();

    double linear_coefficient =
      2.0 * (dot(initial_position, normalised_direction_vector) -
             dot(centre, normalised_direction_vector));

    double constant = initial_position.norm2() + centre.norm2() -
                      2.0 * dot(initial_position, centre) - radius * radius;

    // Find the determinant of this quadratic equation to see how many
    // intersections there are
    double determinant = linear_coefficient * linear_coefficient -
                         4.0 * quadratic_coefficient * constant;

    // If the determinant is negative, there are no intersections
    if (abs(determinant) < 1.0e-8)
    {
      // If the linear coefficient is greater than zero, the distance is
      // negative so there is no intersection
      if (linear_coefficient > 0.0)
      {
        return false;
      }
      else
      {
        // Calculate the distance using the quadratic formula with a zero
        // determinant
        distance = -linear_coefficient / (2.0 * quadratic_coefficient);

        // There is an intersection
        return true;
      }
    }
    else if (determinant < 0.0)
    {
      // There is no intersection
      return false;
    }
    else
    {
      // This Boolean will be used a few times
      const bool linear_coefficient_bigger_than_sqrt_determinant =
        abs(linear_coefficient) > sqrt(determinant);

      if (abs(linear_coefficient) < 1.0e-8)
      {
        // Calculate the distance using the quadratic formula without the linear
        // coefficient
        distance = sqrt(determinant) / (2.0 * quadratic_coefficient);

        // There is an intersection
        return true;
      }
      else if (linear_coefficient > 0.0)
      {
        // If the linear coefficient is negative and has a larger absolute value
        // than the square root of the determinant, the distance will always be
        // negative
        if (linear_coefficient_bigger_than_sqrt_determinant)
        {
          // There is no intersection
          return false;
        }
        else
        {
          // Find the distance using the quadratic formula
          distance = (-linear_coefficient + sqrt(determinant)) /
                     (2.0 * quadratic_coefficient);

          // There is an intersection
          return true;
        }
      }
      else
      // If the linear coefficient is negative
      {
        // If this Boolean is true, there are two positive roots
        if (linear_coefficient_bigger_than_sqrt_determinant)
        {
          // Calculate the smallest positive root of the quadratic equation
          distance = (-linear_coefficient - sqrt(determinant)) /
                     (2.0 * quadratic_coefficient);

          // There is an intersection
          return true;
        }
        else
        {
          // If the linear coefficient is smaller than the square root of the
          // determinant, there is only one positive root
          distance = (-linear_coefficient + sqrt(determinant)) /
                     (2.0 * quadratic_coefficient);

          // There is an intersection
          return true;
        }
      }
    }
  }

private:
  // The centre and radius of the sphere
  Vec3 centre;
  double radius;

  // The reflectivity of the sphere
  Colour reflectivity_over_pi;
}; // End of Sphere


// The class for the observer
class Observer
{
public:
  Observer(const Vec3 &position_,
           const Vec3 &direction_,
           const Vec3 &upward_direction_,
           const double &grid_width_)
  {
    // Define the member variables
    position = position_;
    direction = direction_;
    upward_direction = upward_direction_;
    grid_width = grid_width_;
  }


private:
  // The position of the observer
  Vec3 position;

  // The direction the observer is looking in
  Vec3 direction;

  // The vector defining the upward direction
  Vec3 upward_direction;

  // The size of the width of the "tennis racket"
  double grid_width;

  // The resolution of the observer's camera (e.g 1920x1080 means that there are
  // 1920 columns and 1080 rows of pixels)
  std::vector<unsigned> resolution;
}; // End of Observer


class SceneRender
{
public:
  SceneRender() {}

  // A vector containing pointers to the physical objects in the scene
  std::vector<PhysicalObject *> scene_object_vector_pt;

  // A pointer to the Observer of the scene
  Observer *observer_pt;

  // Random hemisphere vector generator
  Vec3 random_vector_generator(const Vec3 &normal)
  {
    // Create the uniform distribution between -1 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    // Create the Vec3 containing the random Vector in the hemisphere
    Vec3 random_vector;

    // Keep looping until a suitable vector is found
    while (true)
    {
      // Set the components of random_vector to random values inside (-1,1)
      random_vector.x = distribution(generator);
      random_vector.y = distribution(generator);
      random_vector.z = distribution(generator);

      // Find the squared norm
      double modulus = random_vector.norm2();

      // Check the random vector is inside the unit sphere but isn't the zero
      // vector
      if (0 < modulus <= 1.0)
      {
        // Normalise the random vector
        random_vector.normalise();

        // Make sure that the random vector is in the hemisphere defined by the
        // direction of the normal
        if (dot(random_vector, normal) < 0.0)
        {
          random_vector = -random_vector;
        }
        break;
      }
    }

    // Return the random vector
    return random_vector;
  } // End of random_vector_generator

}; // End of Scene

int main()
{
  Vec3 centre(1.0, 1.0, 1.0);
  Colour reflectivity(1.0, 1.0, 1.0);
  Sphere test_sphere(centre, 1.0, reflectivity);

  Vec3 initial_position(0.0, 0.0, 0.0);
  Vec3 direction(1.0, 1.0, 1.0);

  double distance;

  bool test =
    test_sphere.Intersection_Check(initial_position, direction, distance);

  if (test)
  {
    std::cout << distance << std::endl;
  }

  SceneRender test_scene;
  Vec3 normal(-1.0, 3.4, 1.5);
  std::cout << test_scene.random_vector_generator(normal) << std::endl;
}