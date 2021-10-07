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
    return reflectivity_over_pi;
  }

  bool Intersection_Check(const Vec3 &initial_position,
                          const Vec3 &direction_vector,
                          double &distance)
  {
    return false;
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
private:
  // The position of the observer
  Vec3 position;

  // The direction the observer is looking in
  Vec3 direction;

  // The horizontal and vertical viewing angles in radians
  double horizontal_viewing_angle_radians;
  double vertical_viewing_angle_radians;

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
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    Vec3 random_vector;

    while (true)
    {
      random_vector.x = distribution(generator);
      random_vector.y = distribution(generator);
      random_vector.z = distribution(generator);

      double modulus = random_vector.norm2();

      if (0 < modulus <= 1.0)
      {
        random_vector.normalise();
        if (dot(random_vector, normal) < 0.0)
        {
          random_vector = -random_vector;
        }
        break;
      }
    }

    return random_vector;
  }
}; // End of Scene

int main()
{
  Vec3 centre(1.0, 1.0, 1.0);
  Colour reflectivity(1.0, 1.0, 1.0);
  Sphere test_sphere(centre, 1.0, reflectivity);

  SceneRender test_scene;
  Vec3 normal(-1.0, 3.4, 1.5);
  std::cout << test_scene.random_vector_generator(normal) << std::endl;
}