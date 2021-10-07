// Define basic input/output stream objects
#include <iostream>

// Include image write implementation
#include "Image.h"

// Include 3D vectors (Already included in Image.h but kept here for
// readability)
#include "Vec3.h"

// The base class for physical objects present in the scene
class PhysicalObject
{
public:
  // This function returns the light emitted from this object's surface given a
  // position on the surface and the direction of the light.
  virtual double Light_Emitted(const Vec3 &position_vector,
                               const Vec3 &output_direction) = 0;

  // This function returns the Bidirection Reflectance Distribution Function
  // given a position, incident light vector and outgoing light vector.
  virtual double BRDF(const Vec3 &position,
                      const Vec3 &incident_light,
                      const Vec3 &outgoing_light) = 0;

  // This function returns true if a light ray with a given initial position and
  // direction intersects with this object. The third argument (passed by
  // reference) will return the minimum distance the light ray took to intersect
  // with this object.
  virtual bool Intersection_Check(const Vec3 &initial_position,
                                  const Vec3 &direction,
                                  double &distance) = 0;
};

// The derived class of Sphere representing spheres in the scene
class Sphere : public PhysicalObject
{
public:
  Sphere(const Vec3 &Centre_Input, const double &Radius_Input)
  {
    centre = Centre_Input;
    radius = Radius_Input;
  }

  double Light_Emitted(const Vec3 &position_vector,
                       const Vec3 &output_direction)
  {
    return 0.0;
  }

  double BRDF(const Vec3 &position,
              const Vec3 &incident_light,
              const Vec3 &outgoing_light)
  {
    return 0.0;
  }

  bool Intersection_Check(const Vec3 &initial_position,
                          const Vec3 &direction,
                          double &distance)
  {
    return false;
  }

private:
  Vec3 centre;
  double radius;
};

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
};

class Scene
{
public:
  // A vector containing pointers to the physical objects in the scene
  std::vector<PhysicalObject *> scene_object_vector_pt;

  // A pointer to the Observer of the scene
  Observer *observer_pt;
};

int main()
{
  Vec3 vector(1.0, 1.0, 1.0);
  Sphere sphere1(vector, 1.0);
  std::cout << sphere1.Light_Emitted(vector, vector) << std::endl;
}