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
  double virtual Light_Emitted(Vec3 &position_vector,
                               Vec3 &output_direction) = 0;
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

  Vec3 centre;
  double radius;
  double Light_Emitted(Vec3 &position_vector, Vec3 &output_direction)
  {
    return 0.0;
  }
};

int main()
{
  Vec3 vector(1.0, 1.0, 1.0);
  Sphere sphere1(vector, 1.0);
  std::cout << sphere1.Light_Emitted(vector, vector) << std::endl;
}