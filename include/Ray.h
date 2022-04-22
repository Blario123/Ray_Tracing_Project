#ifndef RAY_H
#define RAY_H

#include <iostream>

#include "Vec3.h"

class Ray
{
public:
  Ray(Vec3 &initial_position_, Vec3 &direction_vector_)
  {
    initial_position = initial_position_;
    direction_vector = direction_vector_;
  }

  // Return the initial position of the ray
  Vec3 Get_Initial_Position() const
  {
    return initial_position;
  }

  // Return the direction vector of the ray
  Vec3 Get_Direction_Vector() const
  {
    return direction_vector;
  }

private:
  Vec3 initial_position;
  Vec3 direction_vector;
};

#endif