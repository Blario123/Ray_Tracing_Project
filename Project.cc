#include <iostream>
#include <memory> // Allow usage of unique pointers
#include <thread> // Allow usage of threads
#include <random> // Allow usage of random number generation
#include <chrono> // Allow usage of timing
#include <fstream> // Allow output to file
#include <algorithm> // Allow usage of max

#include "Headers/Image.h" // Include image write implementation
#include "Headers/Vec3.h" // Include 3D vectors (Already included in Image.h but kept
// here for readability)
#include "Headers/Ray.h" // Includes the Ray class

// Create a type "Radiance" to be a Vec3 containing RGB components
typedef Vec3 Radiance;

// Define the values of pi and the reciprocal of pi (To make dividing by pi
// faster)
const double pi = M_PI;
const double pi_reciprocal = M_1_PI;

// The base class for SDF-defined surfaces in the scene. With the current
// implementation of SDFs, we may only allow our light rays to stay in the
// positive regions of the SDFs. Otherwise, we may get arbitrary behaviour.
class SDF
{
public:
  // Do nothing destructor
  virtual ~SDF() {}

  // Pure virtual function for the SDF to be defined in derived classes.
  virtual double SDF_Fct(const Vec3 &position) const = 0;

  // The function describing the inverse of the transformation of the position
  // of the current SDF-defined surface. The function is defined via a function
  // pointer.
  virtual Vec3 Inverse_Transformation(const Vec3 &position) const
  {
    // If there is no inverse transformation defined, don't apply the inverse
    // transformation to the position vector.
    if (Inverse_Transformation_Fct_Pt == 0)
    {
      return position;
    }
    else
    {
      // If there is an inverse transformation defined, apply it to the position
      // vector.
      return Inverse_Transformation_Fct_Pt(position);
    }
  }

  // This function returns the components of the light emitted from this
  // object's surface given a position on the surface and the direction of the
  // light.
  virtual Radiance Light_Emitted(const Vec3 &position)
  {
    if (Light_Emitted_Fct_Pt == 0)
    {
      // If the pointer to the light emitted is a null pointer, return the zero
      // vector
      Radiance zero_vector(0.0, 0.0, 0.0);
      return zero_vector;
    }
    else
    {
      // If Light_Emitted_Fct_Pt points to a function, evaluate this function
      return (*Light_Emitted_Fct_Pt)(position);
    }
  } // End of Light_Emitted


  // This function returns the Bidirection Reflectance Distribution Function
  // given a position, incident light vector and outgoing light vector.
  virtual Radiance BRDF(const Vec3 &position,
                        const Vec3 &incident_light_vector,
                        const Vec3 &outgoing_light_vector)
  {
    if (BRDF_Fct_Pt == 0)
    {
      // If the pointer to the BRDF is a null pointer, return the zero vector
      Radiance zero_vector(0.0, 0.0, 0.0);
      return zero_vector;
    }
    else
    {
      // If BRDF_Fct_Pt points to a function, evaluate this function
      return (*BRDF_Fct_Pt)(
        position, incident_light_vector, outgoing_light_vector);
    }
  } // End of BRDF


  // This function obtains the outward unit normal to the current instance of
  // the SDF-defined defined surface via the central finite difference method.
  virtual Vec3 Outward_Normal(const Vec3 &position,
                              const double &finite_difference_size)
  {
    // The central finite difference of the SDF is taken in all three coordinate
    // directions in order to get the gradient of the SDF. This gradient is a
    // normal to any isosurface at 'position'. This gradient is then normalised
    // to get the unit normal, then the dot product is taken with 'direction'.
    // If the dot product is negative, we flip the signs of the normal and take
    // that as the normal.

    // Initialise storage
    Vec3 normal(0.0, 0.0, 0.0);

    // Find the non-unit normal, division by step size is unnecessary since the
    // normal will be normalised
    normal.x =
      SDF_Fct(
        Vec3(position.x + finite_difference_size, position.y, position.z)) -
      SDF_Fct(
        Vec3(position.x - finite_difference_size, position.y, position.z));

    normal.y =
      SDF_Fct(
        Vec3(position.x, position.y + finite_difference_size, position.z)) -
      SDF_Fct(
        Vec3(position.x, position.y - finite_difference_size, position.z));

    normal.z =
      SDF_Fct(
        Vec3(position.x, position.y, position.z + finite_difference_size)) -
      SDF_Fct(
        Vec3(position.x, position.y, position.z - finite_difference_size));

    // Find the unit normal
    normal.normalise();

    return normal;
  }

  // A function pointer to the inverse transformation of the SDF surface
  Vec3 (*Inverse_Transformation_Fct_Pt)(const Vec3 &position) = 0;

  // A function pointer to the light emitted
  Radiance (*Light_Emitted_Fct_Pt)(const Vec3 &position) = 0;

  // A function pointer to the BRDF
  Radiance (*BRDF_Fct_Pt)(const Vec3 &position,
                          const Vec3 &incident_light_vector,
                          const Vec3 &outgoing_light_vector) = 0;

  // In order to keep consistency in the distance of the global SDF, we must
  // enforce that light rays may only travel in the positive region of an SDF.
  // Therefore, we must be able to choose the inside and outside of an SDF which
  // we do by flipping the sign of an SDF defining a shape.
  bool invert_SDF = false;
};

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// End of SDF ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Start of SphereSDF ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// The class defining an sphere via an SDF.
class SphereSDF : public SDF
{
public:
  // SphereSDF constructor: Set the radius and the inside/outside of the sphere
  SphereSDF(const double &radius_, const bool &invert_SDF_)
  {
    // Set the radius
    radius = radius_;

    // The SDF is positive outside of the sphere unless invert_SDF is true.
    invert_SDF = invert_SDF_;
  }

  // Define the SDF function of a sphere.
  double SDF_Fct(const Vec3 &position) const
  {
    // Apply the inverse transformation to the position vector in order to
    // transform the SDF. Calculate the SDF using this inversely transformed
    // position vector.
    double sdf_value = Inverse_Transformation(position).norm() - radius;

    // If we want to invert the SDF, we flip the negative region with the
    // positive region. This is done by flipping the sign of the SDF.
    if (invert_SDF)
    {
      return -sdf_value;
    }
    else
    {
      return sdf_value;
    }
  }

private:
  // Radius of the sphere
  double radius;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// End of SphereSDF ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Start of PhysicalObject ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// The base class for physical objects present in the scene
class PhysicalObject
{
public:
  // Virtual destructor
  virtual ~PhysicalObject() {}

  // This function returns the components of the light emitted from this
  // object's surface given a position on the surface and the direction of the
  // light.
  virtual Radiance Light_Emitted(const Vec3 &position)
  {
    if (Light_Emitted_Fct_Pt == 0)
    {
      // If the pointer to the light emitted is a null pointer, return the zero
      // vector
      Radiance zero_vector(0.0, 0.0, 0.0);
      return zero_vector;
    }
    else
    {
      // If Light_Emitted_Fct_Pt points to a function, evaluate this function
      return (*Light_Emitted_Fct_Pt)(position);
    }
  } // End of Light_Emitted


  // This function returns the Bidirection Reflectance Distribution Function
  // given a position, incident light vector and outgoing light vector.
  virtual Radiance BRDF(const Vec3 &position,
                        const Vec3 &incident_light_vector,
                        const Vec3 &outgoing_light_vector)
  {
    if (BRDF_Fct_Pt == 0)
    {
      // If the pointer to the BRDF is a null pointer, return the zero vector
      Radiance zero_vector(0.0, 0.0, 0.0);
      return zero_vector;
    }
    else
    {
      // If BRDF_Fct_Pt points to a function, evaluate this function
      return (*BRDF_Fct_Pt)(
        position, incident_light_vector, outgoing_light_vector);
    }
  } // End of BRDF

  // This function returns true if a light ray with a given initial position and
  // direction intersects with this object. The third argument (passed by
  // reference) will return the minimum distance the light ray took to intersect
  // with this object. The fourth argument (passed by reference) will return the
  // normal vector on the object at the intersection, in the direction the light
  // ray came from.
  virtual bool Intersection_Check(const Ray &incident_ray,
                                  double &distance) = 0;

  // Finds the normal to the physical object such that the dot product of the
  // normal with the second argument is positive.
  virtual Vec3 Orientated_Normal(const Vec3 &position,
                                 const Vec3 &direction) = 0;

  // A function pointer to the light emitted
  Radiance (*Light_Emitted_Fct_Pt)(const Vec3 &position) = 0;

  // A function pointer to the BRDF
  Radiance (*BRDF_Fct_Pt)(const Vec3 &position,
                          const Vec3 &incident_light_vector,
                          const Vec3 &outgoing_light_vector) = 0;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////// End of PhysicalObject ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Start of Sphere //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// The derived class of Sphere representing spheres in the scene
class Sphere : public PhysicalObject
{
public:
  // The constructor for a sphere requiring only a centre and a radius
  Sphere(const Vec3 &centre_, const double &radius_)
  {
    // Define the private member data values
    centre = centre_;
    radius = radius_;
  }


  // Returns true if a light ray with an initial position and direction vector
  // intersects with this object. Also returns the closest distance from the
  // light ray source to the sphere along with the normal vector to the sphere
  // at the point of intersection.
  bool Intersection_Check(const Ray &light_ray, double &distance)
  {
    Vec3 initial_position = light_ray.Get_Initial_Position();
    Vec3 direction_vector = light_ray.Get_Direction_Vector();

    // By substituting the equation of a line into the equation of a sphere, we
    // obtain a quadratic equation for the distance a light ray must travel with
    // initial_position initial position and direction_vector direction to
    // intersect with this sphere.

    // Get the coefficients of this quadratic equation
    double quadratic_coefficient = direction_vector.norm2();

    double linear_coefficient = 2.0 * (dot(initial_position, direction_vector) -
                                       dot(centre, direction_vector));

    double constant = initial_position.norm2() + centre.norm2() -
                      2.0 * dot(initial_position, centre) - radius * radius;

    // Find the determinant of this quadratic equation to see how many
    // intersections there are
    double determinant = linear_coefficient * linear_coefficient -
                         4.0 * quadratic_coefficient * constant;


    // If the determinant is not positive, there are no intersections.
    if (determinant > 0.0)
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

    // Return false for no intersections if the determinant is not positive
    return false;
  } // End of Intersection_Check


  // Gives the normal to the sphere at a given point such that the dot product
  // with the second argument is positive
  Vec3 Orientated_Normal(const Vec3 &position, const Vec3 &direction)
  {
    // Find the outward normal of a sphere
    Vec3 normal = position - centre;

    // Normalise the vector
    normal.normalise();

    // Ensure that the normal vector given has a positive dot product with the
    // second argument to this function
    if (dot(normal, direction) < 0.0)
    {
      normal = -normal;
    }

    return normal;
  } // End of Orientated_Normal

private:
  // The centre and radius of the sphere
  Vec3 centre;
  double radius;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////// End of Sphere ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Start of Observer /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// The class for the observer
class Observer
{
public:
  Observer(const Vec3 &position_,
           const Vec3 &camera_direction_,
           const Vec3 &upward_direction_,
           const double &horizontal_field_of_view_angle_,
           const std::vector<unsigned> &resolution_)
  {
    // Define the member variables

    // Position of the camera
    position = position_;

    // Direction the camera is pointing in
    unit_camera_direction = camera_direction_ / camera_direction_.norm();

    // The upwards direction of the camera
    unit_upward_direction = upward_direction_ / upward_direction_.norm();

    // The vector pointing "right" from the point of view of the camera
    unit_sideways_direction =
      cross(unit_camera_direction, unit_upward_direction);

    // The horizontal viewing angle(radians)
    horizontal_field_of_view_angle = horizontal_field_of_view_angle_;

    // The resolution of the camera
    resolution = resolution_;

    // The vertical viewing angle
    vertical_field_of_view_angle =
      2.0 * atan((double(resolution[1]) / double(resolution[0])) *
                 tan(horizontal_field_of_view_angle / 2.0));
  }


  // When given a pixel, this function will return the normalised vector from
  // the camera to the given pixel on the "grid" in front of the camera
  Ray Ray_To_Pixel_XY(const unsigned &x_pixel, const unsigned &y_pixel)
  {
    // Add the vector from the camera to the centre of the "tennis racket" with
    // the vector to the pixel specified from the centre of the "tennis racket"

    Vec3 direction_to_pixel = unit_camera_direction;

    // If one of the dimensions of the resolution is 1, the calculation of
    // 'direction_to_pixel' returns nan when we want 0. So we check the
    // dimensions of the resolution in turn.
    if (resolution[0] != 1)
    {
      direction_to_pixel +=
        tan(horizontal_field_of_view_angle / 2.0) *
        ((2.0 * double(x_pixel)) / (double(resolution[0]) - 1.0) - 1.0) *
        unit_sideways_direction;
    }

    if (resolution[1] != 1)
    {
      direction_to_pixel -=
        tan(vertical_field_of_view_angle / 2.0) *
        ((2.0 * double(y_pixel)) / (double(resolution[1]) - 1.0) - 1.0) *
        unit_upward_direction;
    }

    // Normalise this direction vector
    direction_to_pixel.normalise();

    return Ray(position, direction_to_pixel);
  } // End of Ray_To_Pixel_XY

  // Constant access function to the resolution of the camera
  std::vector<unsigned> Get_Resolution() const
  {
    // Just return the vector containing the resolution
    return resolution;
  } // End of Get_Resolution

private:
  // The position of the observer
  Vec3 position;

  // The direction the observer is looking in (unit vector)
  Vec3 unit_camera_direction;

  // The vector defining the upward direction (unit vector)
  Vec3 unit_upward_direction;

  // The vector resulting from the cross product of the upward direction with
  // the direction of the observer
  Vec3 unit_sideways_direction;

  // The horizontal viewing angle
  double horizontal_field_of_view_angle;

  // The vertical viewing angle (calculated from the horizontal field of view
  // and the resolution in order to correspond to square pixels)
  double vertical_field_of_view_angle;

  // The resolution of the observer's camera (e.g 1920x1080 means that there are
  // 1920 columns and 1080 rows of pixels)
  std::vector<unsigned> resolution;
}; // End of Observer


////////////////////////////////////////////////////////////////////////////////
////////////////////////////// End of Observer /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Start of SceneRender //////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// The class for the scene containing the functions needed to render the scene
// as an image
class SceneRender
{
public:
  // Constructor
  SceneRender(Observer &observer_)
  {
    // Only an observer is required in the construction of a scene
    observer_pt = std::make_unique<Observer>(observer_);

    // Get a seed with a random value
    std::random_device random_seed;

    // Create a uniform distribution between -1 and 1
    generator = std::default_random_engine(random_seed());
    distribution = std::uniform_real_distribution<double>(-1.0, 1.0);
  }

  // Add objects to the scene via a unique pointer
  void Add_Object(std::unique_ptr<PhysicalObject> &object_upt)
  {
    object_pt_vector.push_back(std::move(object_upt));
  }


  // The same function as above but now accepting rvalues
  void Add_Object(std::unique_ptr<PhysicalObject> &&object_upt)
  {
    object_pt_vector.push_back(std::move(object_upt));
  }


  // Add SDF defined objects to the scene via a unique pointer
  void Add_Object(std::unique_ptr<SDF> &object_upt)
  {
    sdf_object_pt_vector.push_back(std::move(object_upt));
  }


  // The same function as above but now accepting rvalues
  void Add_Object(std::unique_ptr<SDF> &&object_upt)
  {
    sdf_object_pt_vector.push_back(std::move(object_upt));
  }


  // Random hemisphere vector generator
  Vec3 Hemisphere_Vector_Generator(const Vec3 &normal)
  {
    // The random vector is instantiated first
    Vec3 random_vector;

    // Keep looping over this algorithm until a suitable vector is found
    while (true)
    {
      // The purpose of this algorithm is to generate a normal distribution of
      // points on a unit hemisphere. The hemisphere we have is such that a
      // vector from the centre to a point on the hemisphere has a positive dot
      // product with the argument "normal".
      // The idea of this algorithm is to generate a point from a uniform
      // distribution inside a 2x2x2 cube, if it is outside the unit sphere
      // (diameter 2), discard the point, otherwise map it to the hemisphere by
      // multiplying the position vector of the point with the appropriate
      // value.

      // Set the components of random_vector to random values inside (-1,1)
      random_vector.x = distribution(generator);
      random_vector.y = distribution(generator);
      random_vector.z = distribution(generator);

      // Find the squared norm
      double modulus = random_vector.norm2();

      // Check the random vector is inside the unit sphere but isn't the zero
      // vector
      if (modulus <= 1.0 && modulus > 0.0)
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
  } // End of Hemisphere_Vector_Generator


  // Take the union of the multiple SDFs provided, by taking the minimum value
  double SDF_Fct(const Vec3 &position) const
  {
    // Allocate storage for the minimum SDF value and the value of the current
    // SDF as we loop through them
    double minimum_distance = 0.0;
    double current_distance = 0.0;

    if (sdf_object_pt_vector.size() > 0)
    {
      // Set the minimum SDF value equal to the SDF of the first SDF provided
      minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
    }

    // Loop through the rest of the SDFs and find the minimum value of the SDFs
    // at "position"
    for (unsigned i = 1; i < sdf_object_pt_vector.size(); i++)
    {
      // If we find a smaller SDF value, set minimum_distance to that value
      double current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
      if (current_distance < minimum_distance)
      {
        minimum_distance = current_distance;
      }
    }

    return minimum_distance;
  }

  // Take the union of the multiple SDFs provided, by taking the minimum value.
  // The second argument returns the index of the smallest SDF in
  // sdf_object_pt_vector. The reason for having sdf_intersection_index, is to
  // allow us to find out which surface light rays are intersecting with so we
  // can find the appropriate BRDF and Light_Emitted function.
  double SDF_Fct(const Vec3 &position, int &sdf_intersection_index) const
  {
    // Allocate storage for the minimum SDF value and the value of the current
    // SDF as we loop through them
    double minimum_distance = 0.0;
    double current_distance = 0.0;

    if (sdf_object_pt_vector.size() > 0)
    {
      // Set the minimum SDF value equal to the SDF of the first SDF provided,
      // and
      // set the value of sdf_intersection_index to the index of this first SDF.
      minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
      sdf_intersection_index = 0;
    }

    // Loop through the rest of the SDFs and find the minimum value of the SDFs
    // at "position". Also set sdf_intersection_index to the corresponding index
    // value.
    for (unsigned i = 1; i < sdf_object_pt_vector.size(); i++)
    {
      // If we find a smaller SDF value, set minimum_distance and
      // sdf_intersection_index
      current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
      if (current_distance < minimum_distance)
      {
        minimum_distance = current_distance;
        sdf_intersection_index = i;
      }
    }

    return minimum_distance;
  }

  // Loop over all the objects in the scene, then all the SDF-defined objects in
  // the scene to find the first intersection of light_ray with the objects in
  // the scene. The two int arguments return the index of the closest object
  // that light_ray intersects with, if either return -1, there has been no
  // intersection with the non SDF-defined objects in the scene or the
  // SDF-defined objects in the scene. The argument
  // threshold_intersection_distance is the distance at which we consider a
  // light_ray to have intersected with an SDF-defined surface in the scene. The
  // argument threshold_no_intersection_distance is the distance at which we say
  // there has been no intersection between light_ray and an SDF-defined
  // surface.
  Vec3 First_Intersection_Point(
    const Ray &light_ray,
    int &object_intersection_index,
    int &sdf_intersection_index,
    const double &threshold_intersection_distance,
    const double &threshold_no_intersection_distance)
  {
    // Store the indices of the closest intersection that light_ray intersects
    // with. Only one of these indices can not be -1 when this function returns
    // since only one object is the closest object that light_ray intersects
    // with. (At least for the purposes of coding this up)
    object_intersection_index = -1;
    sdf_intersection_index = -1;

    // Store the SDF value at the origin of light_ray. These two variables are
    // allocated for use in finding the intersection of light_ray with an SDF.
    double safe_travel_distance =
      std::abs(SDF_Fct(light_ray.Get_Initial_Position()));
    double new_safe_travel_distance = safe_travel_distance;

    // Store the smallest distance from the light ray source to an intersection
    // along with the index of the corresponding object in object_pt_vector.
    double smallest_distance = 0.0;

    // Store the distance of the light ray to the intersection with the current
    // object being looped over.
    double current_distance = 0.0;

    Vec3 intersection_point(0.0, 0.0, 0.0);

    // Loop over every object in the scene to find an intersection
    for (unsigned i = 0; i < object_pt_vector.size(); i++)
    {
      // If an intersection has already been found, check whether this new
      // intersection is closer to the light ray source than the previous
      // closest intersection.
      if (object_intersection_index != -1)
      {
        // Check whether this intersection is closer than the previous closest
        // one. If so, replace smallest_distance
        if (object_pt_vector[i]->Intersection_Check(light_ray,
                                                    current_distance) &&
            current_distance < smallest_distance)
        {
          smallest_distance = current_distance;
          object_intersection_index = i;
        }
      }
      // If an intersection hasn't been found yet, any intersection will be the
      // closest intersection so far.
      else
      {
        if (object_pt_vector[i]->Intersection_Check(light_ray,
                                                    current_distance))
        {
          // This section of code will only be invoked when the first
          // intersection is found.
          smallest_distance = current_distance;
          object_intersection_index = i;
        }
      }
    }

    // The closest intersection of light_ray with a non-SDF defined surface has
    // been found above, now we find the closest intersection of light_ray with
    // an SDF defined surface using the standard iteration method.

    // This variable keeps track of the distance travelled along the ray from
    // the origin of light_ray
    double distance_travelled = 0.0;

    // If the SDF is smaller than the distance to the closest intersection with
    // a non-SDF defined surface, or the light ray hasn't intersected with a
    // non-SDF defined surface, we then work on finding the closest intersection
    // with SDF defined surfaces.

    if (sdf_object_pt_vector.size() > 0 &&
        safe_travel_distance < smallest_distance)
    {
      // Travel along the light_ray while the light_ray has not intersected with
      // an SDF, has not gotten too far away from a surface, or has already
      // passed through a non-SDF defined surface
      while (safe_travel_distance >= threshold_intersection_distance &&
             safe_travel_distance <= threshold_no_intersection_distance &&
             (distance_travelled < smallest_distance ||
              object_intersection_index == -1))
      {
        // Move along the ray by "safe_travel_distance" and evaluate the SDF,
        // the value of the SDF is the distance we travel in the next iteration,
        // so we save it as "new_safe_travel_distance"
        new_safe_travel_distance =
          std::abs(SDF_Fct(light_ray.Get_Initial_Position() +
                           (distance_travelled + safe_travel_distance) *
                             light_ray.Get_Direction_Vector()));

        // Add the distance travelled in this iteration to "distance_travelled"
        distance_travelled += safe_travel_distance;

        // Set the distance we travel in the next iteration
        safe_travel_distance = new_safe_travel_distance;
      }

      // If we have intersected with an SDF, we return the index of the SDF that
      // our ray intersected with and the point of intersection.
      if (safe_travel_distance < threshold_intersection_distance)
      {
        // If we get to this block of code, our light ray will have intersected
        // an SDF before a non-SDF surface. Therefore we first find the
        // intersection point.
        intersection_point =
          light_ray.Get_Initial_Position() +
          distance_travelled * light_ray.Get_Direction_Vector();

        // Set the value of sdf_intersection_index to the index of the SDF of
        // intersection.
        SDF_Fct(intersection_point, sdf_intersection_index);

        // Set object_intersection_index to -1 since our point of first
        // intersection is with an SDF.
        object_intersection_index = -1;

        return intersection_point;
      }
    }

    // If we reach this section of the code, our light ray has either
    // intersected with a non-SDF surface first or has gotten too far away from
    // the SDF surfaces.

    // We return the point of intersection between our ray and a
    //  non-SDF surface if there has been an intersection.
    if (object_intersection_index != -1)
    {
      intersection_point =
        (light_ray.Get_Initial_Position() +
         smallest_distance * light_ray.Get_Direction_Vector());

      return intersection_point;
    }

    // If there has been no intersection at all, return the zero vector.
    return Vec3(0.0, 0.0, 0.0);

  } // End of First_Intersection_Point


  Radiance Light_Out(const Ray &light_ray,
                     unsigned bounces_remaining,
                     const double &threshold_intersection_distance,
                     const double &threshold_no_intersection_distance,
                     const double &finite_difference_size)
  {
    // If the light ray can't bounce off a single object, no light will reach
    // the "observer".
    if (bounces_remaining == 0)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // Store the radiance that light_ray will "carry"
    Radiance resulting_light(0.0, 0.0, 0.0);

    // Store indices of the surfaces that light_ray will intersect with
    int object_intersection_index = -1;
    int sdf_intersection_index = -1;

    // Storage for a normal vector
    Vec3 normal(0.0, 0.0, 0.0);

    // Storage for the BRDF
    Radiance BRDF(0.0, 0.0, 0.0);

    // Return the closest intersection of light_ray with a surface, along with
    // the index of the object
    Vec3 intersection_point =
      First_Intersection_Point(light_ray,
                               object_intersection_index,
                               sdf_intersection_index,
                               threshold_intersection_distance,
                               threshold_no_intersection_distance);

    // If no surface was hit by light_ray, we return zero radiance
    if (object_intersection_index == -1 && sdf_intersection_index == -1)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // Past this point, either object_intersection_index or (exclusive or)
    // sdf_intersection_index is -1.

    // If object_intersection_index is -1, we must have intersected with an SDF,
    // therefore we find the Light_Emitted at the point of intersection.
    if (object_intersection_index == -1)
    {
      resulting_light =
        sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(
          intersection_point);

      // Find the outward normal to the surface hit by light_ray at the point of
      // intersection
      normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(
        intersection_point, finite_difference_size);
    }
    else
    {
      // Otherwise we intersected with a non-SDF surface and we find the light
      // emitted by this surface at the point of intersection.
      resulting_light =
        object_pt_vector[object_intersection_index]->Light_Emitted(
          intersection_point);

      // Find the outward normal to the surface hit by light_ray at the point of
      // intersection
      normal = object_pt_vector[object_intersection_index]->Orientated_Normal(
        intersection_point, -light_ray.Get_Direction_Vector());
    }

    // Rupinder: Fix this so that it works again (Currently this can not be
    // copied and pasted since the code has changed)
    /* Uncomment this and put the above in comments in order to only output the
       contribution from a certain bounce.
        if (bounces_remaining == 1)
        {
          // First, find the light emitted by the object in the direction of
          // light_ray
          resulting_light =
       object_pt_vector[object_intersection_index]->Light_Emitted(
            intersection_point);
        }
        */

    // Find the direction of a random incident ray onto an object, along with
    // its point of origin
    Vec3 new_direction = Hemisphere_Vector_Generator(normal);
    Vec3 new_position = intersection_point + 1.0e-8 * normal;
    Ray new_ray(new_position, new_direction);

    // Find the BRDF of the first object hit by light_ray, the code slightly
    // differs for the different types of surfaces.
    if (object_intersection_index == -1)
    {
      BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(
        intersection_point,
        light_ray.Get_Direction_Vector(),
        new_ray.Get_Direction_Vector());
    }
    else
    {
      BRDF = object_pt_vector[object_intersection_index]->BRDF(
        intersection_point,
        light_ray.Get_Direction_Vector(),
        new_ray.Get_Direction_Vector());
    }

    // Add the radiance from all the objects the light ray hits as it bounces a
    // fixed amount of times
    resulting_light += 2.0 * pi * BRDF * dot(new_direction, normal) *
                       Light_Out(new_ray,
                                 bounces_remaining - 1,
                                 threshold_intersection_distance,
                                 threshold_no_intersection_distance,
                                 finite_difference_size);


    return resulting_light;
  } // End of Light_Out


  Image Render_Image(const unsigned &number_of_bounces,
                     const unsigned &number_of_random_samples,
                     const double &threshold_intersection_distance,
                     const double &threshold_no_intersection_distance,
                     const double &finite_difference_size,
                     const bool &silent = false)
  {
    // Find the resolution of the observer
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Create an image of the correct resolution
    Image image(resolution[0], resolution[1]);

    // Create storage for the radiance at each pixel
    Radiance pixel_radiance(0.0, 0.0, 0.0);

    // Initialise variables used for outputting the current progress to the
    // terminal
    unsigned tenth_percentiles = 0;
    double proportion_done = 0.0;

    // Calculate the radiance of each pixel of the image
    for (unsigned i = 0; i < resolution[0]; i++)
    {
      // This section of code is for outputting the progress to the terminal

      // Don't output progress if 'silent' is true
      if (!silent)
      {
        // Find the proportion of pixels calculated
        proportion_done = double(i + 1) / double(resolution[0]);

        // Every time 10% or more of the pixels have been rendered, print to the
        // terminal
        if (unsigned(10.0 * proportion_done) > tenth_percentiles)
        {
          // tenth_percentiles is used to keep track of the last 10th percent
          // printed
          tenth_percentiles = unsigned(10.0 * proportion_done);

          // The tenth tenth_percentile is 100% so the image will have been
          // rendered
          if (tenth_percentiles == 10)
          {
            std::cout << "The image has been rendered." << std::endl;
          }
          else
          {
            // Update the terminal on the latest 10% done
            std::cout << "Roughly " << 10 * tenth_percentiles
                      << "\% of the pixels have been rendered." << std::endl;
          }
        }
      }

      // Find the radiance at each pixel and set each pixel to this RGB value
      for (unsigned j = 0; j < resolution[1]; j++)
      {
        // Set the radiance at each pixel to zero before calculating the
        // radiance
        pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;

        // Sum up the radiance of multiple light rays "shot out" from the camera
        // in the direction of this pixel
        for (unsigned k = 0; k < number_of_random_samples; k++)
        {
          pixel_radiance += Light_Out(observer_pt->Ray_To_Pixel_XY(i, j),
                                      number_of_bounces,
                                      threshold_intersection_distance,
                                      threshold_no_intersection_distance,
                                      finite_difference_size);
        }

        // Set the RGB value of this pixel to the average radiance over the
        // number of light rays shot out
        image(i, j) = pixel_radiance / number_of_random_samples;
      }
    }

    return image;
  }

  // This function should be used purely in Render_Image_Multithreaded, when a
  // thread is created to do this job, it will find the Radiance of pixels in
  // such a way that it will work in parallel with other threads. The pixel
  // data calculated by all the threads are used in Render_Image_Multithreaded
  // to create the whole picture.
  void Render_Image_Per_Thread(std::vector<Radiance> &partition,
                               const unsigned &number_of_bounces,
                               const unsigned &number_of_random_samples,
                               const double &threshold_intersection_distance,
                               const double &threshold_no_intersection_distance,
                               const double &finite_difference_size,
                               const unsigned &thread_index,
                               const unsigned &number_of_threads)
  {
    // Get the resolution of the image
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Find the number of pixels there are
    unsigned number_of_pixels = resolution[0] * resolution[1];

    // Make sure that the vector of Radiance is empty
    partition.resize(0);

    // Reserve the minimum amount of memory that each vector will need to
    // reduce the cost of push_back
    partition.reserve(number_of_pixels / number_of_threads);

    // If the rows of pixels from top to bottom in an image were laid end to
    // end, pixel_index is the index of a pixel in this "vector" of pixels.
    unsigned pixel_index = thread_index;

    // Initialise variables
    Radiance pixel_radiance(0.0, 0.0, 0.0);
    unsigned pixel_index_i = 0;
    unsigned pixel_index_j = 0;

    while (pixel_index < number_of_pixels)
    {
      // pixel_index_i and pixel_index_j refer to the (i, j)-th pixel of an
      // image, this is calculated from pixel_index and the resolution
      pixel_index_i = pixel_index % resolution[0];
      pixel_index_j = pixel_index / resolution[0];

      // Set the radiance to zero for each pixel before calculating the radiance
      pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;

      // Take the total radiance for a pixel over number_of_random_samples
      // light rays
      for (unsigned i = 0; i < number_of_random_samples; i++)
      {
        pixel_radiance +=
          Light_Out(observer_pt->Ray_To_Pixel_XY(pixel_index_i, pixel_index_j),
                    number_of_bounces,
                    threshold_intersection_distance,
                    threshold_no_intersection_distance,
                    finite_difference_size);
      }

      // Set the RGB value of this pixel to the average radiance over all the
      // rays traced
      partition.push_back(pixel_radiance / number_of_random_samples);

      // Move on to the pixel that is number_of_threads further so that each
      // thread works on a different pixel
      pixel_index += number_of_threads;
    }
  }


  // MULTITHREADING DOES NOT WORK CORRECTLY RIGHT NOW: TO BE FIXED
  // Rupinder: Fix this
  // This function can take in the number of threads to be utilised while
  // rendering the image. This relies on the function of
  // Render_Image_Per_Thread.
  Image Render_Image_Multithreaded(
    const unsigned &number_of_bounces,
    const unsigned &number_of_random_samples,
    const double &threshold_intersection_distance,
    const double &threshold_no_intersection_distance,
    const double &finite_difference_size,
    const unsigned &number_of_threads = std::thread::hardware_concurrency())
  {
    // Create a vector of threads
    std::vector<std::thread> threads;

    // Get the resolution of the image
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Create an image of the correct size
    Image image(resolution[0], resolution[1]);

    // Create the vectors of Radiance that each thread will work on
    std::vector<std::vector<Radiance>> partitions(number_of_threads);

    for (unsigned thread_index = 0; thread_index < number_of_threads;
         thread_index++)
    {
      // Create number_of_threads threads to render the pixels of the image
      threads.push_back(
        std::thread(&SceneRender::Render_Image_Per_Thread,
                    this,
                    std::ref(partitions[thread_index]),
                    std::ref(number_of_bounces),
                    std::ref(number_of_random_samples),
                    std::ref(threshold_intersection_distance),
                    std::ref(threshold_no_intersection_distance),
                    std::ref(finite_difference_size),
                    std::ref(thread_index),
                    std::ref(number_of_threads)));
    }

    // If a thread is joinable, join it
    for (std::thread &t : threads)
    {
      if (t.joinable())
      {
        t.join();
      }
    }

    // Initialise variables
    unsigned pixel_index = 0;
    unsigned pixel_index_i = 0;
    unsigned pixel_index_j = 0;

    // Loop over every pixel in image and find the correct value of Radiance
    // in partitions to assign to this pixel.
    while (pixel_index < resolution[0] * resolution[1])
    {
      // If the rows of pixels from top to bottom of an image were laid end
      // to end, and the index of a pixel in this "vector" of pixels were given
      // as pixel_index, find the (i, j)-th index of this pixel in the image.
      pixel_index_i = pixel_index % resolution[0];
      pixel_index_j = pixel_index / resolution[0];

      // Set the pixel with an index of (pixel_index_i, pixel_index_j) to
      // the correct Radiance value
      image(pixel_index_i, pixel_index_j) =
        partitions[pixel_index % number_of_threads]
                  [pixel_index / number_of_threads];

      // Move on to the next pixel
      pixel_index += 1;
    }

    return image;
  }


  //////////////////////////////////////////////////////////////////////////////
  ///////////////// Equivalent code for Russian Roulette ahead /////////////////
  //////////////////////////////////////////////////////////////////////////////


  // Calculate the radiance coming from the direction of light_ray using the
  // Light Transport Equation with Russian Roulette implemented
  Radiance Light_Out_Russian(const Ray &light_ray,
                             const double &threshold_intersection_distance,
                             const double &threshold_no_intersection_distance,
                             const double &finite_difference_size)
  {
    // Store the radiance that light_ray will "carry"
    Radiance resulting_light(0.0, 0.0, 0.0);

    // Store indices of the surfaces that light_ray will intersect with
    int object_intersection_index = -1;
    int sdf_intersection_index = -1;

    // Storage for a normal vector
    Vec3 normal(0.0, 0.0, 0.0);

    // Storage for the BRDF
    Radiance BRDF(0.0, 0.0, 0.0);

    // Find the closest intersection of light_ray with an object along with the
    // index of the object hit in object_pt_vector
    Vec3 intersection_point =
      First_Intersection_Point(light_ray,
                               object_intersection_index,
                               sdf_intersection_index,
                               threshold_intersection_distance,
                               threshold_no_intersection_distance);

    // If no surface was hit by light_ray, we return zero radiance
    if (object_intersection_index == -1 && sdf_intersection_index == -1)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // Past this point, either object_intersection_index or (exclusive or)
    // sdf_intersection_index is -1.

    // If object_intersection_index is -1, we must have intersected with an SDF,
    // therefore we find the Light_Emitted at the point of intersection.
    if (object_intersection_index == -1)
    {
      resulting_light =
        sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(
          intersection_point);

      // Find the outward normal to the surface hit by light_ray at the point of
      // intersection
      normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(
        intersection_point, finite_difference_size);
    }
    else
    {
      // Otherwise we intersected with a non-SDF surface and we find the light
      // emitted by this surface at the point of intersection.
      resulting_light =
        object_pt_vector[object_intersection_index]->Light_Emitted(
          intersection_point);

      // Find the outward normal to the surface hit by light_ray at the point of
      // intersection
      normal = object_pt_vector[object_intersection_index]->Orientated_Normal(
        intersection_point, -light_ray.Get_Direction_Vector());
    }

    // Find the direction of a random incident ray onto an object, along with
    // its point of origin
    Vec3 new_direction = Hemisphere_Vector_Generator(normal);
    Vec3 new_position = intersection_point + 1.0e-8 * normal;
    Ray new_ray(new_position, new_direction);

    // Find the BRDF of the first object hit by light_ray, the code slightly
    // differs for the different types of surfaces.
    if (object_intersection_index == -1)
    {
      BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(
        intersection_point,
        light_ray.Get_Direction_Vector(),
        new_ray.Get_Direction_Vector());
    }
    else
    {
      BRDF = object_pt_vector[object_intersection_index]->BRDF(
        intersection_point,
        light_ray.Get_Direction_Vector(),
        new_ray.Get_Direction_Vector());
    }

    // Here, we implement the Russian Roulette part of our estimator

    // We choose the termination chance as (1 - the normalised dot product) of
    // the direction of the new ray and the normal at the surface. This is a
    // reasonable function since rays with a more oblique angle to the normal
    // will carry less radiance and contribute less to the overall image so
    // accuracy is not as important.

    // Find the dot product of the new direction of the ray along with the
    // normal to the surface
    double dot_product = dot(new_direction, normal);

    // Create a uniform distribution between 0 and 1.
    std::random_device random_seed;
    std::default_random_engine generator(random_seed());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Enforce the termination probability
    if (distribution(generator) < dot_product)
    {
      // Add the radiance from all the objects the light ray hits as it bounces
      // a fixed amount of times
      resulting_light += 2.0 * pi * BRDF *
                         Light_Out_Russian(new_ray,
                                           threshold_intersection_distance,
                                           threshold_no_intersection_distance,
                                           finite_difference_size);
    }

    return resulting_light;
  } // End of Light_Out_Russian


  Image Render_Image_Russian(const unsigned &number_of_random_samples,
                             const double &threshold_intersection_distance,
                             const double &threshold_no_intersection_distance,
                             const double &finite_difference_size,
                             const bool &silent = false)
  {
    // Find the resolution of the observer
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Create an image of the correct resolution
    Image image(resolution[0], resolution[1]);

    // Create storage for the radiance at each pixel
    Radiance pixel_radiance(0.0, 0.0, 0.0);

    // Initialise variables used for outputting the current progress to the
    // terminal
    unsigned tenth_percentiles = 0;
    double proportion_done = 0.0;

    // Calculate the radiance of each pixel of the image
    for (unsigned i = 0; i < resolution[0]; i++)
    {
      // This section of code is for outputting the progress to the terminal

      // Don't output progress if 'silent' is true
      if (!silent)
      {
        // Find the proportion of pixels calculated
        proportion_done = double(i + 1) / double(resolution[0]);

        // Every time 10% or more of the pixels have been rendered, print to the
        // terminal
        if (unsigned(10.0 * proportion_done) > tenth_percentiles)
        {
          // tenth_percentiles is used to keep track of the last 10th percent
          // printed
          tenth_percentiles = unsigned(10.0 * proportion_done);

          // The tenth tenth_percentile is 100% so the image will have been
          // rendered
          if (tenth_percentiles == 10)
          {
            std::cout << "The image has been rendered." << std::endl;
          }
          else
          {
            // Update the terminal on the latest 10% done
            std::cout << "Roughly " << 10 * tenth_percentiles
                      << "\% of the pixels have been rendered." << std::endl;
          }
        }
      }

      // Find the radiance at each pixel and set each pixel to this RGB value
      for (unsigned j = 0; j < resolution[1]; j++)
      {
        // Set the radiance at each pixel to zero before calculating the
        // radiance
        pixel_radiance.x = 0.0;
        pixel_radiance.y = 0.0;
        pixel_radiance.z = 0.0;

        // Sum up the radiance of multiple light rays "shot out" from the camera
        // in the direction of this pixel
        for (unsigned k = 0; k < number_of_random_samples; k++)
        {
          pixel_radiance +=
            Light_Out_Russian(observer_pt->Ray_To_Pixel_XY(i, j),
                              threshold_intersection_distance,
                              threshold_no_intersection_distance,
                              finite_difference_size);
        }

        // Set the RGB value of this pixel to the average radiance over the
        // number of light rays shot out
        image(i, j) = pixel_radiance / number_of_random_samples;
      }
    }

    return image;
  }

  // This function should be used purely in Render_Image_Russian_Multithreaded,
  // when a thread is created to do this job, it will find the Radiance of
  // pixels in such a way that it will work in parallel with other threads. The
  // pixel data calculated by all the threads are used in
  // Render_Image_Russian_Multithreaded to create the whole picture.
  void Render_Image_Russian_Per_Thread(
    std::vector<Radiance> &partition,
    const unsigned &number_of_random_samples,
    const double &threshold_intersection_distance,
    const double &threshold_no_intersection_distance,
    const double &finite_difference_size,
    const unsigned &thread_index,
    const unsigned &number_of_threads)
  {
    // Get the resolution of the image
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Find the number of pixels there are
    unsigned number_of_pixels = resolution[0] * resolution[1];

    // Make sure that the vector of Radiance is empty
    partition.resize(0);

    // Reserve the minimum amount of memory that each vector will need to
    // reduce the cost of push_back
    partition.reserve(number_of_pixels / number_of_threads);

    // If the rows of pixels from top to bottom in an image were laid end to
    // end, pixel_index is the index of a pixel in this "vector" of pixels.
    unsigned pixel_index = thread_index;

    // Initialise variables
    Radiance pixel_radiance(0.0, 0.0, 0.0);
    unsigned pixel_index_i = 0;
    unsigned pixel_index_j = 0;

    while (pixel_index < number_of_pixels)
    {
      // pixel_index_i and pixel_index_j refer to the (i, j)-th pixel of an
      // image, this is calculated from pixel_index and the resolution
      pixel_index_i = pixel_index % resolution[0];
      pixel_index_j = pixel_index / resolution[0];

      // Set the radiance to zero for each pixel before calculating the radiance
      pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;

      // Take the total radiance for a pixel over number_of_random_samples
      // light rays
      for (unsigned i = 0; i < number_of_random_samples; i++)
      {
        pixel_radiance += Light_Out_Russian(
          observer_pt->Ray_To_Pixel_XY(pixel_index_i, pixel_index_j),
          threshold_intersection_distance,
          threshold_no_intersection_distance,
          finite_difference_size);
      }

      // Set the RGB value of this pixel to the average radiance over all the
      // rays traced
      partition.push_back(pixel_radiance / number_of_random_samples);

      // Move on to the pixel that is number_of_threads further so that each
      // thread works on a different pixel
      pixel_index += number_of_threads;
    }
  }

  // MULTITHREADING DOES NOT WORK CORRECTLY RIGHT NOW: TO BE FIXED
  // Rupinder: Fix this
  // This function can take in the number of threads to be utilised while
  // rendering the image. This relies on the function of
  // Render_Image_Russian_Per_Thread.
  Image Render_Image_Russian_Multithreaded(
    const unsigned &number_of_random_samples,
    const double &threshold_intersection_distance,
    const double &threshold_no_intersection_distance,
    const double &finite_difference_size,
    const unsigned &number_of_threads = std::thread::hardware_concurrency())
  {
    // Create a vector of threads
    std::vector<std::thread> threads;

    // Get the resolution of the image
    std::vector<unsigned> resolution = observer_pt->Get_Resolution();

    // Create an image of the correct size
    Image image(resolution[0], resolution[1]);

    // Create the vectors of Radiance that each thread will work on
    std::vector<std::vector<Radiance>> partitions(number_of_threads);

    for (unsigned thread_index = 0; thread_index < number_of_threads;
         thread_index++)
    {
      // Create number_of_threads threads to render the pixels of the image
      threads.push_back(
        std::thread(&SceneRender::Render_Image_Russian_Per_Thread,
                    this,
                    std::ref(partitions[thread_index]),
                    number_of_random_samples,
                    threshold_intersection_distance,
                    threshold_no_intersection_distance,
                    finite_difference_size,
                    thread_index,
                    number_of_threads));
    }

    // If a thread is joinable, join it
    for (std::thread &t : threads)
    {
      if (t.joinable())
      {
        t.join();
      }
    }

    // Initialise variables
    unsigned pixel_index = 0;
    unsigned pixel_index_i = 0;
    unsigned pixel_index_j = 0;

    // Loop over every pixel in image and find the correct value of Radiance
    // in partitions to assign to this pixel.
    while (pixel_index < resolution[0] * resolution[1])
    {
      // If the rows of pixels from top to bottom of an image were laid end
      // to end, and the index of a pixel in this "vector" of pixels were given
      // as pixel_index, find the (i, j)-th index of this pixel in the image.
      pixel_index_i = pixel_index % resolution[0];
      pixel_index_j = pixel_index / resolution[0];

      // Set the pixel with an index of (pixel_index_i, pixel_index_j) to
      // the correct Radiance value
      image(pixel_index_i, pixel_index_j) =
        partitions[pixel_index % number_of_threads]
                  [pixel_index / number_of_threads];

      // Move on to the next pixel
      pixel_index += 1;
    }

    return image;
  }


private:
  // A vector containing pointers to the physical objects in the scene
  std::vector<std::unique_ptr<PhysicalObject>> object_pt_vector;

  // A vector containing pointers to the physical objects described by SDFs in
  // the scene
  std::vector<std::unique_ptr<SDF>> sdf_object_pt_vector;

  // A pointer to the Observer of the scene
  std::unique_ptr<Observer> observer_pt;

  // Random engine and distribution for use in generating random vectors
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution;
};


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// End of SceneRender ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// Start of Validation Case ////////////////////////////
////////////////////////////////////////////////////////////////////////////////


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

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// End of Validation Case /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Start of Driver Code //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double damping_factor = 0.5;

Radiance red_BRDF(const Vec3 &position,
                  const Vec3 &incident_light_vector,
                  const Vec3 &outgoing_light_vector)
{
  return damping_factor * pi_reciprocal * Radiance(1.0, 0.2, 0.2);
}

Radiance blue_BRDF(const Vec3 &position,
                   const Vec3 &incident_light_vector,
                   const Vec3 &outgoing_light_vector)
{
  return damping_factor * pi_reciprocal * Radiance(0.2, 0.2, 1.0);
}

Radiance white_BRDF(const Vec3 &position,
                    const Vec3 &incident_light_vector,
                    const Vec3 &outgoing_light_vector)
{
  return damping_factor * pi_reciprocal * Radiance(1.0, 1.0, 1.0);
}

Radiance pink_BRDF(const Vec3 &position,
                   const Vec3 &incident_light_vector,
                   const Vec3 &outgoing_light_vector)
{
  return damping_factor * pi_reciprocal * Radiance(1.0, 0.2, 1.0);
}

Radiance purple_BRDF(const Vec3 &position,
                     const Vec3 &incident_light_vector,
                     const Vec3 &outgoing_light_vector)
{
  return damping_factor * pi_reciprocal * Radiance(191.0, 64.0, 191.0) / 255.0;
}

Radiance ceiling_light_emitted(const Vec3 &position)
{
  return Radiance(1.0, 1.0, 1.0);
}


Radiance white_light_emitted(const Vec3 &position)
{
  return Radiance(1.0, 1.0, 1.0);
}


int main()
{
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
}

#endif