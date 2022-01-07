#include <iostream>
#include <memory> // Allow usage of unique pointers
#include <thread> // Allow usage of threads
#include <random> // Allow usage of random number generation
#include <chrono> // Allow usage of timing
#include <fstream> // Allow output to file

#include "Image.h" // Include image write implementation
#include "Vec3.h" // Include 3D vectors (Already included in Image.h but kept
// here for readability)
#include "Ray.h" // Includes the Ray class

// Create a type "Radiance" to be a Vec3 containing RGB components
typedef Vec3 Radiance;

// Define the values of pi and the reciprocal of pi (To make dividing by pi
// faster)
const double pi = M_PI;
const double pi_reciprocal = M_1_PI;

// The base class for physical objects present in the scene
class PhysicalObject
{
public:
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
#ifdef TEST
    if (radius_ <= 0.0)
    {
      throw std::invalid_argument(
        "Can not create a sphere with a non-positive radius");
    }
#endif

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

#ifdef TEST
    // Check that the direction vector of the ray is normalised
    if (direction_vector.norm() - 1.0 > 1.0e-8)
    {
      throw std::invalid_argument(
        "The Ray given has a non-normalised direction vector");
    }
#endif

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

#ifdef TEST
    // Check that the position given is on the surface of the sphere
    if (normal.norm2() - radius * radius > 1.0e-8)
    {
      throw std::invalid_argument(
        "The position given must be on the sphere's surface.");
    }
#endif

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

    // The horizontal viewing angle
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
    Vec3 direction_to_pixel =
      unit_camera_direction +
      (2.0 * tan(horizontal_field_of_view_angle / 2.0) / resolution[0]) *
        (x_pixel - (resolution[0] / 2.0 - 0.5)) * unit_sideways_direction +
      (2.0 * tan(vertical_field_of_view_angle / 2.0) / resolution[1]) *
        ((resolution[1] / 2.0 - 0.5) - y_pixel) * unit_upward_direction;

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


  // Random hemisphere vector generator
  Vec3 Hemisphere_Vector_Generator(const Vec3 &normal)
  {
    // Get a seed with a random value
    std::random_device random_seed;

    // Create a uniform distribution between -1 and 1
    std::default_random_engine generator(random_seed());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

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
  } // End of Hemisphere_Vector_Generator


  // Loop over all the objects in the scene and check which object was
  // intersected with first by the argument Ray "ray".
  Vec3 First_Intersection_Point(const Ray &light_ray,
                                int &object_intersection_index)
  {
    // If there is no intersection between light_ray and an object, the index
    // will stay as -1.
    object_intersection_index = -1;

    // A Boolean showing whether light_ray has been found to intersect with an
    // object in the scene. Used for knowing when to update smallest_distance.
    bool found_an_intersection = false;

    // Store the smallest distance from the light ray source to an intersection
    // along with the index of the corresponding object in object_pt_vector.
    double smallest_distance = 0.0;

    // Store the distance of the light ray to the intersection with the current
    // object being looped over.
    double current_distance = 0.0;

    // Loop over every object in the scene to find an intersection
    for (unsigned i = 0; i < object_pt_vector.size(); i++)
    {
      // If an intersection has already been found, check whether this new
      // intersection is closer to the light ray source than the previous
      // closest intersection.
      if (found_an_intersection)
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

          // An intersection has been found
          found_an_intersection = true;
        }
      }
    }

    // If there hasn't been an intersection, return the zero vector
    if (object_intersection_index == -1)
    {
      return Vec3(0.0, 0.0, 0.0);
    }

    // Calculate the position of closest intersection between a light ray and
    // any object in the scene.
    Vec3 vector = light_ray.Get_Initial_Position() +
                  smallest_distance * light_ray.Get_Direction_Vector();

    return vector;
  } // End of First_Intersection_Point


  // Calculate the radiance coming from the direction of light_ray using the
  // Light Transport Equation considering the light rays take bounces_remaining
  // number of bounces.
  Radiance Light_Out(const Ray &light_ray, unsigned bounces_remaining)
  {
    // If the light ray can't bounce off a single object, no light will reach
    // the "observer".
    if (bounces_remaining == 0)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // Initialise the index of the first object in object_pt_vector hit by
    // light_ray
    int index_of_object_hit = 0;

    // Find the closest intersection of light_ray with an object along with the
    // index of the object hit in object_pt_vector
    Vec3 intersection_point =
      First_Intersection_Point(light_ray, index_of_object_hit);

    // If no object was hit according to First_Intersection_Point,
    // index_of_object_hit will be -1 and therefore no light will be seen
    if (index_of_object_hit == -1)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // First, find the light emitted by the object in the direction of light_ray
    Radiance resulting_light =
      object_pt_vector[index_of_object_hit]->Light_Emitted(intersection_point);

    // Find the normal to the surface hit by light_ray at the point of
    // intersection
    Vec3 normal = object_pt_vector[index_of_object_hit]->Orientated_Normal(
      intersection_point, -light_ray.Get_Direction_Vector());

    // Find the direction of a random incident ray onto an object, along with
    // its point of origin
    Vec3 new_direction = Hemisphere_Vector_Generator(normal);
    Vec3 new_position = intersection_point + 1.0e-8 * normal;
    Ray new_ray(new_position, new_direction);

    // Find the BRDF of the first object hit by light_ray
    Vec3 brdf = object_pt_vector[index_of_object_hit]->BRDF(
      intersection_point,
      light_ray.Get_Direction_Vector(),
      new_ray.Get_Direction_Vector());

    // Add the radiance from all the objects the light ray hits as it bounces a
    // fixed amount of times
    resulting_light += (2.0 * pi * brdf) * dot(new_direction, normal) *
                       Light_Out(new_ray, bounces_remaining - 1);


    return resulting_light;
  } // End of Light_Out


  // Calculate the radiance coming from the direction of light_ray using the
  // Light Transport Equation with Russian Roulette implemented
  Radiance Light_Out_Russian(const Ray &light_ray)
  {
    // Initialise the index of the first object in object_pt_vector hit by
    // light_ray
    int index_of_object_hit = 0;

    // Find the closest intersection of light_ray with an object along with the
    // index of the object hit in object_pt_vector
    Vec3 intersection_point =
      First_Intersection_Point(light_ray, index_of_object_hit);

    // If no object was hit according to First_Intersection_Point,
    // index_of_object_hit will be -1 and therefore no light will be seen
    if (index_of_object_hit == -1)
    {
      return Radiance(0.0, 0.0, 0.0);
    }

    // First, find the light emitted by the object in the direction of light_ray
    Radiance resulting_light =
      object_pt_vector[index_of_object_hit]->Light_Emitted(intersection_point);

    // Find the normal to the surface hit by light_ray at the point of
    // intersection
    Vec3 normal = object_pt_vector[index_of_object_hit]->Orientated_Normal(
      intersection_point, -light_ray.Get_Direction_Vector());

    // Find the direction of a random incident ray onto an object, along with
    // its point of origin
    Vec3 new_direction = Hemisphere_Vector_Generator(normal);
    Vec3 new_position = intersection_point + 1.0e-8 * normal;
    Ray new_ray(new_position, new_direction);

    // Choose the termination chance of a bounce based on the reflected
    // direction
    double reflected_angle = dot(new_direction, normal);

    // Create a uniform distribution between 0 and 1.
    std::random_device random_seed;
    std::default_random_engine generator(random_seed());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Find the BRDF of the first object hit by light_ray
    Vec3 brdf = object_pt_vector[index_of_object_hit]->BRDF(
      intersection_point,
      light_ray.Get_Direction_Vector(),
      new_ray.Get_Direction_Vector());

    // Enforce a termination probability on the ray.
    if (distribution(generator) < reflected_angle)
    {
      // Add the radiance from all the objects the light ray hits as it bounces
      // a fixed amount of times
      resulting_light += (2.0 * pi * brdf) * Light_Out_Russian(new_ray);
    }

    return resulting_light;
  } // End of Light_Out


  // This function should be used purely in Render_Image_Multithreaded, when a
  // thread is created to do this job, it will find the Radiance of pixels in
  // such a way that it will work in parallel with other threads. The pixel
  // data calculated by all the threads are used in Render_Image_Multithreaded
  // to create the whole picture.
  void Render_Image_Per_Thread(std::vector<Radiance> &partition,
                               const unsigned number_of_bounces,
                               const unsigned number_of_random_samples,
                               const unsigned thread_index,
                               const unsigned number_of_threads)
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
                    number_of_bounces);
      }

      // Set the RGB value of this pixel to the average radiance over all the
      // rays traced
      partition.push_back(pixel_radiance / number_of_random_samples);

      // Move on to the pixel that is number_of_threads further so that each
      // thread works on a different pixel
      pixel_index += number_of_threads;
    }
  }


  // This function should be used purely in Render_Image_Multithreaded_Russian,
  // when a thread is created to do this job, it will find the Radiance of
  // pixels in such a way that it will work in parallel with other threads. The
  // pixel data calculated by all the threads are used in
  // Render_Image_Multithreaded_Russian to create the whole picture.
  void Render_Image_Per_Thread_Russian(std::vector<Radiance> &partition,
                                       const unsigned number_of_random_samples,
                                       const unsigned thread_index,
                                       const unsigned number_of_threads)
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
          observer_pt->Ray_To_Pixel_XY(pixel_index_i, pixel_index_j));
      }

      // Set the RGB value of this pixel to the average radiance over all the
      // rays traced
      partition.push_back(pixel_radiance / number_of_random_samples);

      // Move on to the pixel that is number_of_threads further so that each
      // thread works on a different pixel
      pixel_index += number_of_threads;
    }
  }


  // This function can take in the number of threads to be utilised while
  // rendering the image. This relies on the function of
  // Render_Image_Per_Thread.
  Image Render_Image_Multithreaded(
    const unsigned &number_of_bounces,
    const unsigned &number_of_random_samples,
    const unsigned &number_of_threads = std::thread::hardware_concurrency())
  {
#ifdef TEST
    // If no random samples are taken, no image will be produced and a division
    // by zero may occur so we throw an error
    if (number_of_random_samples == 0)
    {
      throw std::invalid_argument(
        "The number of random samples may not be chosen as zero");
    }

    // Check if the number of threads given as an argument is less than or
    // equal to the number of threads available. If not, reduce
    // number_of_threads to the maximum possible value.
    if (number_of_threads > std::thread::hardware_concurrency())
    {
      throw std::invalid_argument(
        "You don't have the number of threads specified in "
        "Render_Image_Multithreaded available.");
    }
#endif

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
      threads.push_back(std::thread(&SceneRender::Render_Image_Per_Thread,
                                    this,
                                    std::ref(partitions[thread_index]),
                                    number_of_bounces,
                                    number_of_random_samples,
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


  // This function can take in the number of threads to be utilised while
  // rendering the image. This relies on the function of
  // Render_Image_Per_Thread_Russian.
  Image Render_Image_Multithreaded_Russian(
    const unsigned &number_of_random_samples,
    const unsigned &number_of_threads = std::thread::hardware_concurrency())
  {
#ifdef TEST
    // If no random samples are taken, no image will be produced and a division
    // by zero may occur so we throw an error
    if (number_of_random_samples == 0)
    {
      throw std::invalid_argument(
        "The number of random samples may not be chosen as zero");
    }

    // Check if the number of threads given as an argument is less than or
    // equal to the number of threads available. If not, reduce
    // number_of_threads to the maximum possible value.
    if (number_of_threads > std::thread::hardware_concurrency())
    {
      throw std::invalid_argument(
        "You don't have the number of threads specified in "
        "Render_Image_Multithreaded available.");
    }
#endif

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
        std::thread(&SceneRender::Render_Image_Per_Thread_Russian,
                    this,
                    std::ref(partitions[thread_index]),
                    number_of_random_samples,
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


  // Render the image and export it to filename.
  Image Render_Image(const unsigned &number_of_bounces,
                     const unsigned &number_of_random_samples)
  {
#ifdef TEST
    // If no random samples are taken, no image will be produced and a division
    // by zero may occur so we throw an error
    if (number_of_random_samples == 0)
    {
      throw std::invalid_argument(
        "The number of random samples may not be chosen as zero");
    }
#endif
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
          pixel_radiance +=
            Light_Out(observer_pt->Ray_To_Pixel_XY(i, j), number_of_bounces);
        }

        // Set the RGB value of this pixel to the average radiance over the
        // number of light rays shot out
        image(i, j) = pixel_radiance / number_of_random_samples;
      }
    }

    return image;
  }

private:
  // A vector containing pointers to the physical objects in the scene
  std::vector<std::unique_ptr<PhysicalObject>> object_pt_vector;

  // A pointer to the Observer of the scene
  std::unique_ptr<Observer> observer_pt;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////// End of SceneRender ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Validation case
// Create a light emitted function for the validation case (Only the x_value
// shall be non-zero for convenience)
Radiance validation_light_emitted(const Vec3 &position)
{
  // Return this random RGB radiance value
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
  return Vec3(0.25 * pi_reciprocal, 0.0, 0.0);
}

Radiance red_BRDF(const Vec3 &position,
                  const Vec3 &incident_light_vector,
                  const Vec3 &outgoing_light_vector)
{
  return pi_reciprocal * Radiance(1.0, 0.0, 0.0);
}

Radiance blue_BRDF(const Vec3 &position,
                   const Vec3 &incident_light_vector,
                   const Vec3 &outgoing_light_vector)
{
  return pi_reciprocal * Radiance(0.0, 0.0, 1.0);
}

Radiance white_BRDF(const Vec3 &position,
                    const Vec3 &incident_light_vector,
                    const Vec3 &outgoing_light_vector)
{
  return pi_reciprocal * Radiance(1.0, 1.0, 1.0);
}

Radiance pink_BRDF(const Vec3 &position,
                   const Vec3 &incident_light_vector,
                   const Vec3 &outgoing_light_vector)
{
  return pi_reciprocal * Radiance(1.0, 0.0, 1.0);
}

Radiance purple_BRDF(const Vec3 &position,
                     const Vec3 &incident_light_vector,
                     const Vec3 &outgoing_light_vector)
{
  return pi_reciprocal * Radiance(191.0, 64.0, 191.0) / 255.0;
}

Radiance ceiling_light_emitted(const Vec3 &position)
{
  if (position.x < 0.385 && position.x > 0.315 && position.y < 0.035 &&
      position.y > -0.035)
  {
    return Radiance(3.0, 3.0, 3.0);
  }
  else
  {
    return Radiance(0.0, 0.0, 0.0);
  }
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
                               std::vector<unsigned>(2, 1.0));

  SceneRender validation_scene(validation_observer);

  // Create a single sphere for the validation case
  Sphere validation_sphere(Vec3(0.0, 0.0, 0.0), 1.0);

  // Set the light emitted by this validation sphere
  validation_sphere.Light_Emitted_Fct_Pt = validation_light_emitted;

  // Set the BRDF of the validation sphere to a Lambertian BRDF
  validation_sphere.BRDF_Fct_Pt = validation_BRDF;

  // Add the validation sphere to the validation scene
  validation_scene.Add_Object(std::make_unique<Sphere>(validation_sphere));

  // Find the analytical solution to the LTE for the validation case
  double result = validation_light_emitted(Vec3()).x /
                  (1.0 - pi * validation_BRDF(Vec3(), Vec3(), Vec3()).x);

  // Store variables used in working out the variance of the "empirical" method
  double variance = 0.0;
  unsigned number_of_variance_samples = 100;
  unsigned max_camera_samples = 1000;

  // Open a file to store results in
  std::ofstream validation_results_file;
  validation_results_file.open("Validation_results.dat");

  // Output the name of each column at the top of the file
  validation_results_file.width(20);
  validation_results_file << "No. camera samples";
  validation_results_file.width(20);
  validation_results_file << "Variance" << std::endl;
  validation_results_file << std::fixed;

  // Loop over a range of number of camera samples
  for (unsigned i = 1; i < (max_camera_samples / 10); i++)
  {
    // Reset the variance after each value of number of camera samples worked on
    variance = 0.0;

    // Find the variance over number_of_variance_samples number of samples
    for (unsigned j = 0; j < number_of_variance_samples; j++)
    {
      // Inside of sum of variance formula
      variance += pow(
        validation_scene.Render_Image_Multithreaded_Russian(i)(0, 0).x - result,
        2);
    }
    // Divide by number of samples to get the variance
    variance /= number_of_variance_samples;

    // Output the number of camera samples used and the respective variance in
    // separate columns
    validation_results_file.width(20);
    validation_results_file << 10 * i;
    validation_results_file.width(20);
    validation_results_file.precision(16);
    validation_results_file << variance << std::endl;
  }

  validation_results_file.close();

  std::vector<unsigned> resolution;
  resolution.push_back(128);
  resolution.push_back(108);

  Observer observer(Vec3(0.0, 0.0, 1.0),
                    Vec3(1.0, 0.0, 0.0),
                    Vec3(0.0, 0.0, 1.0),
                    pi / 3.0,
                    resolution);

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

  // scene.Render_Image_Multithreaded_Russian(1).Save("Cornell_Box_3.png");
}