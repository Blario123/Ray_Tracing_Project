#ifndef PHYSICALOBJECT_H
#define PHYSICALOBJECT_H

#include "RayTracingProject.h"

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

#endif //PHYSICALOBJECT_H
