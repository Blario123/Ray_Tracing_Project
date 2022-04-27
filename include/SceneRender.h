#ifndef SCENERENDER_H
#define SCENERENDER_H

#include "SDF.h"

class SceneRender {
public:
	// Constructor
	SceneRender(Observer &observer_) {
		// Only an observer is required in the construction of a scene
		observer_pt = std::make_unique<Observer>(observer_);
		
		// Get a seed with a random value
		std::random_device random_seed;
		
		// Create a uniform distribution between 0 and 1
		generator = std::default_random_engine(random_seed());
		distribution = std::uniform_real_distribution<double>(0.0, 1.0);
	}
	
	SceneRender(DOFObserver &observer_) {
		// Only an observer is required in the construction of a scene
		observer_pt = std::make_unique<DOFObserver>(observer_);
		
		// Get a seed with a random value
		std::random_device random_seed;
		
		// Create a uniform distribution between 0 and 1
		generator = std::default_random_engine(random_seed());
		distribution = std::uniform_real_distribution<double>(0.0, 1.0);
	}
	
	// Add objects to the scene via a unique pointer
	void Add_Object(std::unique_ptr<PhysicalObject> &object_upt) {
		object_pt_vector.push_back(std::move(object_upt));
	}
	
	// The same function as above but now accepting rvalues
	void Add_Object(std::unique_ptr<PhysicalObject> &&object_upt) {
		object_pt_vector.push_back(std::move(object_upt));
	}
	
	// Add SDF defined objects to the scene via a unique pointer
	void Add_Object(std::unique_ptr<SDF> &object_upt) {
		sdf_object_pt_vector.push_back(std::move(object_upt));
	}
	
	// The same function as above but now accepting rvalues
	void Add_Object(std::unique_ptr<SDF> &&object_upt) {
		sdf_object_pt_vector.push_back(std::move(object_upt));
	}
	
	// Random hemisphere vector generator
	Vec3 Hemisphere_Vector_Generator(const Vec3 &normal) {
		// The random vector is instantiated first
		Vec3 random_vector;
		
		// Keep looping over this algorithm until a suitable vector is found
		while (true) {
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
			random_vector.x = 2.0 * distribution(generator) - 1.0;
			random_vector.y = 2.0 * distribution(generator) - 1.0;
			random_vector.z = 2.0 * distribution(generator) - 1.0;
			
			// Find the squared norm
			double modulus = random_vector.norm2();
			
			// Check the random vector is inside the unit sphere but isn't the zero
			// vector
			if (modulus <= 1.0 && modulus > 0.0) {
				// Normalise the random vector
				random_vector.normalise();
				
				// Make sure that the random vector is in the hemisphere defined by the
				// direction of the normal
				if (dot(random_vector, normal) < 0.0) {
					random_vector = -random_vector;
				}
				break;
			}
		}
		
		// Return the random vector
		return random_vector;
	} // End of Hemisphere_Vector_Generator
	
	// Take the union of the multiple SDFs provided, by taking the minimum value
	double SDF_Fct(const Vec3 &position) const {
		// Allocate storage for the minimum SDF value and the value of the current
		// SDF as we loop through them
		double minimum_distance = 0.0;
		double current_distance = 0.0;
		
		if (sdf_object_pt_vector.size() > 0) {
			// Set the minimum SDF value equal to the SDF of the first SDF provided
			minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
		}
		
		// Loop through the rest of the SDFs and find the minimum value of the SDFs
		// at "position"
		for (unsigned i = 1; i < sdf_object_pt_vector.size(); i++) {
			// If we find a smaller SDF value, set minimum_distance to that value
			double current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
			if (current_distance < minimum_distance) {
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
	double SDF_Fct(const Vec3 &position, int &sdf_intersection_index) const {
		// Allocate storage for the minimum SDF value and the value of the current
		// SDF as we loop through them
		double minimum_distance = 0.0;
		double current_distance = 0.0;
		
		if (sdf_object_pt_vector.size() > 0) {
			// Set the minimum SDF value equal to the SDF of the first SDF provided,
			// and
			// set the value of sdf_intersection_index to the index of this first SDF.
			minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
			sdf_intersection_index = 0;
		}
		
		// Loop through the rest of the SDFs and find the minimum value of the SDFs
		// at "position". Also set sdf_intersection_index to the corresponding index
		// value.
		for (unsigned i = 1; i < sdf_object_pt_vector.size(); i++) {
			// If we find a smaller SDF value, set minimum_distance and
			// sdf_intersection_index
			current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
			if (current_distance < minimum_distance) {
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
			const double &threshold_no_intersection_distance) {
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
		for (unsigned i = 0; i < object_pt_vector.size(); i++) {
			// If an intersection has already been found, check whether this new
			// intersection is closer to the light ray source than the previous
			// closest intersection.
			if (object_intersection_index != -1) {
				// Check whether this intersection is closer than the previous closest
				// one. If so, replace smallest_distance
				if (object_pt_vector[i]->Intersection_Check(light_ray,
															current_distance) &&
					current_distance < smallest_distance) {
					smallest_distance = current_distance;
					object_intersection_index = i;
				}
			}
				// If an intersection hasn't been found yet, any intersection will be the
				// closest intersection so far.
			else {
				if (object_pt_vector[i]->Intersection_Check(light_ray,
															current_distance)) {
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
			safe_travel_distance < smallest_distance) {
			// Travel along the light_ray while the light_ray has not intersected with
			// an SDF, has not gotten too far away from a surface, or has already
			// passed through a non-SDF defined surface
			while (safe_travel_distance >= threshold_intersection_distance &&
				   safe_travel_distance <= threshold_no_intersection_distance &&
				   (distance_travelled < smallest_distance ||
					object_intersection_index == -1)) {
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
			if (safe_travel_distance < threshold_intersection_distance) {
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
		if (object_intersection_index != -1) {
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
					   const double &finite_difference_size) {
		// If the light ray can't bounce off a single object, no light will reach
		// the "observer".
		if (bounces_remaining == 0) {
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
		if (object_intersection_index == -1 && sdf_intersection_index == -1) {
			return Radiance(0.0, 0.0, 0.0);
		}
		
		// Past this point, either object_intersection_index or (exclusive or)
		// sdf_intersection_index is -1.
		
		// If object_intersection_index is -1, we must have intersected with an SDF,
		// therefore we find the Light_Emitted at the point of intersection.
		if (object_intersection_index == -1) {
			resulting_light =
					sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(
							intersection_point);
			
			// Find the outward normal to the surface hit by light_ray at the point of
			// intersection
			normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(
					intersection_point, finite_difference_size);
		} else {
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
		if (object_intersection_index == -1) {
			BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(
					intersection_point,
					light_ray.Get_Direction_Vector(),
					new_ray.Get_Direction_Vector());
		} else {
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
					   const bool &silent = false) {
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
		for (unsigned i = 0; i < resolution[0]; i++) {
			// This section of code is for outputting the progress to the terminal
			
			// Don't output progress if 'silent' is true
			if (!silent) {
				// Find the proportion of pixels calculated
				proportion_done = double(i + 1) / double(resolution[0]);
				
				// Every time 10% or more of the pixels have been rendered, print to the
				// terminal
				if (unsigned(10.0 * proportion_done) > tenth_percentiles) {
					// tenth_percentiles is used to keep track of the last 10th percent
					// printed
					tenth_percentiles = unsigned(10.0 * proportion_done);
					
					// The tenth tenth_percentile is 100% so the image will have been
					// rendered
					if (tenth_percentiles == 10) {
						std::cout << "The image has been rendered." << std::endl;
					} else {
						// Update the terminal on the latest 10% done
						std::cout << "Roughly " << 10 * tenth_percentiles
								  << "\% of the pixels have been rendered." << std::endl;
					}
				}
			}
			
			// Find the radiance at each pixel and set each pixel to this RGB value
			for (unsigned j = 0; j < resolution[1]; j++) {
				// Set the radiance at each pixel to zero before calculating the
				// radiance
				pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;
				
				// Sum up the radiance of multiple light rays "shot out" from the camera
				// in the direction of this pixel
				for (unsigned k = 0; k < number_of_random_samples; k++) {
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
								 const unsigned &number_of_threads) {
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
		
		while (pixel_index < number_of_pixels) {
			// pixel_index_i and pixel_index_j refer to the (i, j)-th pixel of an
			// image, this is calculated from pixel_index and the resolution
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			
			// Set the radiance to zero for each pixel before calculating the radiance
			pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;
			
			// Take the total radiance for a pixel over number_of_random_samples
			// light rays
			for (unsigned i = 0; i < number_of_random_samples; i++) {
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
	
	
	// This function can take in the number of threads to be utilised while
	// rendering the image. This relies on the function of
	// Render_Image_Per_Thread.
	Image Render_Image_Multithreaded(
			const unsigned &number_of_bounces,
			const unsigned &number_of_random_samples,
			const double &threshold_intersection_distance,
			const double &threshold_no_intersection_distance,
			const double &finite_difference_size,
			const unsigned &number_of_threads = std::thread::hardware_concurrency()) {
		// Create a vector of threads
		std::vector<std::thread> threads;
		
		// Get the resolution of the image
		std::vector<unsigned> resolution = observer_pt->Get_Resolution();
		
		// Create an image of the correct size
		Image image(resolution[0], resolution[1]);
		
		// Create the vectors of Radiance that each thread will work on
		std::vector<std::vector<Radiance>> partitions(number_of_threads);
		
		for (unsigned thread_index = 0; thread_index < number_of_threads;
			 thread_index++) {
			// Create number_of_threads threads to render the pixels of the image
			threads.push_back(std::thread(&SceneRender::Render_Image_Per_Thread,
										  this,
										  std::ref(partitions[thread_index]),
										  number_of_bounces,
										  number_of_random_samples,
										  threshold_intersection_distance,
										  threshold_no_intersection_distance,
										  finite_difference_size,
										  thread_index,
										  number_of_threads));
		}
		
		// If a thread is joinable, join it
		for (std::thread &t: threads) {
			if (t.joinable()) {
				t.join();
			}
		}
		
		// Initialise variables
		unsigned pixel_index = 0;
		unsigned pixel_index_i = 0;
		unsigned pixel_index_j = 0;
		
		// Loop over every pixel in image and find the correct value of Radiance
		// in partitions to assign to this pixel.
		while (pixel_index < resolution[0] * resolution[1]) {
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
	
	// Calculate the radiance coming from the direction of light_ray using the
	// Light Transport Equation with Russian Roulette implemented
	Radiance Light_Out_Russian(const Ray &light_ray,
							   const double &threshold_intersection_distance,
							   const double &threshold_no_intersection_distance,
							   const double &finite_difference_size) {
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
		if (object_intersection_index == -1 && sdf_intersection_index == -1) {
			return Radiance(0.0, 0.0, 0.0);
		}
		
		// Past this point, either object_intersection_index or (exclusive or)
		// sdf_intersection_index is -1.
		
		// If object_intersection_index is -1, we must have intersected with an SDF,
		// therefore we find the Light_Emitted at the point of intersection.
		if (object_intersection_index == -1) {
			resulting_light =
					sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(
							intersection_point);
			
			// Find the outward normal to the surface hit by light_ray at the point of
			// intersection
			normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(
					intersection_point, finite_difference_size);
		} else {
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
		if (object_intersection_index == -1) {
			BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(
					intersection_point,
					light_ray.Get_Direction_Vector(),
					new_ray.Get_Direction_Vector());
		} else {
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
		
		// Enforce the termination probability
		if (distribution(generator) < dot_product) {
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
							   const bool &silent = false) {
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
		for (unsigned i = 0; i < resolution[0]; i++) {
			// This section of code is for outputting the progress to the terminal
			
			// Don't output progress if 'silent' is true
			if (!silent) {
				// Find the proportion of pixels calculated
				proportion_done = double(i + 1) / double(resolution[0]);
				
				// Every time 10% or more of the pixels have been rendered, print to the
				// terminal
				if (unsigned(10.0 * proportion_done) > tenth_percentiles) {
					// tenth_percentiles is used to keep track of the last 10th percent
					// printed
					tenth_percentiles = unsigned(10.0 * proportion_done);
					
					// The tenth tenth_percentile is 100% so the image will have been
					// rendered
					if (tenth_percentiles == 10) {
						std::cout << "The image has been rendered." << std::endl;
					} else {
						// Update the terminal on the latest 10% done
						std::cout << "Roughly " << 10 * tenth_percentiles
								  << "\% of the pixels have been rendered." << std::endl;
					}
				}
			}
			
			// Find the radiance at each pixel and set each pixel to this RGB value
			for (unsigned j = 0; j < resolution[1]; j++) {
				// Set the radiance at each pixel to zero before calculating the
				// radiance
				pixel_radiance.x = 0.0;
				pixel_radiance.y = 0.0;
				pixel_radiance.z = 0.0;
				
				// Sum up the radiance of multiple light rays "shot out" from the camera
				// in the direction of this pixel
				for (unsigned k = 0; k < number_of_random_samples; k++) {
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
			const unsigned &number_of_threads) {
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
		
		while (pixel_index < number_of_pixels) {
			// pixel_index_i and pixel_index_j refer to the (i, j)-th pixel of an
			// image, this is calculated from pixel_index and the resolution
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			
			// Set the radiance to zero for each pixel before calculating the radiance
			pixel_radiance.x = pixel_radiance.y = pixel_radiance.z = 0.0;
			
			// Take the total radiance for a pixel over number_of_random_samples
			// light rays
			for (unsigned i = 0; i < number_of_random_samples; i++) {
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

	// This function can take in the number of threads to be utilised while
	// rendering the image. This relies on the function of
	// Render_Image_Russian_Per_Thread.
	Image Render_Image_Russian_Multithreaded(
			const unsigned &number_of_random_samples,
			const double &threshold_intersection_distance,
			const double &threshold_no_intersection_distance,
			const double &finite_difference_size,
			const unsigned &number_of_threads = std::thread::hardware_concurrency()) {
		// Create a vector of threads
		std::vector<std::thread> threads;
		
		// Get the resolution of the image
		std::vector<unsigned> resolution = observer_pt->Get_Resolution();
		
		// Create an image of the correct size
		Image image(resolution[0], resolution[1]);
		
		// Create the vectors of Radiance that each thread will work on
		std::vector<std::vector<Radiance>> partitions(number_of_threads);
		
		for (unsigned thread_index = 0; thread_index < number_of_threads;
			 thread_index++) {
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
		for (std::thread &t: threads) {
			if (t.joinable()) {
				t.join();
			}
		}
		
		// Initialise variables
		unsigned pixel_index = 0;
		unsigned pixel_index_i = 0;
		unsigned pixel_index_j = 0;
		
		// Loop over every pixel in image and find the correct value of Radiance
		// in partitions to assign to this pixel.
		while (pixel_index < resolution[0] * resolution[1]) {
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

#endif //SCENERENDER_H
