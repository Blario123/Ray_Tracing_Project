#include <iostream>

#include "Vec3.h"
#include "Image.h"

int main() {
    // Examples of Vec3 usage
    Vec3 a = {1.0, 2.0, 3.0}, b(4.0, 5.0, 6.0), c({7, 8, 9});
    std::cout << a << " " << b << " " << c << std::endl;
    a.x = 0;
    b.y = 0;
    c.z = 0;
    std::cout << a << " " << b << " " << c << std::endl;
    a += b;
    b -= c;
    c *= a;
    std::cout << a << " " << b << " " << c << std::endl;
    std::cout << a.norm() << " " << a.norm2() << std::endl;
    a.normalise();
    std::cout << a << " " << a.norm() << " " << a.norm2() << std::endl;
    std::cout << -dot(a, b) << " " << dot(-a, b) << " " << cross(a, b) << std::endl;

    // Examples of image class usage
    Image img(60, 40); // 60 x 40 pixel image
    // Image x-coordinage increases from 0 to 59, left to right
    // y-coordinate increases from 0, to 39, downwards from the top of the image
    img(3, 5) = {1, 0.5, 0.5}; // Set pixel at (3,5) to pink
    img(3, 6) = Vec3(0.5, 1, 0.5); // Set pixel at (3,6) to light green
    img.Save("test.png"); // save image as test.png

    return 0;
}
