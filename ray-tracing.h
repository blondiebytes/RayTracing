#include <list>
#include <cstdlib>

using namespace std;

class Color;
class Light;
class Figure;
class Plane;
class Sphere;

class Vec
{
  private:
    double x;
    double y;
    double z;

  public:
    Vec();
	Vec::Vec(double x, double y, double z);
    Vec(ifstream& ifs);
	double dot(Vec v);

};


class Color
{
  friend Color operator*(double num, const Color& c);

  protected:
   double red, green, blue;

  public:
	Color();
    Color(ifstream& ifs);
    Color(double r, double g, double b);
	Color add(Color c);
	Color scale(Color c);

};

class Light
{
  public:
    Light(ifstream& ifs);
	Vec getPosition();
  private:
    Vec position;
    Color shading;
    double c0, c1, c2;
};

class Ray {
	public:
		Ray(Vec a, Vec b);
	private:
		Vec p0; 
		Vec p1;
	//	Color color;

};


class Figure
{
 protected:
   Color ambient;
   Color diffuse;
   Color specular;
   Color reflectivity;
   Color transmissivity;
   double shininess;
   double indexOfRefraction;
   int rFlag, tFlag;

 public:
   Figure();
   void initFigure(ifstream& ifs);
};



class Plane : public Figure
{
  private:
    Vec abcVector;
    double dScalar;
    Vec direction1;
    Vec direction2;
  public:
    Plane(ifstream& ifs);
};

class Sphere : public Figure
{
  private:
    Vec center;
    double radius;
  public:
    Sphere(ifstream& ifs);
};

