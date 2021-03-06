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
	Vec(double x, double y, double z);
    Vec(ifstream& ifs);
	double dot(Vec v);
	Vec normalize();
	double getLength();
	Vec divide(double v);
	Vec subtract(Vec v);
	Vec add(Vec v);
	Vec scale(double t);

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
	double getRed();
	double getBlue();
	double getGreen();
	Color add(Color c);
	Color multiply(Color c);
	Color scale(double n);
	bool isEqual(Color c);
	string toString();
	Color cut();

};

class Light
{
  public:
    Light(ifstream& ifs);
	Vec getPosition();
	Color getShading();
  private:
    Vec position;
    Color shading;
    double c0, c1, c2;
};

class Ray {
	public:
		Ray(Vec a, Vec b);
		Vec getP0();
		Vec getP1();
		Plane* getTravelingThroughPlane();
		Sphere* getTravelingThroughSphere();
		bool travelingThoughAnything();
	private:
		Vec p0; 
		Vec p1;
		Plane* travelingThroughPlane;
		Sphere* travelingThroughSphere;
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
   Color getTransmissivity();
   Color getReflectivity();
   Color getColorAmbient();
   Color getColorDiffuse();
   Color getColorSpecular();
   double getShininess();
   virtual double intersection(const Ray& r, double minT, double maxT) const = 0;
   virtual Vec* getNormal(Vec* i) = 0;
   
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
	Vec * getNormal(Vec * i);
	double Plane::intersection(const Ray& r, double minT, double maxT) const;
};

class Sphere : public Figure
{
  private:
    Vec center;
    double radius;
  public:
    Sphere(ifstream& ifs);
	Vec * getNormal(Vec * i);
	double Sphere::intersection(const Ray& r, double minT, double maxT) const;
};

Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth);

pair<double, Figure*> nearestIntersection(const Ray& r,
	double minT, double maxT,
	bool mayBeTransparent = true);

