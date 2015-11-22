#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <list>
#include <iostream>

using namespace std;

#include "ray-tracing.h"

list<Figure*> shapeList;
list<Light*> lightList;



double cameraX, cameraY, cameraZ;
int horizontalResolution, verticalResolution;
double zCoor;
double minX, maxX;
double minY, maxY;
double win_h = maxX - minX;
double win_w = maxY - minY;

int maxDepth;
Color backgroundColor;
Color ambient;

Vec::Vec() : x(0.0),y(0.0),z(0.0) {}

Vec::Vec(double x, double y, double z) : x(x), y(y), z(z) {}

Vec::Vec(ifstream& ifs)
{
  ifs >> x >> y >> z;
}

double Vec::dot(Vec v) {
	return ((x * v.x) + (y * v.y) + (z * v.z));
}

 Vec Vec::normalize() {
	 return this->divide(this->getLength());
}

Vec Vec::divide(double b) {
	return Vec(x / b, y / b, z / b);
}

Vec Vec::subtract(Vec v) {
	return Vec(x - v.x, y - v.y, z - v.z);
}

double Vec::getLength() {
	return sqrt((x * x) + (y * y) + (z * z));
}


Color::Color(): red(0.0), green(0.0), blue(0.0) {};

Color::Color(ifstream& ifs)
{
  ifs >> red >> green >> blue;
}

Color::Color(double r, double g, double b) : red(r), green(g), blue(b)
{}

Color Color::add(Color c) {
	return Color(red + c.red, green + c.green, blue + c.blue);
}

Color Color::scale(Color c) {
	return Color(red * c.red, green * c.green, blue * c.blue);
}

bool Color::isEqual(Color c) {
	return red == c.red && green == c.green && blue == c.blue;
}



Figure::Figure(){}

void Figure::initFigure(ifstream& ifs)
{
 ambient = Color(ifs);
 diffuse = Color(ifs);
 specular = Color(ifs);
 reflectivity = Color(ifs);
 transmissivity = Color(ifs);
 ifs >> shininess >> indexOfRefraction >> rFlag >> tFlag;
}

Color Figure::getTransmissivity() {
	return transmissivity;
}

Color Figure::getReflectivity() {
	return reflectivity;
}

Color Figure::getColorAmbient() {
	return ambient;
}

Light::Light(ifstream& ifs) : position(ifs), shading(ifs)
{
  ifs >> c0 >> c1 >> c2;
}

Vec Light::getPosition() {
	return position;
}


Sphere::Sphere(ifstream& ifs) : center(ifs)
{
  ifs >> radius;
  initFigure(ifs);
}

Plane::Plane(ifstream& ifs) : abcVector(ifs)
{
  ifs >> dScalar;
  initFigure(ifs);
  direction1 = Vec(ifs);
  direction2 = Vec(ifs);
}

Ray::Ray(Vec a, Vec b) {
	p0 = a;
	p1 = b;
}

Vec Ray::getP0() {
	return p0;
}

Vec Ray::getP1() {
	return p1;
}

void parseSceneFile(char* sceneName)
{
  double bgr, bgg, bgb;
  double ar, ag, ab;
  ifstream ifs;
  assert (sceneName != 0);
  ifs.open(sceneName);
  ifs >> cameraX;
  ifs >> cameraY;
  ifs >> cameraZ;
  ifs >> zCoor;
  ifs >> minX >> maxX;
  ifs >> minY >> maxY;
  ifs >> bgr >> bgg >> bgb;
  backgroundColor = Color(bgr,bgg,bgb);
  ifs >> ar >> ag >> ab;
  ambient = Color(ar,ag,ab);
  ifs >> maxDepth;
  ifs >> horizontalResolution >> verticalResolution;
  int numLights, numSpheres, numPlanes;
  ifs >> numLights;
  for (int i=0; i<numLights; ++i) lightList.push_front(new Light(ifs));
  ifs >> numSpheres;
  for (int i=0; i<numSpheres; ++i) shapeList.push_front(new Sphere(ifs));
  ifs >> numPlanes;
  for (int i=0; i<numPlanes; ++i) shapeList.push_front(new Plane(ifs));
  ifs.close();
}

void RT_algorithm()
{
	Vec origin = Vec();
	for (int u = 0; u < win_w; u++)
	{
		for (int v = 0; v < win_h; v++)
		{
			Vec p = Vec(u, v, 0); // ?? is z 0 here or 1 ??
			Ray R = Ray(origin, p);
			//Set color of (u, v) to = RT_trace(R, 1);
		}
	}
}


Color RT_trace(const Ray& r, double depth)
{
	double epsilon = .01;
	double maxT = 5;
	pair<double, Figure*> pair = nearestIntersection(r, epsilon, maxT);
	double i = pair.first;
	Figure* figure = pair.second;

	if (pair.second != nullptr) {
		
		Ray n =
			; // surface normal for figure at i
			// figure out if it's entering the object
		return RT_shade(figure, r, i, normal, entering, depth);
	}
	else {
		return backgroundColor;
	}
	
}

pair<double, Figure*> nearestIntersection(const Ray& r,
	double minT, double maxT,
	bool mayBeTransparent = true)
{
	for (int i = 0; i < maxT; i = i + minT) // Because minT is the epsilon? What is maxT / minT ??
	{
		for each(Figure* fig in shapeList)
		pair<double, Figure*> intersection = fig->intersection(r, i, maxT); //intersection(r, minT, maxT) // mad because we aren't calling intersection on anything
		if (intersection != NULL) {
			return intersection;
		}
	}

}

double Sphere::intersection(const Ray& r, double minT, double maxT) const
{

}

double Plane::intersection(const Ray& r, double minT, double maxT) const
{
	double t = (dScalar - ((r->getP0())->dot(abcVector))) / (((r->getP1())->subtract(r->getP0))->dot(abcVector);

	if (t == 0) {
		return NULL; // will check this higher up b/c NULL = 0;
	}
	else {
		return t;
	}
}

// Return the color of OBJ at I as seen along RAY including effects of ambient and directional light sources, and reflected and transmitted rays
Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth)
{
	Color newColor = obj->getColorAmbient(); // Ambient coor term for OBJ at I ??
	newColor = newColor.add(RT_lights(obj, ray, i, normal));
	if (depth < maxDepth) {
		newColor = newColor.add(RT_reflect(obj, ray, i, normal, depth));
		newColor = newColor.add(RT_transmit(obj, ray, i, normal, entering, depth));
	}
	return newColor;
}

// Return color of OBJ at INT as seen along Ray, including diffuse and specular effects of each diretional light source
Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{
	Color newColor = Color();
	for each (Light* light in lightList) {
		Ray l_ray = Ray(i, light->getPosition());
		Vec lightDirection = (((light->getPosition())->subtract(i)).normalize()); // Let lightDirection be the direction of the L_RAY --> do this by normalizing L_ray
		double dotProduct = lightDirection.dot(normal);
		if (dotProduct > 0) {
			// adding diffuse and specular terms for light striking 
			Color diffuseColor = diffuseShade(obj, light, dotProduct);
			Color specularColor = specularShade(obj, normal, light,
				lightDirection, dotProduct,
				ray);
			Color lightColors = diffuseColor.add(specularColor);
			// OBJ along L_Ray at i, scaled by the opacity of objects intersecting L_RAY between i and light
			Color scaledLightColors = lightColors.scale(obj->getTransmissivity());
			newColor = newColor.add(scaledLightColors);
		}
	}
	return newColor;
}

Color specularShade(Figure* obj, const Vec& normal,
	Light* light, const Vec& lightDirection, double dotProduct,
	const Ray& ray)
{
	...
}


Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
	...
}
.;

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth)
{
	if ((obj->getTransmissivity()).isEqual(Color())) {
		Ray transmittedRay = Ray(i, normal);
		Vec transmittedDirection = ((normal->subtract(i)).normalize());// transmittedRay normalized
		double dotProduct = transmittedDirection.dot(normal);
		// If no total internal reflection
		if (dotProduct < 0) {
			Color newColor = RT_trace(transmittedRay, depth + 1);
			newColor = newColor.scale(obj->getTransmissivity());
			return newColor;
		}
	}
	else {
		// return the zero color
		return Color();
	}
}


Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, double depth)
{
	if (((obj->getReflectivity()).isEqual(Color()))) {
		Ray reflectedRay = Ray(i, normal);
		Color newColor = RT_trace(reflectedRay, depth + 1);
		//Scale newColor according to the inter-object reflection coefficients kr = (krr, krg, krb) of object OBJ
		newColor = newColor.scale(obj->getReflectivity());
		return newColor;
	}
	else {
		// return the zero color Vector
		return Color();
	}
}

int main(int, char *argv[])
{
    parseSceneFile(argv[1]);
	initializeImage();
	RT_algorithm();
	writeImageFile();
	return 0;
}

