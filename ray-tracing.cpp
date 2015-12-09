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

Vec Vec::add(Vec v) {
	return Vec(x + v.x, y + v.y, z + v.z);
}

Vec Vec::scale(double t) {
	return Vec(x * t, y * t, z * t);
}


Color::Color(): red(0.0), green(0.0), blue(0.0) {};

Color::Color(ifstream& ifs)
{
  ifs >> red >> green >> blue;
}

Color::Color(double r, double g, double b) : red(r), green(g), blue(b)
{}

double Color::getRed() {
	return red;
}

double Color::getBlue() {
	return blue;
}

double Color::getGreen() {
	return green;
}

Color Color::add(Color c) {
	return Color(red + c.red, green + c.green, blue + c.blue);
}

Color Color::multiply(Color c) {
	return Color(red * c.red, green * c.green, blue * c.blue);
}

bool Color::isEqual(Color c) {
	return red == c.red && green == c.green && blue == c.blue;
}

string Color::toString() {
	int r = (int)round(red * 255);
	int g = (int)round(green * 255);
	int b = (int)round(blue * 255);
	char toReturn[12];
	sprintf_s(toReturn, "%d %d %d", r, g, b);
	return string(toReturn);

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

//Vec* Figure::getNormal(Vec* i) {
//	return NULL;
//}

//double Figure::intersection(const Ray& r, double minT, double maxT) const {
//	return NULL;
//}

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



Vec* Sphere::getNormal(Vec* i) {
	// normal = line from the center to the point
	return new Vec(center.subtract(*(i)));
}


double Sphere::intersection(const Ray& r, double minT, double maxT) const
{
	Ray ray = r;
	Vec h = (ray.getP1()).subtract(ray.getP0());
	Vec k = (ray.getP0()).subtract(ray.getP1());
	
	// Solve the equation.... using the quadratic formula 
	// t^2(h*h) + 2t(h*k) + k*k - r^2 = 0;
	// do - first : (-b - squareroot(b^2 - 4ac))/2a
	double t = (-(h.dot(k)) - sqrt((h.dot(k) * h.dot(k)) - (4 * h.dot(h) * k.dot(k)))) / (2 * k.dot(k));
	if (t < minT || t > maxT) {
		// then do + if the first invalid : (-b + squareroot(b^2 - 4ac))/2a
		t = (-(h.dot(k)) + sqrt((h.dot(k) * h.dot(k)) - (4 * h.dot(h) * k.dot(k)))) / (2 * k.dot(k));
		if (t < minT || t > maxT) {
			// then return null if none of them worked
			return NULL;
		}
	}
	return t;
}

Vec* Plane::getNormal(Vec* i) {
	return &abcVector;
}

double Plane::intersection(const Ray& r, double minT, double maxT) const
{
	Ray ray = r;
	double t = (dScalar - ((ray.getP0()).dot(abcVector))) / (((ray.getP1()).subtract(ray.getP0())).dot(abcVector));

	if (t == 0) {
		return NULL; // will check this higher up b/c NULL = 0;
	}
	else {
		return t;
	}
}
	
Color image[1000][1000];

//void initializeImage() {
//	for (int u = 0; u < horizontalResolution; u++) {
//		image[u] = new Color[verticalResolution];
//	}
//}

void writeImageFile() {
	ofstream myfile;
	myfile.open("picture.ppm");
	char* c = new char[horizontalResolution * verticalResolution * 12];
	int index = 0;
	myfile << "P3" << endl;
	//cout << "P3" << endl;
	myfile << "500 500" << endl;
	//cout << "500 500" << endl;
	myfile << "255" << endl;
	//cout << "255" << endl;
	for (int i = 0; i < horizontalResolution; i++) {
		for (int j = 0; j < verticalResolution; j++) {
			Color color = image[i][j];
			//string s = color.toString();
			int r = (int)round(color.getRed() * 255);
			int g = (int)round(color.getGreen() * 255);
			int b = (int)round(color.getBlue() * 255);
			myfile << r << " " << g << " " << b << endl;
			//cout << r << " " << g << " " << b << endl;
			//myfile.write(s, 12 + index);
			//for (int k = 0; k < 13; k++) {
			//	c[index + k] = s[k];
			//}
			//index = index + 12;
		}
	}
	//myfile.write(c, horizontalResolution * verticalResolution * 12);
	myfile.close();
}

Color specularShade(Figure* obj, const Vec& normal,
	Light* light, const Vec& lightDirection, double dotProduct,
	const Ray& ray)
{
	//...
	return Color();
}


Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
	///...
	return Color();
}

// Return color of OBJ at INT as seen along Ray, including diffuse and specular effects of each diretional light source
Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{
	Color newColor = Color();
	for each (Light* light in lightList) {
		Ray l_ray = Ray(i, light->getPosition());
		Vec lightDirection = (((light->getPosition()).subtract(i)).normalize()); // Let lightDirection be the direction of the L_RAY --> do this by normalizing L_ray
		double dotProduct = lightDirection.dot(normal);
		if (dotProduct > 0) {
			// adding diffuse and specular terms for light striking 
			Color diffuseColor = diffuseShade(obj, light, dotProduct);
			Color specularColor = specularShade(obj, normal, light,
				lightDirection, dotProduct,
				ray);
			Color lightColors = diffuseColor.add(specularColor);
			// OBJ along L_Ray at i, scaled by the opacity of objects intersecting L_RAY between i and light
			Color scaledLightColors = lightColors.multiply(obj->getTransmissivity());
			newColor = newColor.add(scaledLightColors);
		}
	}
	return newColor;
}

// Return the color of OBJ at I as seen along RAY including effects of ambient and directional light sources, and reflected and transmitted rays
Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth)
{
	Color newColor = (obj->getColorAmbient()).multiply(ambient); // Ambient coor term for OBJ at I ??
																 //newColor = newColor.add(RT_lights(obj, ray, i, normal));
																 //if (depth < maxDepth) {
																 //newColor = newColor.add(RT_reflect(obj, ray, i, normal, depth));
																 //newColor = newColor.add(RT_transmit(obj, ray, i, normal, entering, depth));
																 //}
	return newColor;
}

pair<double, Figure*> nearestIntersection(const Ray& r,
	double minT, double maxT,
	bool mayBeTransparent = true)
{

	double currentT = 1000000;
	Figure* currentFig = NULL;
	for each(Figure* fig in shapeList) {
		double intersection = fig->intersection(r, minT, maxT);
		if (intersection != NULL) {
			if (intersection < currentT && intersection > minT && intersection < maxT) {
				currentT = intersection;
				currentFig = fig;
			}
		}
	}

	return pair<double, Figure*>(currentT, currentFig);
}

// Return the color of the first object point hit by RAY
// Return the background color if RAY hits no object
Color RT_trace(const Ray& r, double depth)
{
	double epsilon = .01;
	double maxT = 5;
	pair<double, Figure*> pair = nearestIntersection(r, epsilon, maxT);
	double t = pair.first;
	Figure* figure = pair.second;
	Ray ray = r;

	if (pair.second != nullptr) {
		// Point of Intersection = p0 + t*(P1-P0)
		// add and scale for vector ->
		Vec pOI = (ray.getP0()).add((ray.getP1().subtract(ray.getP0())).scale(t));
		Vec* n = figure->getNormal(&pOI); // surface normal for figure at t
										  // figure out if it's entering the object ???
										  // if the ray is not traveling through an object, then you must be entering
										  // ray will have a pointer to the object that it is traveling through if there is one. 
										  // set in rt_refract
		return RT_shade(figure, r, pOI, *n, false, depth);
	}
	else {
		return backgroundColor;
	}

}

void RT_algorithm()
{
	Vec origin = Vec();
	for (int u = 0; u < horizontalResolution; u++)
	{
		for (int v = 0; v < verticalResolution; v++)
		{
			Vec p = Vec(u, v, 0);
			Ray R = Ray(origin, p);
			image[u][v] = RT_trace(R, 1);
		}
	}
}

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth)
{
	if ((obj->getTransmissivity()).isEqual(Color())) {
		Vec norm = normal;
		Ray transmittedRay = Ray(i, norm);
		Vec transmittedDirection = ((norm.subtract(i)).normalize());// transmittedRay normalized
		double dotProduct = transmittedDirection.dot(norm);
		// If no total internal reflection
		if (dotProduct < 0) {
			Color newColor = RT_trace(transmittedRay, depth + 1);
			newColor = newColor.multiply(obj->getTransmissivity());
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
		newColor = newColor.multiply(obj->getReflectivity());
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
	//initializeImage();
	RT_algorithm();
	writeImageFile();
	return 0;
}

