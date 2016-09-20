//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"
#include <cmath>
using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	if (!checkRobust()) return;

	// ***********written by Bingchen*****************

	// Draw the curve follow the time line, from time = 0 to time = end time, with step size "window"
	float time = 0;
	float endTime = controlPoints[controlPoints.size() - 1].time;

	Point startPoint = controlPoints[0].position;
	Point nextPoint;

	float timeWindow = 0.5;
	while (time <= endTime) {
		if (calculatePoint(nextPoint, time))
		{
			DrawLib::drawLine(startPoint, nextPoint, curveColor, curveThickness);
			startPoint = nextPoint;
		}
		else
		{
			std::cout << "Sorry I cannot find next point at time: " << time << std::endl;
		}
		time += window;
	}
	// ***********written by Bingchen*****************

	// Robustness: make sure there is at least two control point: start and end points
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve

	return;
#endif
}

// ***********written by Bingchen*****************
//define a structure that is used to sort the cp vector in ascending order
bool compareTimeAscending(CurvePoint cp1, CurvePoint cp2) {
	return (cp1.time < cp2.time);
}

bool compareSameTime(CurvePoint cp1, CurvePoint cp2) {
	return (cp1.time == cp2.time);
}
// ***********written by Bingchen*****************


// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	// ***********written by Bingchen*****************

	std::sort(controlPoints.begin(), controlPoints.end(), compareTimeAscending);
	controlPoints.erase(std::unique(controlPoints.begin(), controlPoints.end(), compareSameTime), controlPoints.end());
	// ***********written by Bingchen*****************

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	// ***********written by Bingchen*****************

	// check how many points do we have in the array(vector),
	//if less than 2 then we cannot draw curve, so return false
	return !(controlPoints.size() < 2);

	// ***********written by Bingchen*****************

	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	// ***********written by Bingchen*****************

	// It is used to get the index of next controlPoint in the vextor of controlPoints
	// by the given time.
	if (controlPoints[0].time == time) {
		nextPoint = 0;
		return true;
	}
	for (int i = 0; i < controlPoints.size() - 1; ++i) {
		if (controlPoints[i].time < time && controlPoints[i + 1].time >= time) {
			nextPoint = i + 1;
			return true;
		}
	}

	return false;
	// ***********written by Bingchen*****************

}


// ***********written by Bingchen*****************
// preparation for the 4 Fh(t) = T*Mh, just refer to the formula on slides
float get_h1(float t) {
	return (2 * t * t * t - 3 * t * t + 1);
}

float get_h2(float t) {
	return (-2 * t * t * t + 3 * t * t);
}

float get_h3(float t) {
	return (t * t * t - 2 * t * t + t);
}

float get_h4(float t) {
	return (t * t * t - t * t);
}
// ***********written by Bingchen*****************



// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	// ***********written by Bingchen*****************
	if (nextPoint == 0) {
		newPosition = controlPoints[nextPoint].position;
		return newPosition;
	}
	// Calculate time interval, and normal time required for later curve calculations
	// interval time is seen as 1 unit of time , which is t0 ~ t1
	// so then we have to normalize every single "time" into this interval
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;

	// Calculate position at t = time on Hermite curve
	Point p0 = controlPoints[nextPoint - 1].position;
	Point p1 = controlPoints[nextPoint].position;
	Vector r0 = controlPoints[nextPoint - 1].tangent;
	Vector r1 = controlPoints[nextPoint].tangent;

	float h1 = get_h1(normalTime);
	float h2 = get_h2(normalTime);
	float h3 = get_h3(normalTime);
	float h4 = get_h4(normalTime);
	//	float _dist = dist(p0, p1);

	// Q = T * M * G , F = T * M , F is _h, so we simple use all _hs time G
	// G  = [P0 , P1, R0, R1]
	newPosition.x = p0.x * h1 + p1.x * h2 + r0.x * h3 * intervalTime + r1.x * h4 * intervalTime;
	newPosition.y = p0.y * h1 + p1.y * h2 + r0.y * h3 * intervalTime + r1.y * h4 * intervalTime;
	newPosition.z = p0.z * h1 + p1.z * h2 + r0.z * h3 * intervalTime + r1.z * h4 * intervalTime;

	// ***********written by Bingchen*****************

	// Calculate position at t = time on Hermite curve

	// Return result
	return newPosition;
}



// ***********written by Bingchen **************
// we have to calculate both tangents and boundaries
// 1. first we calculate the tangent

float getTangent(float t0, float t1, float t2, float p0, float p1, float p2) {
	return ((t1 - t0) / (t2 - t0)) * ((p2 - p1) / (t2 - t1)) + ((t2 - t1) / (t2 - t0)) * ((p1 - p0) / (t1 - t0));
}

Vector getPointTangent(float t0, float t1, float t2, Point p0, Point p1, Point p2) {
	Vector s;
	s.x = getTangent(t0, t1, t2, p0.x, p1.x, p2.x);
	s.y = getTangent(t0, t1, t2, p0.y, p1.y, p2.y);
	s.z = getTangent(t0, t1, t2, p0.z, p1.z, p2.z);
	return s;
}

// 2. then we calculate the boundaries
float getBoundary(float t0, float t1, float t2, float p0, float p1, float p2) {
	return ((t2 - t0) / (t2 - t1)) * ((p1 - p0) / (t1 - t0)) - ((t1 - t0) / (t2 - t1)) * ((p2 - p0) / (t2 - t0));
}

Vector getPointBoundary(float t0, float t1, float t2, Point p0, Point p1, Point p2) {
	Vector s;
	s.x = getBoundary(t0, t1, t2, p0.x, p1.x, p2.x);
	s.y = getBoundary(t0, t1, t2, p0.y, p1.y, p2.y);
	s.z = getBoundary(t0, t1, t2, p0.z, p1.z, p2.z);
	return s;
}

Point getNewPosition(float normalTime, float intervalTime, Point p0, Point p1, Vector r0, Vector r1) {
	Point newPoint;

	float h1 = get_h1(normalTime);
	float h2 = get_h2(normalTime);
	float h3 = get_h3(normalTime);
	float h4 = get_h4(normalTime);
	//	float _dist = dist(p0, p1);

	// same as Hermit curve here
	newPoint.x = p0.x * h1 + p1.x * h2 + r0.x * h3 * intervalTime + r1.x * h4 * intervalTime;
	newPoint.y = p0.y * h1 + p1.y * h2 + r0.y * h3 * intervalTime + r1.y * h4 * intervalTime;
	newPoint.z = p0.z * h1 + p1.z * h2 + r0.z * h3 * intervalTime + r1.z * h4 * intervalTime;


	return newPoint;
}

// create a new non-exist point when there at no [0] and [end]
Point reversePoint(Point positive) {
	Point result;
	result.x = -positive.x;
	result.y = -positive.y;
	result.z = -positive.z;

	return result;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;


	int points_size = controlPoints.size();
	if (nextPoint == 0)
		return controlPoints[0].position;

	// Calculate time interval, and normal time required for later curve calculations
	float intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	float normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;

	// switch between 3 cases, 1: start point, where we don't have 4 constrain points,
	// 2: normal point, 3: last points, where we also don't have 4 points.
	if (nextPoint == 1) {
		Point p0 = controlPoints[nextPoint - 1].position;
		Point p1 = controlPoints[nextPoint].position;
		Point p2 = controlPoints[nextPoint + 1].position;
		//Point p3 = controlPoints[nextPoint + 2].position;

		float t0 = controlPoints[nextPoint -1 ].time;
		float t1 = controlPoints[nextPoint ].time;
		float t2 = controlPoints[nextPoint + 1].time;
		//float t3 = controlPoints[nextPoint + 2].time;

		Vector s0 = getPointTangent(t0, t1, t2, p0, p1, p2);
		Vector s1 = getPointBoundary(t0, t1, t2, p0, p1, p2);
		//		std::cout << "s0 = " << s0 << " s1 = " << s1 << std::endl;
		newPosition = getNewPosition(normalTime, intervalTime, p0, p1, s0, s1);

	}
	else if (nextPoint == (points_size - 1)) {
		Point p0 = controlPoints[nextPoint - 2].position;
		Point p1 = controlPoints[nextPoint - 1].position;
		Point p2 = controlPoints[nextPoint].position;
		float t0 = controlPoints[nextPoint - 2].time;
		float t1 = controlPoints[nextPoint - 1].time;
		float t2 = controlPoints[nextPoint].time;
		Vector s1 = getPointTangent(t0, t1, t2, p0, p1, p2);
		Vector s2 = getPointBoundary(t2, t1, t0, p2, p1, p0);
		//		std::cout << "s_n-1 = " << s1 << " s_n = " << s2 << std::endl;
		newPosition = getNewPosition(normalTime, intervalTime, p1, p2, s1, s2);

	}
	else {
		Point p0 = controlPoints[nextPoint - 2].position;
		Point p1 = controlPoints[nextPoint - 1].position;
		Point p2 = controlPoints[nextPoint].position;
		Point p3 = controlPoints[nextPoint + 1].position;
		float t0 = controlPoints[nextPoint - 2].time;
		float t1 = controlPoints[nextPoint - 1].time;
		float t2 = controlPoints[nextPoint].time;
		float t3 = controlPoints[nextPoint + 1].time;
		Vector s1 = getPointTangent(t0, t1, t2, p0, p1, p2);
		Vector s2 = getPointBoundary(t1, t2, t3, p1, p2, p3);
		//		std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;
		newPosition = getNewPosition(normalTime, intervalTime, p1, p2, s1, s2);
	}

	// ***********written by Bingchen***************** 

	// Return result
	return newPosition;
}