/*
@Bingchen Liu
@Joshua
@Micheal
*/

#include "obstacles/GJK_EPA.h"

#define THRESHOLD 0.0002

SteerLib::GJK_EPA::GJK_EPA()
{
}

//*********used for GJK below, written by bingchen**********************

// _shape is a set that contains all vectors (from origin to the control point) in a shape
void getShapeCenter(Util::Vector& c, const std::vector<Util::Vector>& _shape)
{
	c.x = 0;
	c.y = 0;
	c.z = 0;
	for (std::vector<Util::Vector>::const_iterator it = _shape.begin(); it != _shape.end(); ++it)
	{
		c.x += it->x;
		c.y += it->y;
		c.z += it->z;
	}
	if (!_shape.empty())
	{
		c.x /= _shape.size();
		c.y /= _shape.size();
		c.z /= _shape.size();
	}
	return;
}

// iterate each two edges in the shape, determine if it is convex based on the angles of each two edges
bool isConvex(const std::vector<Util::Vector>& _shape) {
	if (_shape.size() <= 2) {
		std::cerr << "no enough points in shape" << std::endl;
		return false;
	}

	for (int current = 0; current < _shape.size(); current++) {
		int next = (current + 1 == _shape.size()) ? 0 : (current + 1);
		int last = (current == 0) ? (_shape.size() - 1) : (current - 1);

		Util::Vector edge1 = _shape[next] - _shape[current];
		Util::Vector edge2 = _shape[current] - _shape[last];
		//cross product larger than 0 means the outer angle is less than 180, so it is not convex
		if (Util::cross(edge1, edge2).y > 0 || Util::cross(edge1, edge2).x > 0 || Util::cross(edge1, edge2).z > 0) {
			return false;
		}
	}
	return true;
}

bool detectConvex(const std::vector<Util::Vector>& _shape) {
	if (_shape.size() <= 2) {
		std::cerr << "not enough vertices in shape!!!" << std::endl;
		return false;
	}

	for (int i = 0; i < _shape.size(); ++i) {
		int next = i + 1;
		if (next == _shape.size()) {
			next = 0;
		}
		Util::Vector v1 = _shape[i];
		Util::Vector v2 = _shape[next];
		Util::Vector edge12 = v2 - v1;
		std::vector<Util::Vector> x_multiply_list;
		for (int j = 0; j < _shape.size(); ++j) {
			if (j == i || j == next) continue;
			Util::Vector v3 = _shape[j];
			Util::Vector edge13 = v3 - v1;
			Util::Vector x_multiply = Util::cross(edge12, edge13);
			x_multiply_list.push_back(x_multiply);
		}

		for (int j = 1; j < x_multiply_list.size(); ++j) {
			if (x_multiply_list[0] * x_multiply_list[j] < 0) {
				//			std::cout << "It's not a convex!!!" << std::endl;
				return false;
			}
		}
		x_multiply_list.clear();
	}
	return true;
}
// find the point in each shape follow a given direction, the return point used when calculate Minkowski Difference
Util::Vector getFarthestPointInDirection(const std::vector<Util::Vector>& _shape, Util::Vector direction)
{
	float distance = _shape[0].x * direction.x + _shape[0].y * direction.y + _shape[0].z * direction.z;
	int index = 0;
	for (int current = 0; current < _shape.size(); ++current) {
		float tmp_distance = _shape[current].x * direction.x + _shape[current].y * direction.y + _shape[current].z * direction.z;
		if (tmp_distance > distance) {
			distance = tmp_distance;
			index = current;
		}
	}
	return _shape[index];
}

// it will return a point in Minkowski Difference
Util::Vector support(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, Util::Vector direction)
{
	Util::Vector A_point = getFarthestPointInDirection(_shapeA, direction);

	direction *= -1;    // reverse the direction 
	Util::Vector B_point = getFarthestPointInDirection(_shapeB, direction);

	return (A_point - B_point); // shape A - shape B 
}

// read the last two points in the list, and get the direction that used to find the next point
// the direction is vertical to the line consists the last two points
Util::Vector getDirection(const std::vector<Util::Vector>& simplexList)
{
	if (simplexList.size() < 2)
	{
		std::cout << "Error: getting direction" << std::endl;
		return Util::Vector(0, 0, 0);
	}
	std::vector<Util::Vector>::const_reverse_iterator rit = simplexList.rbegin();
	Util::Vector lastSimplex(rit->x, rit->y, rit->z);
	rit++;
	Util::Vector secondLastSimplex(rit->x, rit->y, rit->z);

	Util::Vector AB = lastSimplex - secondLastSimplex;
	Util::Vector AO = lastSimplex * (-1);
	Util::Vector direction = AO*(AB*AB) - AB*(AB*AO); // (AB X AO) X AB, 
													  //Note that the following triple product expansion is used:
													  //(A x B) x C = B(C.dot(A)) �C A(C.dot(B)) to evaluate the triple product.;
													  //this can get the direction vector that perpendicular to AB and points to origin
	return direction;
}

// check when the simplex is a line (AB)
bool onSegment(Util::Vector A, Util::Vector B, Util::Vector Origin)
{
	if (A == Origin || B == Origin) {
		return true;
	}
	if (A.x == 0 && B.x == 0) {
		if ((A.z * B.z) < 0) {
			return true;
		}
	}
	else if (A.z == 0 && B.z == 0) {
		if ((A.x * B.x) < 0) {
			return true;
		}
	}
	else if (A.x != B.x) {
		float k = (A.z - B.z) / (A.x - B.x);
		if ((A.z - k*A.x) == 0) {
			return true;
		}
	}
	return false;
}

// determine if a simplex (triangle) covers the origin 
bool containsOrigin(std::vector<Util::Vector>& simplexList, Util::Vector Origin)
{
	if (simplexList.size() != 3)
	{
		std::cout << "Error: check containsOrigin" << std::endl;
		return false;
	}

	Util::Vector C = simplexList[0];
	Util::Vector B = simplexList[1]; // the second last added simplex
	Util::Vector A = simplexList[2]; // the last added simplex

	Util::Vector AB = B - A;
	Util::Vector AC = C - A;
	Util::Vector AO = A * (-1); // 0 - a
	Util::Vector prep_AB = AC*(AB*AB) - AB*(AB*AC); // (AB X AC) X AB, the vector that perpendicular to AB;  
	Util::Vector prep_AC = AB*(AC*AC) - AC*(AC*AB); // (AC X AB) X AC;  

	if (onSegment(A, B, Origin) || onSegment(A, C, Origin))
	{
		return true;
	}

	if (prep_AB * AO <= 0)
	{// not in the R4 area, remove C
		simplexList.erase(simplexList.begin());
		return false;
	}
	else if (prep_AC * AO <= 0)
	{// not in the R3 area, remove b
		simplexList.erase(simplexList.begin() + 1);
		return false;
	}
	return true;

}

//*********used for GJK above, written by bingchen**********************

//*********used for EPA below, written by bingchen**********************
bool getShortestEdge(const std::vector<Util::Vector>& simplexList, SteerLib::Edge& shortestEdge)
{
	float shortestDistance = 10000;

	for (int i = 0; i < simplexList.size(); ++i)
	{
		int next = i + 1;
		if (next == simplexList.size()) {
			next = 0;
		}

		Util::Vector AB = simplexList[next] - simplexList[i];
		Util::Vector AO = -1 * simplexList[i];
		Util::Vector perp;
		//check if A, B, O are on same line
		if (fabs(fabs(AB*AO) - AB.norm() * AO.norm()) < THRESHOLD) {
			perp.x = -1 * AB.z;
			perp.y = 0;
			perp.z = AB.x;
		}
		else {
			perp = AO*(AB*AB) - AB*(AB*AO);
		}

		//printVector("perp", perp);
		//printVector("AB", AB);
		perp = Util::normalize(perp);
		float distance = fabs(AO*perp);
		if (shortestDistance > distance) {
			shortestEdge.p1 = simplexList[i];
			shortestEdge.p2 = simplexList[next];
			shortestEdge.perp = -1 * perp;
			shortestEdge.index = next;
			shortestEdge.distance = distance;
			shortestDistance = distance;
		}
	}

	if (shortestDistance == FLT_MAX) {
		return false;
	}
	else {
		return true;
	}

}


bool EPA(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector> simplexList)
{
	while (1) {
		SteerLib::Edge shortestEdge;
		if (!getShortestEdge(simplexList, shortestEdge)) std::cerr << "No shortestEdge found!!!" << std::endl;

		Util::Vector simplex = support(_shapeA, _shapeB, shortestEdge.perp);
		float difference = fabs(simplex * shortestEdge.perp - shortestEdge.distance);
		if (difference < THRESHOLD) {
			return_penetration_depth = shortestEdge.distance;
			return_penetration_vector = shortestEdge.perp;
			return true;
		}
		else {
			simplexList.insert(simplexList.begin() + shortestEdge.index, simplex);
		}
	}
	return false;
}
//*********used for EPA above, written by bingchen**********************

//*********used for decomposition below, written by bingchen**********************

// used to detect how the programe compute the angle between two edges
// clockwise or not between two edges determines whether the angle larger than 180 or not
bool detectClockwise(const std::vector<Util::Vector>& _shape) {

	int rightMostIndex = -1;
	float rightMostX = _shape[0].x;
	for (int i = 0; i < _shape.size(); ++i) {
		if (_shape[i].x > rightMostX) {
			rightMostIndex = i;
		}
	}
	rightMostX = _shape[rightMostIndex].x;

	int last = (rightMostIndex - 1) % _shape.size();
	int next = (rightMostIndex + 1) % _shape.size();

	// the angle that on the right most point of a shape must be less than 180(PI)
	Util::Vector X_last = _shape[rightMostIndex] - _shape[last];
	Util::Vector X_next = _shape[next] - _shape[rightMostIndex];
	Util::Vector cross_product = Util::cross(X_last, X_next);

	if (cross_product.y < 0) {
		return true;
	}
	else if (cross_product.y > 0) {
		return false;
	}
	else {
		std::cerr << "detect clockwise failed!!" << std::endl;
		return true;
	}
}

//get all the points in a concave angle of the shape
std::vector<int>  getConcavePoints(const std::vector<Util::Vector>& _shape, bool clockwise) {
	std::vector<int> concavePoints;
	for (int current = 0; current < _shape.size(); ++current) {

		int last = (current - 1) % _shape.size();
		int next = (current + 1) % _shape.size();
		Util::Vector X_last = _shape[current] - _shape[last];
		Util::Vector X_next = _shape[next] - _shape[current];
		Util::Vector cross_product = Util::cross(X_last, X_next);
		if (clockwise) {
			if (cross_product.y > 0) {
				concavePoints.push_back(current);
			}
		}
		else {
			if (cross_product.y < 0) {
				concavePoints.push_back(current);
			}
		}
	}
	return concavePoints;
}


//***** this is specifically for polygons2 test below , it is wrong so not used *****
bool isInOneLine(const Util::Vector& point1, const Util::Vector& point2, const Util::Vector& point3) {

	float x1 = point1.x - point2.x;
	float z1 = point1.z - point2.z;

	float x2 = point1.x - point3.x;
	float z2 = point1.z - point3.z;

	if ((x1 * z2) == (x2 * z1)) {
		/*
		std::cerr << "point1: " << " x: " << point1.x << " y: " << point1.y << " z: " << point1.z << std::endl;
		std::cerr << "point2: " << " x: " << point2.x << " y: " << point2.y << " z: " << point2.z << std::endl;
		std::cerr << "point3: " << " x: " << point3.x << " y: " << point3.y << " z: " << point3.z << std::endl;
		*/
		return true;
	}

	return false;
}

bool getIntersectPoint(const Util::Vector& point1, const Util::Vector& point2, const Util::Vector& point3, const Util::Vector& point4, Util::Vector& intersectPoint) {
	// y1 = k1 * x + b1  point1, point2     k1 = (point1.z - point2.z) / ( point1.x - point2.x ) 
	// y2 = k2 * x + b2
	float k1 = (point1.z - point2.z) / (point1.x - point2.x);
	float b1 = point1.z - k1 * point1.x;

	float k2 = (point3.z - point4.z) / (point3.x - point4.x);
	float b2 = point3.z - k1 * point3.x;

	float intersectX = (b2 - b1) / (k1 - k2);
	float intersectZ = intersectX * k1 + b1;

	// is parallel, then no intersect point
	if (k1 == k2) {
		return false;
	}
	/*
	if (point1.x == point2.x) {
	if (point3.x != point4.x ) {
	intersectX = point1.x;
	intersectZ = k2 * intersectX + b2;
	}
	}
	else if (point3.x == point4.x) {
	if (point1.x != point2.x) {
	intersectX = point3.x;
	intersectZ = k1 * intersectX + b1;
	}
	}*/

	//std::cerr <<"intersect x: " << intersectX << " intersect z: " << intersectZ << std::endl;

	Util::Vector newPoint{ intersectX, 0, intersectZ };

	if (intersectX <= std::max(point1.x, point2.x) && intersectX >= std::min(point1.x, point2.x)) {
		intersectPoint = newPoint;
		return true;
	}

	return false;
}
//***** this is specifically for polygons2 test above *****

//removing the concave angle by using the point in that angle and its last two points to 
//create a new triangle, then do same operation on the shape built by the rest points recursively 
std::vector< std::vector<Util::Vector> > decompositeShape(const std::vector<Util::Vector>& _shape) {
	std::vector< std::vector<Util::Vector>> decompositedShapeList;
	std::vector< Util::Vector> tmpShape = _shape;

	while (!detectConvex(tmpShape)) {

		std::vector<int> concavePointIndexs = getConcavePoints(tmpShape, detectClockwise(tmpShape));
		std::vector<int> removedPoints;

		// create new triangles in each point which has a concave angle
		for (int i = 0; i < concavePointIndexs.size(); i++) {
			int last_1 = (concavePointIndexs[i] == 0) ? tmpShape.size() - 1 : concavePointIndexs[i] - 1;
			int last_2 = (last_1 == 0) ? tmpShape.size() - 1 : last_1 - 1;
			std::vector< Util::Vector> tmp_triangle{ tmpShape[last_2], tmpShape[last_1], tmpShape[concavePointIndexs[i]] };

			decompositedShapeList.push_back(tmp_triangle);

			removedPoints.push_back(last_1);
		}

		//create new shape with rest points, exclude the points that used to create triangle,
		// only one piont in each triangle is emoved 
		std::vector< Util::Vector> newTmpShape;
		for (int i = 0; i < tmpShape.size(); i++) {
			if (std::find(removedPoints.begin(), removedPoints.end(), i) == removedPoints.end()) {
				newTmpShape.push_back(tmpShape[i]);
			}
		}
		tmpShape = newTmpShape;

		if (tmpShape.size() <= 3) {
			break;

		}
	}
	decompositedShapeList.push_back(tmpShape);
	return decompositedShapeList;

}

//*********used for decomposition above, written by bingchen,**********************

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	//***** decomposition part below, written by bingchen, called when shape A or B is not a convex ****
	/*
	if (!isConvex(_shapeA)) {
	std::cerr << "A is not convex" << std::endl;
	}
	if (!isConvex(_shapeB)) {
	std::cerr << "B is not convex" << std::endl;
	}
	*/
	if (!detectConvex(_shapeA) || !detectConvex(_shapeB)) {
		//std::cerr << "shape A or B is not a convex!!" << std::endl;
		std::vector<std::vector<Util::Vector>> decomp_A = decompositeShape(_shapeA);
		std::vector<std::vector<Util::Vector>> decomp_B = decompositeShape(_shapeB);

		for (int i = 0; i < decomp_A.size(); ++i) {
			for (int j = 0; j < decomp_B.size(); ++j) {
				if (intersect(return_penetration_depth, return_penetration_vector, decomp_A[i], decomp_B[j])) {
					return true;
				}
			}
		}
		return false;
	}
	//***** decomposition part above, written by bingchen, called when shape A or B is not a convex ****

	Util::Vector Origin(0, 0, 0);
	std::vector<Util::Vector> simplexList;
	// compute the shape Center
	Util::Vector centerShapeA(0, 0, 0);
	Util::Vector centerShapeB(0, 0, 0);
	getShapeCenter(centerShapeA, _shapeA);
	getShapeCenter(centerShapeB, _shapeB);

	Util::Vector direction = centerShapeA - centerShapeB;
	if (direction == Origin)
	{
		return true;
	}
	Util::Vector simplex = support(_shapeA, _shapeB, direction);
	simplexList.push_back(simplex);

	direction *= -1; // reverse direction
	simplex = support(_shapeA, _shapeB, direction);
	simplexList.push_back(simplex);        // get the first two simplex. build a line

	int index = 0;
	while (true)
	{
		index++;
		direction = getDirection(simplexList);

		if (direction == Origin) // if direction = 0,0,0 // origin is on the line AB
		{
			//std::cout << " direction == Origin" << std::endl;
			if (simplexList.size() != 2) {
				std::cout << "Error: Direction Origin " << std::endl;
				break;
			}

			std::vector<Util::Vector>::iterator it = simplexList.begin();
			Util::Vector p1(it->x, it->y, it->z);
			it++;
			Util::Vector p2(it->x, it->y, it->z);
			if (onSegment(p1, p2, Origin))
			{
				EPA(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplexList);
				return true;
			}
			else
				return false;
		}

		simplex = support(_shapeA, _shapeB, direction);
		simplexList.push_back(simplex);

		if (simplex * direction <= 0) {
			//the new simplex is not past the origin in the direction
			//impossibly contain the origin 
			return false;
		}
		else {
			if (containsOrigin(simplexList, Origin))
			{
				EPA(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplexList);
				return true;
			}
		}
	}

	return false; // There is no collision
}