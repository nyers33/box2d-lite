/*
* Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
*
* Permission to use, copy, modify, distribute and sell this software
* and its documentation for any purpose is hereby granted without fee,
* provided that the above copyright notice appear in all copies.
* Erin Catto makes no representations about the suitability 
* of this software for any purpose.  
* It is provided "as is" without express or implied warranty.
*/

#ifndef BODY_H
#define BODY_H

#include "MathUtils.h"

struct Body
{
	Body();
	void Set(const Vec2& dim, float m);
	float GetH() { return dim.x; }
	float GetR() { return dim.y; }

	void AddForce(const Vec2& f)
	{
		force += f;
	}

	Vec2 position;
	float rotation;

	Vec2 velocity;
	float angularVelocity;

	Vec2 force;
	float torque;

	Vec2 dim; // dimensions (height & radius)

	float friction;
	float mass, invMass;
	float I, invI;
};

#endif
