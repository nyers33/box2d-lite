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

#include "box2d-lite/Body.h"

Body::Body()
{
	position.Set(0.0f, 0.0f);
	rotation = 0.0f;
	velocity.Set(0.0f, 0.0f);
	angularVelocity = 0.0f;
	force.Set(0.0f, 0.0f);
	torque = 0.0f;
	friction = 0.2f;

	dim.Set(1.0f, 0.5f);
	mass = FLT_MAX;
	invMass = 0.0f;
	I = FLT_MAX;
	invI = 0.0f;
}

void Body::Set(const Vec2& d, float m)
{
	position.Set(0.0f, 0.0f);
	rotation = 0.0f;
	velocity.Set(0.0f, 0.0f);
	angularVelocity = 0.0f;
	force.Set(0.0f, 0.0f);
	torque = 0.0f;
	friction = 0.2f;

	dim = d;
	mass = m;

	float h = GetH();
	float r = GetR();
	float rho = mass / (2.0f * h * r + r * r * k_pi);
	
	if (mass < FLT_MAX)
	{
		invMass = 1.0f / mass;
		I = rho * r * (12.0f * powf(h, 3.0f) + 3.0f * powf(h, 2.0f) * r * k_pi + 64.0f * h * powf(r, 2.0f) + 6.0f * powf(r, 3.0f) * k_pi) / 12.0f;
		invI = 1.0f / I;
	}
	else
	{
		invMass = 0.0f;
		I = FLT_MAX;
		invI = 0.0f;
	}
}
