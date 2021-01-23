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

#include "box2d-lite/Arbiter.h"
#include "box2d-lite/Body.h"

#include <limits>
#include <cmath>
#include <cfloat>

// Capsule feature numbering:
//
//         ^ y
//         |
//         1
//      /-----\
//     /       \
//    /         \
//   |           |
//   |           |
//   |     0     | --> x
//   |           |
//   |           |
//    \         /
//     \       /
//      \-----/
//         -1

enum FeatureNumbers
{
	NEG_Y = -1,
	BODY,
	POS_Y
};

void qrDecompGS(Mat22 mat, Mat22 matOut[2])
{
	Vec2 a1 = mat.col1;
	Vec2 a2 = mat.col2;

	Vec2 v1 = a1;
	float v1Norm = v1.Length();
	Vec2 e1 = (v1Norm > 1e-4) ? Normalize(v1) : Vec2(0.0f, 0.0f);
	Vec2 v2 = a2 - Dot(a2, e1) * e1;
	float v2Norm = v2.Length();
	Vec2 e2 = (v2Norm > 1e-4) ? Normalize(v2) : Vec2(0.0f,0.0f);

	matOut[0].col1 = e1;
	matOut[0].col2 = e2;

	matOut[1].col1.x = Dot(a1, e1);
	matOut[1].col1.y = 0.0f;
	matOut[1].col2.x = Dot(a2, e1);
	matOut[1].col2.y = Dot(a2, e2);
}

float sqDistPointSeg(Vec2 a, Vec2 b, Vec2 c)
{
	Vec2 ab = b - a;
	Vec2 ac = c - a;
	Vec2 bc = c - b;
	float e = Dot(ac, ab);
	float f = Dot(ab, ab);
	if (e < 0.0f)
	{
		return Dot(ac, ac);
	}
	else if ( e >= f )
	{
		return Dot(bc, bc);
	}
	else
	{
		return Dot(ac, ac) - e * e / f;
	}
}

Vec2 closestPtPointSeg(Vec2 a, Vec2 b, Vec2 c)
{
	Vec2 ab = b - a;
	float t = Dot(c - a, ab) / Dot(ab, ab);
	t = Clamp(t, 0.0f, 1.0f);
	return a + t * ab;

}

Vec2 parallelogramContainsPt(Vec2 o, Vec2 a, Vec2 b, Vec2 pt)
{
	Vec2 p = pt - o;
	float denomAbs = Abs(a.x * b.y - a.y * b.x);

	float mu, lambda;

	mu = (p.x * b.y - p.y * b.x) / (a.x * b.y - a.y * b.x);
	lambda = (p.x * a.y - p.y * a.x) / (a.y * b.x - a.x * b.y);

	if ( denomAbs > 1e-4 && (0.0f <= mu && mu <= 1.0f && 0.0f <= lambda && lambda <= 1.0f) )
	{
		// inside
		return pt;
	}
	else
	{
		// outside
		Vec2 edgeList[4][2] = { { Vec2(0.0f, 0.0f), Vec2(a.x, a.y) }, { Vec2(a.x, a.y), Vec2(a.x, a.y) + Vec2(b.x, b.y) }, { Vec2(a.x, a.y) + Vec2(b.x, b.y), Vec2(b.x, b.y) }, { Vec2(b.x, b.y), Vec2(0.0f, 0.0f) } };
		float minDistToEdge = std::numeric_limits<float>::max();
		int iEdge = -1;
		for (int i = 0; i < 4; i++)
		{
			float distToEdge = sqDistPointSeg(edgeList[i][0], edgeList[i][1], p);
			if (distToEdge <= minDistToEdge)
			{
				minDistToEdge = distToEdge;
				iEdge = i;
			}
		}

		if (iEdge == -1)
		{
			int alma = 666;
		}
		return closestPtPointSeg(edgeList[iEdge][0], edgeList[iEdge][1], p) + o;
	}

}

// The normal points from A to B
int Collide(Contact* contacts, Body* bodyA, Body* bodyB)
{
	Vec2 posA = bodyA->position;
	Vec2 posB = bodyB->position;
	Mat22 RotA(bodyA->rotation), RotB(bodyB->rotation);

	Vec2 p1 = posA - 0.5f * bodyA->GetH() * RotA.col2; // seg0p0
	Vec2 u1 = posA + 0.5f * bodyA->GetH() * RotA.col2; // seg0p1

	Vec2 p2 = posB - 0.5f * bodyB->GetH() * RotB.col2; // seg1p0
	Vec2 u2 = posB + 0.5f * bodyB->GetH() * RotB.col2; // seg1p1

	Vec2 s1 = u1 - p1;
	Vec2 s2 = u2 - p2;

	float r1 = bodyA->GetR();
	float r2 = bodyB->GetR();

	Mat22 qrMat[2];
	qrDecompGS(Mat22(s2, -s1), qrMat);

	Vec2 parallelogramOrigin = qrMat[1] * Vec2(0.0f, 0.0f) + qrMat[0].Transpose() * (p2 - p1);
	Vec2 parallelogramSideA = qrMat[1] * Vec2(1.0f, 0.0f) + qrMat[0].Transpose() * (p2 - p1) - parallelogramOrigin;
	Vec2 parallelogramSideB = qrMat[1] * Vec2(0.0f, 1.0f) + qrMat[0].Transpose() * (p2 - p1) - parallelogramOrigin;

	Vec2 opt = parallelogramContainsPt(parallelogramOrigin, parallelogramSideA, parallelogramSideB, Vec2(0.0f, 0.0f));
	float sqDist = Dot(opt, opt) + Dot((p2 - p1), (p2 - p1)) - Dot((p2 - p1), qrMat[0] * qrMat[0].Transpose() * (p2 - p1));

	if (sqDist < 0.0f) sqDist = 0.0f;

	if ((bodyA->GetR() + bodyB->GetR()) * (bodyA->GetR() + bodyB->GetR()) < sqDist)
	{
		return 0;
	}

	Vec2 rx = opt - qrMat[0].Transpose() * (p2 - p1);

	float parallelTest = Abs(Abs(Dot(s1, s2)) - s1.Length() * s2.Length());

	// line vs line nearest in point - not parallel
	if (Abs(qrMat[1].col1.x) > 1e-4 && Abs(qrMat[1].col2.y) > 1e-4 && parallelTest > 1e-4)
	{
		float x1 = rx.y / qrMat[1].col2.y;
		float x2 = (qrMat[1].col2.y * rx.x - qrMat[1].col2.x * rx.y) / (qrMat[1].col1.x * qrMat[1].col2.y);

		if (x2 < 0.0f)
			x2 = 0.0f;

		FeaturePair fp;
		if (x1 <= 0.0f + 1e-4)
		{
			fp.e.inFeature = NEG_Y;
		}
		else if (x1 >= 1.0f - 1e-4)
		{
			fp.e.inFeature = POS_Y;
		}
		else
		{
			fp.e.inFeature = BODY;
		}
		if (x2 <= 0.0f + 1e-4)
		{
			fp.e.outFeature = NEG_Y;
		}
		else if (x2 >= 1.0f - 1e-4)
		{
			fp.e.outFeature= POS_Y;
		}
		else
		{
			fp.e.outFeature = BODY;
		}
		Vec2 nContact = Normalize((p2 + x2 * s2) - (p1 + x1 * s1));
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = 0.5f * (((p1 + x1 * s1) + r1 * nContact) + ((p2 + x2 * s2) - r2 * nContact));
		contacts[0].feature = fp;
		return 1;
	}

	float s1LenSqd = Dot(s1, s1);
	float s2LenSqd = Dot(s2, s2);
	
	// point vs point
	if (s1LenSqd < 1e-4 && s2LenSqd < 1e-4)
	{
		FeaturePair fp;
		fp.e.inFeature = BODY;
		fp.e.outFeature = BODY;
		Vec2 nContact = Normalize(p2 - p1);
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = 0.5f * ((p1 + r1 * nContact) + (p2 - r2 * nContact));
		contacts[0].feature = fp;
		return 1;
	}

	// line vs line - parallel lines
	if (s1LenSqd * s2LenSqd > 1e-4)
	{
		float x1[2];
		float x2[2];
		FeaturePair fp;
		Vec2 nContact;

		if (s1LenSqd > s2LenSqd)
		{
			x1[0] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, p2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[1] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, u2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);

			x2[0] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, p1 - p2 + x1[0] * s1) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, p1 - p2 + x1[1] * s1) / sqrtf(s2LenSqd), 0.0f), 1.0f);

		}
		else
		{
			x2[0] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, p1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, u1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);

			x1[0] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, p2 - p1 + x2[0] * s2) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[1] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, p2 - p1 + x2[1] * s2) / sqrtf(s1LenSqd), 0.0f), 1.0f);
		}

		nContact = (p2 + 0.5f * (x2[0] + x2[1]) * s2) - (p1 + 0.5f * (x1[0] + x1[1]) * s1);
		float nContactNorm = nContact.Length();
		if (nContactNorm > 1e-4)
			nContact *= (1.0f / nContactNorm);
		else
			nContact = Vec2(0.0f, 1.0f);
		for (int i = 0; i < 2; ++i)
		{
			contacts[i].separation = sqrtf(sqDist) - r1 - r2;
			contacts[i].normal = nContact;
			contacts[i].position = (p1 + x1[i] * s1) + 0.5f * (r1 - r2 + sqrtf(sqDist)) * nContact;
			fp.e.inFeature = BODY;
			fp.e.outFeature = BODY;
			contacts[i].feature = fp;
		}
		if (((x1[0] <= 0.0f + 1e-4 && x1[1] <= 0.0f + 1e-4) || (x1[0] >= 1.0f - 1e-4 && x1[1] >= 1.0f - 1e-4)) && ((x2[0] <= 0.0f + 1e-4 && x2[1] <= 0.0f + 1e-4) || (x2[0] >= 1.0f - 1e-4 && x2[1] >= 1.0f - 1e-4)))
		{
			fp.e.inFeature = BODY;
			fp.e.outFeature = BODY;
			contacts[0].feature = fp;
			return 1;
		}
		else
			return 2;
	}

	// line vs point
	{
		float x1[2];
		float x2[2];
		FeaturePair fp;
		Vec2 nContact;

		if (s1LenSqd > s2LenSqd)
		{
			x1[0] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, p2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[1] = Min(Max(Dot((1.0f / sqrtf(s1LenSqd)) * s1, u2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);

			x2[0] = 0.0f;
			x2[1] = 1.0f;

			if (x1[0] == 0.0f)
			{
				fp.e.inFeature = NEG_Y;
			}
			else if (x1[0] == 1.0f)
			{
				fp.e.inFeature = POS_Y;
			}
			else
			{
				fp.e.inFeature = BODY;
			}
			fp.e.outFeature = BODY;
		}
		else
		{
			x1[0] = 0.0f;
			x1[1] = 1.0f;

			x2[0] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, p1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = Min(Max(Dot((1.0f / sqrtf(s2LenSqd)) * s2, u1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);

			if (x2[0] == 0.0f)
			{
				fp.e.outFeature = NEG_Y;
			}
			else if (x2[0] == 1.0f)
			{
				fp.e.outFeature = POS_Y;
			}
			else
			{
				fp.e.outFeature = BODY;
			}
			fp.e.inFeature = BODY;
		}

		nContact = Normalize((p2 + Vec2(x2[0] * s2.x, x2[1] * s2.y)) - (p1 + Vec2(x1[0] * s1.x, x1[1] * s1.y)));
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = (p1 + Vec2(x1[0] * s1.x, x1[1] * s1.y)) + 0.5f * (r1 - r2 + sqrtf(sqDist)) * nContact;
		contacts[0].feature = fp;
		return 1;
	}
}
