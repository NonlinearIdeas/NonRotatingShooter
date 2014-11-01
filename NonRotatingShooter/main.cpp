/********************************************************************
 * File   : main.cpp
 * Project: NonRotatingShooter
 *
 ********************************************************************
 * Created on 11/1/14 By Nonlinear Ideas Inc.
 * Copyright (c) 2014 Nonlinear Ideas Inc. All rights reserved.
 ********************************************************************
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any
 * damages arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any
 * purpose, including commercial applications, and to alter it and
 * redistribute it freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must
 *    not claim that you wrote the original software. If you use this
 *    software in a product, an acknowledgment in the product
 *    documentation would be appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and
 *    must not be misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source
 *    distribution.
 */

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <ostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <float.h>

using namespace std;

typedef signed char int8;
typedef signed short int16;
typedef signed int int32;
typedef signed long long int64;

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

typedef double float64;

#define Min(a,b) (a<b?a:b)
#define Max(a,b) (a<b?b:a)

// This 2D Vector is based on code that is in the
// Box2D physics engine.  It has been copied here
// so this example project is not dependent on the
// entire Box2D code base.

/// A 2D column vector.
struct Vec2
{
   /// Default constructor does nothing (for performance).
   Vec2() :
   x(0.0),
   y(0.0)
   {
   }
   
   /// Construct using coordinates.
   Vec2(float64 x_, float64 y_) : x(x_), y(y_) {}
   
   /// Set this vector to all zeros.
   void SetZero() { x = 0.0; y = 0.0; }
   
   /// Set this vector to some specified coordinates.
   void Set(float64 x_, float64 y_) { x = x_; y = y_; }
   
   /// Negate this vector.
   Vec2 operator -() const { Vec2 v; v.Set(-x, -y); return v; }
   
   /// Read from and indexed element.
   float64 operator () (int32 i) const
   {
      return (&x)[i];
   }
   
   /// Write to an indexed element.
   float64& operator () (int32 i)
   {
      return (&x)[i];
   }
   
   /// Add a vector to this vector.
   void operator += (const Vec2& v)
   {
      x += v.x; y += v.y;
   }
   
   /// Subtract a vector from this vector.
   void operator -= (const Vec2& v)
   {
      x -= v.x; y -= v.y;
   }
   
   /// Multiply this vector by a scalar.
   void operator *= (float64 a)
   {
      x *= a; y *= a;
   }
   
   /// Get the length of this vector (the norm).
   float64 Length() const
   {
      return sqrt(x * x + y * y);
   }
   
   /// Get the length squared. For performance, use this instead of
   /// Vec2::Length (if possible).
   float64 LengthSquared() const
   {
      return x * x + y * y;
   }
   
   // Calculate the dot product between this vec2 and
   // another one.
   float64 Dot(const Vec2& other) const
   {
      return x*other.x + y*other.y;
   }
   
   inline static Vec2 FromPolar(float64 radius, float64 angleRads)
   {
      return Vec2(radius*cos(angleRads),radius*sin(angleRads));
   }
   
   /// Convert this vector into a unit vector. Returns the length.
   float64 Normalize()
   {
      float64 length = Length();
      if (length < DBL_MIN*2)
      {
         return 0.0;
      }
      float64 invLength = 1.0 / length;
      x *= invLength;
      y *= invLength;
      
      return length;
   }
   
   /// Get the skew vector such that dot(skew_vec, other) == cross(vec, other)
   Vec2 Skew() const
   {
      return Vec2(-y, x);
   }
   
   float64 AngleRads() const { return atan2(y,x); }
   float64 AngleDegs() const { return atan2(y,x)*180/M_PI; }
   
   inline friend std::ostream& operator<<(std::ostream& out, const Vec2& s)
   {
      out << "(" << s.x << "," << s.y << ")";
      return out;
   }
   
   std::string ToString() const
   {
      char buffer[1024];
      sprintf(buffer,"(%lf,%lf)",x,y);
      return std::string(buffer);
   }
   
   float64 x, y;
};

/// Perform the dot product on two vectors.
inline float64 b2Dot(const Vec2& a, const Vec2& b)
{
   return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline float64 b2Cross(const Vec2& a, const Vec2& b)
{
   return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a vector and a scalar. In 2D this produces
/// a vector.
inline Vec2 b2Cross(const Vec2& a, float64 s)
{
   return Vec2(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a vector. In 2D this produces
/// a vector.
inline Vec2 b2Cross(float64 s, const Vec2& a)
{
   return Vec2(-s * a.y, s * a.x);
}

/// Add two vectors component-wise.
inline Vec2 operator + (const Vec2& a, const Vec2& b)
{
   return Vec2(a.x + b.x, a.y + b.y);
}

/// Subtract two vectors component-wise.
inline Vec2 operator - (const Vec2& a, const Vec2& b)
{
   return Vec2(a.x - b.x, a.y - b.y);
}

inline Vec2 operator * (float64 s, const Vec2& a)
{
   return Vec2(s * a.x, s * a.y);
}

inline bool operator == (const Vec2& a, const Vec2& b)
{
   return a.x == b.x && a.y == b.y;
}

inline float64 b2Distance(const Vec2& a, const Vec2& b)
{
   Vec2 c = a - b;
   return c.Length();
}

inline float64 b2DistanceSquared(const Vec2& a, const Vec2& b)
{
   Vec2 c = a - b;
   return b2Dot(c, c);
}


/* Calculate the future position of a moving target so that
 * a projectile launched immediately can intercept (collide)
 * with it.
 *
 * Some situations where this might be useful for an AI to
 * make this calculation.
 *
 * 1. Shooting a projectile at a moving target.
 * 2. Launching a football or soccer ball for to a player.
 * 3. Figuring out the best position to jump towards in
 *    a platform game.
 *
 *
 * The output value, solution, is the position that the
 * intercept will occur at and the location that the
 * projectile should be launched towards.
 *
 * The function will return false if a solution cannot
 * be found.  Consider the case of a target moving away
 * from the shooter faster than the speed of the
 * projectile and you will see at least one case where
 * this can calculation may fail.
 */
bool CalculateInterceptShotPosition(const Vec2& pShooter,
                                    const Vec2& pTarget0,
                                    const Vec2& vTarget,
                                    float64 sProjectile,
                                    Vec2& solution
                                    )
{
   // This formulation uses the quadratic equation to solve
   // the intercept position.
   Vec2 R = pTarget0 - pShooter;
   float64 a = vTarget.x*vTarget.x + vTarget.y*vTarget.y - sProjectile*sProjectile;
   float64 b = 2*(R.x*vTarget.x + R.y*vTarget.y);
   float64 c = R.x*R.x + R.y*R.y;
   float64 tBullet = 0;
   
   
   // If the target and the shooter have already collided, don't bother.
   if(R.LengthSquared() < 2*DBL_MIN)
   {
      return false;
   }
   
   // If the squared velocity of the target and the bullet are the same, the equation
   // collapses to tBullet*b = -c.  If they are REALLY close to each other (float tol),
   // you could get some weirdness here.  Do some "is it close" checking?
   if(fabs(a) < 2*DBL_MIN)
   {
      // If the b value is 0, we can't get a solution.
      if(fabs(b) < 2*DBL_MIN)
      {
         return false;
      }
      tBullet = -c/b;
   }
   else
   {
      
      // Calculate the discriminant to figure out how many solutions there are.
      float64 discriminant = b*b - 4 * a * c;
      if(discriminant < 0)
      {  // All solutions are complex.
         return false;
      }
      
      if (discriminant > 0)
      {  // Two solutions.  Pick the smaller one.
         // Calculate the quadratic.
         float64 quad = sqrt(discriminant);
         float64 tBullet1 = (-b + quad)/(2*a);
         float64 tBullet2 = (-b - quad)/(2*a);
         if(tBullet1 < tBullet2 && tBullet1 >= 0)
         {
            tBullet = tBullet1;
         }
         else
         {
            tBullet = tBullet2;
         }
      }
      else
      {
         tBullet = -b / (2*a);
      }
   }
   // If the time is negative, we can't get there from here.
   if(tBullet < 0)
   {
      return false;
   }
   // Calculate the intercept position.
   solution = pTarget0 + tBullet*vTarget;
   
   return true;
}

int main(int argc, const char * argv[])
{
   // Initial Conditions
   Vec2 pShooter(0.0,0.0);       // m
   Vec2 pTarget0(2.0,2.0);       // m
   Vec2 vTarget(0.0,-2.0);        // m/s
   float64 sProjectile = 2.0;    // m/s
   float64 tMin = 0.0;           // Min Firing Time
   float64 tMax = 10.0;          // Max Firing Time
   
   const int32 STEPS_PER_SEC = 100;
   Vec2 solution;
   
   bool foundSolution = CalculateInterceptShotPosition(pShooter,
                                                       pTarget0,
                                                       vTarget,
                                                       sProjectile,
                                                       solution);
   
   if(!foundSolution)
   {
      cout << "----------------------------------------------" << endl;
      cout << " NO SOLUTION FOUND!!!" << endl;
      cout << "----------------------------------------------" << endl;
   }
   else
   {
      // Shooting variables.
      Vec2 vProjectile = solution - pShooter;
      vProjectile.Normalize();
      vProjectile *= sProjectile;
      
      const float64 COLL_DIST = 0.05;
      const float64 COLL_DIST_SQ = COLL_DIST*COLL_DIST;
      
      cout << "----------------------------------------------" << endl;
      cout << " Starting Simulation [" << tMin << "," << tMax << "] seconds" << endl;
      cout << " Intercept Position: " << solution.ToString() << endl;
      cout << " Projectile Velocity: " << vProjectile.ToString() << endl;
      cout << "----------------------------------------------" << endl;
      
      float64 fireTime = 0.0;
      bool hitTarget = false;
      
      float64 now = 0;
      
      cout << endl;
      cout << "----------------------------------------------" << endl;
      cout << " FIRING PROJETILE " << endl;
      cout << "----------------------------------------------" << endl;
      fireTime = now;
      
      for(int idx = 0; idx <= tMax*STEPS_PER_SEC && !hitTarget; ++idx)
      {
         float64 now = idx * 1.0/STEPS_PER_SEC;
         // Update the position of the target.
         Vec2 pTarget = pTarget0 + now*vTarget;
         Vec2 pProjectile = pShooter + (now - fireTime)*vProjectile;
         printf("Step[%05d] [%.2f sec] ",idx,now);
         // Update the position of the rotation.
         cout << "pProjectile: " << pProjectile.ToString() << " ";
         cout << "pTarget: " << pTarget.ToString() << " ";
         cout << "DIST: " << (pTarget-pProjectile).Length() << " meters";
         cout << endl;
         if(!hitTarget && (pProjectile - pTarget).LengthSquared() < COLL_DIST_SQ)
         {
            hitTarget = true;
            cout << "----------------------------------------------" << endl;
            cout << " HIT TARGET!!!!! " << endl;
            cout << "----------------------------------------------" << endl;
         }
      }
   }
}

