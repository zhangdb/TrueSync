/* Copyright (C) <2009-2011> <Thorben Linneweber, Jitter Physics>
* 
*  This software is provided 'as-is', without any express or implied
*  warranty.  In no event will the authors be held liable for any damages
*  arising from the use of this software.
*
*  Permission is granted to anyone to use this software for any purpose,
*  including commercial applications, and to alter it and redistribute it
*  freely, subject to the following restrictions:
*
*  1. The origin of this software must not be misrepresented; you must not
*      claim that you wrote the original software. If you use this software
*      in a product, an acknowledgment in the product documentation would be
*      appreciated but is not required.
*  2. Altered source versions must be plainly marked as such, and must not be
*      misrepresented as being the original software.
*  3. This notice may not be removed or altered from any source distribution. 
*/

namespace TrueSync.Physics3D {

    #region public class ContactSettings
    public class ContactSettings
    {
        public enum MaterialCoefficientMixingType { TakeMaximum, TakeMinimum, UseAverage }

        internal FP maximumBias = 10;
        internal FP bias = 25 * FP.EN2;
        internal FP minVelocity = FP.EN3;
        internal FP allowedPenetration = FP.EN3;
        internal FP breakThreshold = FP.EN3;

        internal MaterialCoefficientMixingType materialMode = MaterialCoefficientMixingType.UseAverage;

        public FP MaximumBias { get { return maximumBias; } set { maximumBias = value; } }

        public FP BiasFactor { get { return bias; } set { bias = value; } }

        public FP MinimumVelocity { get { return minVelocity; } set { minVelocity = value; } }

        public FP AllowedPenetration { get { return allowedPenetration; } set { allowedPenetration = value; } }

        public FP BreakThreshold { get { return breakThreshold; } set { breakThreshold = value; } }

        public MaterialCoefficientMixingType MaterialCoefficientMixing { get { return materialMode; } set { materialMode = value; } }
    }
    #endregion


    /// <summary>
    /// </summary>
    public class Contact : IConstraint
    {
        public ContactSettings settings;

		public RigidBody body1, body2;

		public TSVector normal, tangent;

		public TSVector realRelPos1, realRelPos2;
		public TSVector relativePos1, relativePos2;
		public TSVector p1, p2;

		public FP accumulatedNormalImpulse = FP.Zero;
		public FP accumulatedTangentImpulse = FP.Zero;

		public FP penetration = FP.Zero;
		public FP initialPen = FP.Zero;

		public FP staticFriction, dynamicFriction, restitution;
		public FP friction = FP.Zero;

		public FP massNormal = FP.Zero, massTangent = FP.Zero;
		public FP restitutionBias = FP.Zero;

		public bool newContact = false;

		public bool body1IsMassPoint;
		public bool body2IsMassPoint;

		public FP lostSpeculativeBounce = FP.Zero;
		public FP speculativeVelocity = FP.Zero;

        /// <summary>
        /// A contact resource pool.
        /// </summary>
        public static readonly ResourcePool<Contact> Pool =
            new ResourcePool<Contact>();

		public FP lastTimeStep = FP.PositiveInfinity;

        #region Properties
        public FP Restitution
        {
            get { return restitution; }
            set { restitution = value; }
        }

        public FP StaticFriction
        {
            get { return staticFriction; }
            set { staticFriction = value; }
        }

        public FP DynamicFriction
        {
            get { return dynamicFriction; }
            set { dynamicFriction = value; }
        }

        /// <summary>
        /// The first body involved in the contact.
        /// </summary>
        public RigidBody Body1 { get { return body1; } }

        /// <summary>
        /// The second body involved in the contact.
        /// </summary>
        public RigidBody Body2 { get { return body2; } }

        /// <summary>
        /// The penetration of the contact.
        /// </summary>
        public FP Penetration { get { return penetration; } }

        /// <summary>
        /// The collision position in world space of body1.
        /// </summary>
        public TSVector Position1 { get { return p1; } }

        /// <summary>
        /// The collision position in world space of body2.
        /// </summary>
        public TSVector Position2 { get { return p2; } }

        /// <summary>
        /// The contact tangent.
        /// </summary>
        public TSVector Tangent { get { return tangent; } }

        /// <summary>
        /// The contact normal.
        /// </summary>
        public TSVector Normal { get { return normal; } }
        #endregion

        /// <summary>
        /// Calculates relative velocity of body contact points on the bodies.
        /// </summary>
        public TSVector CalculateRelativeVelocity()
        {
            var v2 = TSVector.Cross(body2.angularVelocity, relativePos2) + body2.linearVelocity;
            var v1 = TSVector.Cross(body1.angularVelocity, relativePos1) + body1.linearVelocity;
            return v2 - v1;
        }

        /// <summary>
        /// Solves the contact iteratively.
        /// </summary>
        public void Iterate()
        {
            //body1.linearVelocity = JVector.Zero;
            //body2.linearVelocity = JVector.Zero;
            //return;
          
            if (body1.isStatic && body2.isStatic) return;

            var dv = body2.linearVelocity - body1.linearVelocity;

            if (!body1IsMassPoint)
            {
                dv -= TSVector.Cross(body1.angularVelocity, relativePos1);
            }

            if (!body2IsMassPoint)
            {
                dv += TSVector.Cross(body2.angularVelocity, relativePos2);
            }

            // this gets us some performance
            if (dv.sqrMagnitude < settings.minVelocity * settings.minVelocity) return;

          
            FP vn = TSVector.Dot(normal, dv);
            FP normalImpulse = massNormal * (-vn + restitutionBias + speculativeVelocity);

            FP oldNormalImpulse = accumulatedNormalImpulse;
            accumulatedNormalImpulse = oldNormalImpulse + normalImpulse;
            if (accumulatedNormalImpulse < FP.Zero) accumulatedNormalImpulse = FP.Zero;
            normalImpulse = accumulatedNormalImpulse - oldNormalImpulse;

            FP vt = TSVector.Dot(dv, tangent);
            FP maxTangentImpulse = friction * accumulatedNormalImpulse;
            FP tangentImpulse = massTangent * -vt;

            FP oldTangentImpulse = accumulatedTangentImpulse;
            accumulatedTangentImpulse = oldTangentImpulse + tangentImpulse;
            if (accumulatedTangentImpulse < -maxTangentImpulse) accumulatedTangentImpulse = -maxTangentImpulse;
            else if (accumulatedTangentImpulse > maxTangentImpulse) accumulatedTangentImpulse = maxTangentImpulse;

            tangentImpulse = accumulatedTangentImpulse - oldTangentImpulse;

            
            // Apply contact impulse
            TSVector impulse = normal * normalImpulse + tangent * tangentImpulse;
            ApplyImpulse(impulse);

        }

        public FP AppliedNormalImpulse { get { return accumulatedNormalImpulse; } }
        public FP AppliedTangentImpulse { get { return accumulatedTangentImpulse; } }

        /// <summary>
        /// The points in wolrd space gets recalculated by transforming the
        /// local coordinates. Also new penetration depth is estimated.
        /// </summary>
        public void UpdatePosition()
        {
            if (body1IsMassPoint)
            {
                p1 = realRelPos1 + body1.position;
            }
            else
            {
                p1 = body1.position + body1.Orientation.Multiply(realRelPos1);
            }

            if (body2IsMassPoint)
            {
                p2 = realRelPos2 + body2.position;
            }
            else
            {
                p2 = body2.position + body2.Orientation.Multiply(realRelPos2);
            }
            penetration = TSVector.Dot(p1 - p2, normal);
        }

        /// <summary>
        /// An impulse is applied an both contact points.
        /// </summary>
        /// <param name="impulse">The impulse to apply.</param>
        public void ApplyImpulse(TSVector impulse)
        {
            #region INLINE - HighFrequency
            //JVector temp;

            if (!body1.isStatic)
            {
                body1.ApplyImpulse(-impulse, relativePos1);
            }

            if (!body2.isStatic)
            {

                body2.ApplyImpulse(impulse, relativePos2);

            }
            #endregion
        }

        public static FP CalcK(RigidBody obj, TSVector relativePos, TSVector _normal)
        {
            if (!obj.IsStatic)
            {
                var v = TSVector.Cross(relativePos, _normal);
                v = obj.invInertiaWorld.TransposedMultiply(v);
                v = TSVector.Cross(v, relativePos);
                return obj.inverseMass + TSVector.Dot(v, _normal);
            }
            return 0;
        }

        /// <summary>
        /// PrepareForIteration has to be called before <see cref="Iterate"/>.
        /// </summary>
        /// <param name="timestep">The timestep of the simulation.</param>
        public void PrepareForIteration(FP timestep)
        {
            var dv = CalculateRelativeVelocity();

            var kn = CalcK(body1, relativePos1, normal) +
                     CalcK(body2, relativePos2, normal);

            massNormal = FP.One / kn;

            tangent = dv - TSVector.Dot(dv, normal) * normal;
            tangent.Normalize();


            var kt = CalcK(body1, relativePos1, tangent) +
                     CalcK(body2, relativePos2, tangent);

            massTangent = FP.One / kt;
            
            restitutionBias = lostSpeculativeBounce;

            speculativeVelocity = FP.Zero;

            FP relNormalVel = TSVector.Dot(normal, dv);

            if (Penetration > settings.allowedPenetration)
            {
                restitutionBias = settings.bias * (FP.One / timestep) * TSMath.Max(FP.Zero, Penetration - settings.allowedPenetration);
                restitutionBias = TSMath.Clamp(restitutionBias, FP.Zero, settings.maximumBias);
              //  body1IsMassPoint = body2IsMassPoint = false;
            }
      

            FP timeStepRatio = timestep / lastTimeStep;
            accumulatedNormalImpulse *= timeStepRatio;
            accumulatedTangentImpulse *= timeStepRatio;

            {
                // Static/Dynamic friction
                FP relTangentVel = -TSVector.Dot(tangent, dv);
                FP tangentImpulse = massTangent * relTangentVel;
                FP maxTangentImpulse = -staticFriction * accumulatedNormalImpulse;

                if (tangentImpulse < maxTangentImpulse)
                    friction = dynamicFriction;
                else
                    friction = staticFriction;
            }

            // Simultaneos solving and restitution is simply not possible
            // so fake it a bit by just applying restitution impulse when there
            // is a new contact.
            if (relNormalVel < -FP.One && newContact)
            {
                restitutionBias = TSMath.Max(-restitution * relNormalVel, restitutionBias);
            }

            // Speculative Contacts!
            // if the penetration is negative (which means the bodies are not already in contact, but they will
            // be in the future) we store the current bounce bias in the variable 'lostSpeculativeBounce'
            // and apply it the next frame, when the speculative contact was already solved.
            if (penetration < -settings.allowedPenetration)
            {
                speculativeVelocity = penetration / timestep;

                lostSpeculativeBounce = restitutionBias;
                restitutionBias = FP.Zero;
            }
            else
            {
                lostSpeculativeBounce = FP.Zero;
            }

            var impulse = normal * accumulatedNormalImpulse + tangent * accumulatedTangentImpulse;
            ApplyImpulse(impulse);

            lastTimeStep = timestep;

            newContact = false;
        }

       


        /// <summary>
        /// Initializes a contact.
        /// </summary>
        /// <param name="body1">The first body.</param>
        /// <param name="body2">The second body.</param>
        /// <param name="point1">The collision point in worldspace</param>
        /// <param name="point2">The collision point in worldspace</param>
        /// <param name="n">The normal pointing to body2.</param>
        /// <param name="penetration">The estimated penetration depth.</param>
        public void Initialize(RigidBody body1, RigidBody body2, ref TSVector point1, ref TSVector point2, ref TSVector n,
            FP penetration, bool newContact, ContactSettings settings)
        {
           
            this.body1 = body1;  this.body2 = body2;
            this.normal = n; normal.Normalize();
            this.p1 = point1; this.p2 = point2;

            this.newContact = newContact;

            TSVector.Subtract(p1, body1.position, out relativePos1);
            TSVector.Subtract(p2, body2.position, out relativePos2);

            //TSVector.Transform(relativePos1, body1.invOrientation, out realRelPos1);
            //TSVector.Transform(relativePos2, body2.invOrientation, out realRelPos2);
            realRelPos1 = body1.Orientation.Multiply(relativePos1);
            realRelPos2 = body2.Orientation.Multiply(relativePos2);

            this.initialPen = penetration;
            this.penetration = penetration;

            body1IsMassPoint = body1.isParticle;
            body2IsMassPoint = body2.isParticle;

            // Material Properties
            if (newContact)
            {

                accumulatedNormalImpulse = FP.Zero;
                accumulatedTangentImpulse = FP.Zero;

                lostSpeculativeBounce = FP.Zero;

                switch (settings.MaterialCoefficientMixing)
                {
                    case ContactSettings.MaterialCoefficientMixingType.TakeMaximum:
                        staticFriction = TSMath.Max(body1.staticFriction, body2.staticFriction);
                        dynamicFriction = TSMath.Max(body1.staticFriction, body2.staticFriction);
                        restitution = TSMath.Max(body1.restitution, body2.restitution);
                        break;
                    case ContactSettings.MaterialCoefficientMixingType.TakeMinimum:
                        staticFriction = TSMath.Min(body1.staticFriction, body2.staticFriction);
                        dynamicFriction = TSMath.Min(body1.staticFriction, body2.staticFriction);
                        restitution = TSMath.Min(body1.restitution, body2.restitution);
                        break;
                    case ContactSettings.MaterialCoefficientMixingType.UseAverage:
                        staticFriction = (body1.staticFriction + body2.staticFriction) * FP.Half;
                        dynamicFriction = (body1.staticFriction + body2.staticFriction) * FP.Half;
                        restitution = (body1.restitution + body2.restitution) * FP.Half;
                        break;
                }

            }

            this.settings = settings;
        }
    }
}