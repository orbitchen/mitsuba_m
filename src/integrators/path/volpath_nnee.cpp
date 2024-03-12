/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

const bool enableRIS = true;

static float mediumPhaseG = 0.0f;
static Point lightPosition = Point(-0.4, -1.6, 0.24);
static Spectrum lightIntensity = Spectrum(40.0);

void __debug_print_point (Point a, const char* name)
{
    printf("%s: %f, %f, %f\n",name,a.x,a.y,a.z);
};

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

/*!\plugin{volpath}{Extended volumetric path tracer}
 * \order{4}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *        which the implementation will start to use the ``russian roulette''
 *        path termination criterion. \default{\code{5}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See
 *        page~\pageref{sec:strictnormals} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 *
 * This plugin provides a volumetric path tracer that can be used to
 * compute approximate solutions of the radiative transfer equation.
 * Its implementation makes use of multiple importance sampling to
 * combine BSDF and phase function sampling with direct illumination
 * sampling strategies. On surfaces, it behaves exactly
 * like the standard path tracer.
 *
 * This integrator has special support for \emph{index-matched} transmission
 * events (i.e. surface scattering events that do not change the direction
 * of light). As a consequence, participating media enclosed by a stencil shape (see
 * \secref{shapes} for details) are rendered considerably more efficiently when this
 * shape has \emph{no}\footnote{this is what signals to Mitsuba that the boundary is
 * index-matched and does not interact with light in any way. Alternatively,
 * the \pluginref{mask} and \pluginref{thindielectric} BSDF can be used to specify
 * index-matched boundaries that involve some amount of interaction.} BSDF assigned
 * to it (as compared to, say, a \pluginref{dielectric} or \pluginref{roughdielectric} BSDF).
 *
 * \remarks{
 *    \item This integrator will generally perform poorly when rendering
 *      participating media that have a different index of refraction compared
 *      to the surrounding medium.
 *    \item This integrator has poor convergence properties when rendering
 *      caustics and similar effects. In this case, \pluginref{bdpt} or
 *      one of the photon mappers may be preferable.
 * }
 */

class VolumetricPathTracerNNEE : public MonteCarloIntegrator {
public:
    VolumetricPathTracerNNEE(const Properties &props) : MonteCarloIntegrator(props) { }

    /// Unserialize from a binary data stream
    VolumetricPathTracerNNEE(Stream *stream, InstanceManager *manager)
     : MonteCarloIntegrator(stream, manager) { }

    /// @brief sample on an arbitrary (positional) distribution. 
    /// The distribution don't have to be inside medium.
    /// @param mRec medium sample record. records t, sigmaA, sigmaS, time, medium, p and transmittance.
    /// @param scene scene representation.
    /// @return positional pdf.
    Float sampleSpaceScatteringPoint(MediumSamplingRecord &mRec, const Scene* scene) const
    {
        // TODO
    }

    Float sampleSurfacePdf(MediumSamplingRecord &mRec, RadianceQueryRecord &rRec, const Scene *scene) const
    {
        // return mRec.pdfFailure;
        return mRec.pdfFailure;
    }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        MediumSamplingRecord mRec;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        Float eta = 1.0f;

        {
            auto tempEmitter = scene->getEmitters()[0].get();
            const Transform &trafo = tempEmitter->getWorldTransform()->eval(0.0f);
            lightPosition = trafo(Point(0,0,0));
            lightIntensity = tempEmitter->getIntensity();
        }

        // printf("begin li!\n");

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);

        // ASSUME that the camera is not inside medium. nothing to do here

        Spectrum throughput(1.0f);
        bool scattered = false;

        // Is current scattering point the first scattering point inside medium?
        bool enterMedium = false;

        // record infomation about previous NNEE scattering point
        struct scatteringPointNNEE
        {
            Point sampleCentre;
            Vector PhaseAxis; // REVERSED direction. i.e., outwards sampleCentre
            Float sampleG;

            Point scatteringPoint;
            Float scatteringPdf; // measured by d(omega) * d(distance)

            Ray ray; // original sample ray. could be used by equiangular

            Point surfacePoint;
            Spectrum di; // direct illumination. i.e., Tr * Le * G for point light source
            Spectrum LeDotG; // L_e * G term

            bool selected = false;
        };

        std::vector<scatteringPointNNEE> previousScatter;
        std::vector<scatteringPointNNEE> previousCasualScatter;
        const Float scatterTrMin = 1e-6;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {

            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
            
            // if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) 
            if (rRec.medium)
            {
                // OUT_DEBUG("medium intersection");
                /* Sample the integral
                \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
                */

                mediumPhaseG = rRec.medium->getPhaseFunction()->getMeanCosine();
                
                if (rRec.depth >= m_maxDepth && m_maxDepth != -1) // No more scattering events allowed
                    break;

                // throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                const PhaseFunction *phase = mRec.getPhaseFunction();

                /* ==================================================================== */
                /*                          Luminaire sampling                          */
                /* ==================================================================== */

                /* Estimate the single scattering component if this is requested */
                DirectSamplingRecord dRec(mRec.p, mRec.time);
                if (enterMedium)
                {
                    enterMedium = false;

                    if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
                        int interactions = m_maxDepth - rRec.depth - 1;

                        Spectrum value = scene->sampleAttenuatedEmitterDirect(
                                dRec, rRec.medium, interactions,
                                rRec.nextSample2D(), rRec.sampler);

                        if (!value.isZero()) {
                            const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                            /* Evaluate the phase function */
                            PhaseFunctionSamplingRecord pRec(mRec, -ray.d, dRec.d);
                            Float phaseVal = phase->eval(pRec);

                            if (phaseVal != 0) {
                                /* Calculate prob. of having sampled that direction using
                                phase function sampling */
                                Float phasePdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                        ? phase->pdf(pRec) : (Float) 0.0f;

                                /* Weight using the power heuristic */
                                const Float weight = miWeight(dRec.pdf, phasePdf);
                                Li += throughput * value * phaseVal * weight;
                            }
                        }
                    }
                }

                // Next-Next Event Estimation (NNEE) & MIS
                /* ========================================================================== */
                /*                         Phase function sampling (NNEE)                     */
                /* ========================================================================== */

                // printf("new loop begin scattering point: %f\n",mRec.p.x);

                // const g
                const Float newPhaseFunctionG = 0.8;

                /* note:
                    mediumSampleRecord -> distance sampling, records scattering info of current scattering point
                    phaseSampleRecord -> direction sampling, records scattering direction & pdf, all direction outwards
                */

                // return phase value & phase pdf (same value)
                auto newPhaseSample = [&newPhaseFunctionG] (Vector dir, PhaseFunctionSamplingRecord& pRec, Sampler* sampler) -> Float
                {
                    Float m_g = newPhaseFunctionG;

                    Point2 sample(sampler->next2D());

                    Float cosTheta;
                    if (std::abs(m_g) < Epsilon) {
                        cosTheta = 1 - 2*sample.x;
                    } else {
                        Float sqrTerm = (1 - m_g * m_g) / (1 - m_g + 2 * m_g * sample.x);
                        cosTheta = (1 + m_g * m_g - sqrTerm * sqrTerm) / (2 * m_g);
                    }

                    Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
                        sinPhi, cosPhi;

                    math::sincos(2*M_PI*sample.y, &sinPhi, &cosPhi);

                    pRec.wo = Frame(-dir).toWorld(Vector(
                        sinTheta * cosPhi,
                        sinTheta * sinPhi,
                        cosTheta
                    ));

                    Float temp = 1.0f + m_g*m_g + 2.0f * m_g * dot(dir, pRec.wo);
                    return INV_FOURPI * (1 - m_g*m_g) / (temp * std::sqrt(temp));
                }; 

                auto newPhasePdfBiVec = [] (Vector dir, Vector wo, Float g) -> Float
                {
                    Float m_g= g;
                    Float temp = 1.0f + m_g*m_g + 2.0f * m_g * dot(dir, wo);
                    return INV_FOURPI * (1 - m_g*m_g) / (temp * std::sqrt(temp));
                };

                auto newPhasePdf = [newPhaseFunctionG] (Vector dir, PhaseFunctionSamplingRecord& pRec) -> Float
                {
                    Float m_g= newPhaseFunctionG;
                    Float temp = 1.0f + m_g*m_g + 2.0f * m_g * dot(dir, pRec.wo);
                    return INV_FOURPI * (1 - m_g*m_g) / (temp * std::sqrt(temp));
                };

                std::vector<PhaseFunctionSamplingRecord> pRecordArray; // direction sample record
                std::vector<Ray> rayRecordArray; // point & direction from origin point -> sample point
                std::vector<Intersection> mItsArray; // intersection record
                std::vector<Float> mSolidAnglePdfArray; // direction pdf record
                std::vector<MediumSamplingRecord> mRecordArray; // point sample record
                std::vector<Intersection> diItsArray; // direct illumination intersection (at the surface of medium) record
                std::vector<Spectrum> diArray; // direct illumination record (L_e * G * Tr)
                std::vector<Spectrum> LeDotGArray; // Tr for di

                // printf("NNEE begin\n");

                // sample original phase function
                {
                    Float phasePdf;
                    PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                    Float phaseVal = phase->sample(pRec, phasePdf, rRec.sampler);
                    if (phaseVal == 0)
                        break;
                    // phaseVal = phase / p
                    // throughput *= phaseVal;

                    /* Trace a ray in this direction */
                    Ray initialRay = ray;
                    ray = Ray(mRec.p, pRec.wo, ray.time);
                    ray.mint = Epsilon;
                    rayRecordArray.push_back(ray);
                    ray.mint = 0;

                    Spectrum value(0.0f);
                    Intersection firstIts;
                    Intersection itsNew;
                    // rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
                    //     m_maxDepth - rRec.depth - 1, ray, itsNew, dRec, value, &firstIts);
                    scene->rayIntersect(ray,itsNew);

                    /* If a luminaire was hit, estimate the local illumination and
                    weight using the power heuristic */
                    // NO contributions here -> Assume point light source
                    // if (!value.isZero() && (rRec.type & RadianceQueryRecord::EDirectMediumRadiance)) {
                    //     const Float emitterPdf = scene->pdfEmitterDirect(dRec);
                    //     Li += throughput * phaseVal * value * miWeight(phasePdf, emitterPdf);
                    // }

                    ray = initialRay;

                    pRecordArray.push_back(pRec);
                    mItsArray.push_back(itsNew);
                    mSolidAnglePdfArray.push_back(phasePdf);
                }

                // printf("sample origin phase function done\n");

                Vector towardsLight = normalize(lightPosition - mRec.p);

                // sample heuristic phase function
                {
                    PhaseFunctionSamplingRecord pRecNew(mRec, -ray.d);
                    Float phaseVal = newPhaseSample(-towardsLight,pRecNew, rRec.sampler);

                    Intersection itsNew;

                    Ray initialRay = ray;
                    ray = Ray(mRec.p, pRecNew.wo, ray.time);
                    ray.mint = Epsilon;
                    rayRecordArray.push_back(ray);
                    ray.mint = 0;

                    Spectrum value(0.0f);
                    Intersection firstIts;
                    // rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
                    //     m_maxDepth - rRec.depth - 1, ray, itsNew, dRec, value, &firstIts);
                    scene->rayIntersect(ray,itsNew);

                    ray = initialRay;

                    pRecordArray.push_back(pRecNew);
                    mItsArray.push_back(itsNew);
                    mSolidAnglePdfArray.push_back(phaseVal);
                }

                // printf("sample heuristic phase function done\n");

                // ==============================
                // estimate NNEE contributions
                // ==============================

                auto medium = mRec.medium;

                auto multipleMisWeight = [] (std::vector<Float> pdfs, Float currPdf) -> Float
                {
                    Float sum = 0.0f;
                    for (auto& item: pdfs) sum += item * item;
                    return (currPdf * currPdf) / sum;
                };

                auto calculateDi = [&diItsArray, &scene, &rRec, this] (const MediumSamplingRecord& mRec, const Medium* m, Sampler* sampler, Spectrum& LeDotG) -> Spectrum
                {
                    // di includes transmittance
                    int inter = m_maxDepth;
                    diItsArray.emplace_back(Intersection());
                    Spectrum tr = scene->evalTransmittance(mRec.p, false, lightPosition, false, 0.0f, m, inter, sampler, &diItsArray[diItsArray.size()-1]);

                    LeDotG = lightIntensity * (1.0f / distanceSquared(mRec.p, lightPosition));
                    return tr * lightIntensity * (1.0f / distanceSquared(mRec.p, lightPosition));
                };

                auto scatteringJacobi = [] (const Point& originSamplePoint, const Point& scatteringPoint, const Point& currentSamplePoint) -> Float
                {
                    // return distanceSquared(originSamplePoint, scatteringPoint) / distanceSquared(scatteringPoint, currentSamplePoint);
                    return distanceSquared(scatteringPoint, currentSamplePoint) / distanceSquared(originSamplePoint, scatteringPoint);
                };
                
                auto reuseMisPdf = [newPhaseFunctionG, newPhasePdfBiVec, medium, scatteringJacobi] 
                (const scatteringPointNNEE& sample, const Point& scatteringPoint, const Point& currentSamplePoint, const Ray& ray) -> Float
                {
                    // Vector sampleDir = normalize(scatteringPoint - sample.sampleCentre);
                    Float sampleDistance = distance(scatteringPoint, sample.sampleCentre);

                    if (sampleDistance == 0.0f) return 1.0f;
                    Float phasePdf = newPhasePdfBiVec(sample.PhaseAxis, normalize(scatteringPoint - sample.sampleCentre), sample.sampleG);

                    MediumSamplingRecord temp;
                    temp.t = sampleDistance;
                    Float distancePdf = medium->pdfDistanceMultipleScattering(temp);
                    Float J = scatteringJacobi(sample.sampleCentre, scatteringPoint, currentSamplePoint);

                    return phasePdf * distancePdf * J;
                };
                int NNEESampleCount = 2;

                // select previous NNEE samples
                std::vector<scatteringPointNNEE*> selectedScatteringPoint;
                for(auto& item: previousScatter){
                    if (item.scatteringPoint != mRec.p && medium->evalTransmittance(Ray(ray, 0.0f, distance(item.scatteringPoint, mRec.p))).max() > scatterTrMin) {
                        selectedScatteringPoint.push_back(&item);
                    }

                    if(selectedScatteringPoint.size() >= 6) break; // pick at most 6 samples
                }

                Spectrum nneeContribution = Spectrum(0.0f);

                // fulfill info & calculate current (2 scattering point) contributions
                for (int i=0;i< NNEESampleCount; i++){

                    // sample multiple distances
                    MediumSamplingRecord mRecNext;
                    // rRec.medium->sampleDistanceMultipleScattering(
                    //     Ray(Ray(mRec.p, pRecordArray[i].wo, ray.time), 0, mItsArray[i].t), 
                    //     mRecNext, 
                    //     rRec.sampler);
                    rRec.medium->sampleDistanceMultipleScattering(
                        Ray(Ray(mRec.p, pRecordArray[i].wo, ray.time), 0, mItsArray[i].t), 
                        mRecNext, 
                        rRec.sampler);

                    mRecordArray.push_back(mRecNext);

                    // if(i==0)
                    //     printf("phase/distance scattering point: %f\n", mRecNext.p.x);

                    // if(i==1)
                    //     printf("proxy scattering point: %f\n", mRecNext.p.x);

                    // do intersection & calculate L_e
                    Spectrum LeDotG;
                    Spectrum di = calculateDi(mRecNext, mRec.medium, rRec.sampler,LeDotG);
                    diArray.push_back(di);
                    LeDotGArray.push_back(LeDotG);

                    // mix samples using MIS
                    std::vector<Float> pdfs;

                    pdfs.push_back(mSolidAnglePdfArray[0]);
                    pdfs.push_back(mSolidAnglePdfArray[1]);

                    pdfs[0]*=mRec.medium->pdfDistanceMultipleScattering(mRecNext);
                    pdfs[1]*=mRec.medium->pdfDistanceMultipleScattering(mRecNext);

                    for(auto& item: selectedScatteringPoint){
                        pdfs.push_back(reuseMisPdf(*item, mRecNext.p, mRec.p, rayRecordArray[i]));
                    }

                    Float misWeight = multipleMisWeight(pdfs, pdfs[i]);

                    Float phase2 = mRec.getPhaseFunction()->eval(-pRecordArray[i].wo, normalize(lightPosition - mRecNext.p));
                    // thp * phase1 * tr1 * phase2 * tr2 * Le * mis / (phase 1 pdf * distance 1 pdf)
                    // phase1 / phase 1 pdf == 1.0f
                    // di == tr2 * Le
                    // Spectrum contribution = throughput * 1.0f * mRec.sigmaS * mRecNext.transmittance * phase2 * di * misWeight / (mRecNext.pdfSuccess);
                    Spectrum contribution = throughput * 1.0f * mRec.sigmaS * mRecNext.transmittance * phase2 * di * misWeight / (mRecNext.pdfSuccess);
                    nneeContribution += contribution;
                }

                // calculate reused (indefinite scattering point) contributions
                int selectedScatteringPointCount = selectedScatteringPoint.size();
                for(int si = 0; si < selectedScatteringPointCount; si++){

                    auto& item = *selectedScatteringPoint[si];
                    // if ((item.scatteringPoint - mRec.p).length()==0.0f) printf("v.length == 0\n");
                    Vector scatterDir = normalize(item.scatteringPoint - mRec.p);
                    Float scatterDistance = distance(item.scatteringPoint, mRec.p);
                    std::vector<Float> pdfs;

                    Float phaseValue = mRec.getPhaseFunction()->eval(-ray.d, scatterDir);
                    Spectrum tr = medium->evalTransmittance(Ray(ray, 0.0f, scatterDistance));

                    pdfs.push_back(phaseValue);
                    pdfs.push_back(newPhasePdfBiVec(-towardsLight, scatterDir, item.sampleG));

                    MediumSamplingRecord temp;
                    temp.t = scatterDistance;
                    pdfs[0]*=mRec.medium->pdfDistanceMultipleScattering(temp);
                    pdfs[1]*=mRec.medium->pdfDistanceMultipleScattering(temp);
                    

                    for(auto& ss: selectedScatteringPoint){
                        pdfs.push_back(reuseMisPdf(*ss, item.scatteringPoint, mRec.p, item.ray));
                    }

                    Float misWeight = multipleMisWeight(pdfs, pdfs[2+si]);

                    Float phase2 = mRec.getPhaseFunction()->eval(-scatterDir, normalize(lightPosition - item.scatteringPoint));
                    
                    nneeContribution += throughput * phaseValue * mRec.sigmaS * tr * phase2 * item.di * misWeight / (item.scatteringPdf * scatteringJacobi(item.sampleCentre, item.scatteringPoint, mRec.p));                
                }

                // add current NNEE sample to list
                if (mRecordArray[1].p != mRec.p)
                {
                    scatteringPointNNEE newScattering;
                    newScattering.sampleCentre = mRec.p;
                    newScattering.PhaseAxis = -towardsLight;
                    newScattering.sampleG = newPhaseFunctionG;
                    newScattering.scatteringPdf = mSolidAnglePdfArray[1] * mRecordArray[1].pdfSuccess;
                    newScattering.scatteringPoint = mRecordArray[1].p;
                    newScattering.di = diArray[1];
                    newScattering.ray = rayRecordArray[1];
                    newScattering.LeDotG = LeDotGArray[1];
                    newScattering.surfacePoint = diItsArray[1].p;
                    previousScatter.push_back(newScattering);

                    if (enableRIS)
                    {
                        // previousScatter.push_back(newScattering);
                        selectedScatteringPoint.push_back(&previousScatter[previousScatter.size()-1]);
                    }

                    // printf("add new NNEE point: %f with scattering centre: %f\n",newScattering.scatteringPoint.x, newScattering.sampleCentre.x);
                }

                if (enableRIS && mRecordArray[0].p != mRec.p)
                {
                    scatteringPointNNEE newScattering;
                    newScattering.sampleCentre = mRec.p;
                    newScattering.PhaseAxis = -ray.d;
                    newScattering.sampleG = mediumPhaseG;
                    newScattering.scatteringPdf = mSolidAnglePdfArray[0] * mRecordArray[0].pdfSuccess;
                    newScattering.scatteringPoint = mRecordArray[0].p;
                    newScattering.di = diArray[0];
                    newScattering.ray = rayRecordArray[0];
                    newScattering.LeDotG = LeDotGArray[0];
                    newScattering.surfacePoint = diItsArray[0].p;
                    previousCasualScatter.push_back(newScattering);

                    // printf("add casual NNEE point: %f with scattering centre: %f\n",newScattering.scatteringPoint.x, newScattering.sampleCentre.x);
                }

                // printf("for end\n");
                // reset intersection & record
                // TODO: resample

                auto bssrdf = [] (const Float albedo, const Float sigma_t,  const Float d) -> Float
                {
                    Float s = 1.85 - albedo + 7.0f * (albedo - 0.8) * (albedo - 0.8) * abs(albedo - 0.8);
                    Float w = albedo * s;
                    Float l = 1.0f / sigma_t;
                    // printf("s: %f, w: %f, l: %f, albedo: %f, sigma_t: %f, d: %f\n", s,w,l,albedo,sigma_t,d);
                    return w * (exp(- s * d / l) + exp(- s * d / (3.0f * l))) / (8.0f * M_PI * l * d);

                    // printf("bssrdf: %f\n", w * (exp(- s * d / l) + exp(- s * d / (3.0f * l))) / (8.0f * M_PI * l * d));
                    // return 1.0f;
                };

                auto pHat = [bssrdf, medium] (const scatteringPointNNEE& sample, const Point& curScatteringP, const Vector& rayDir) ->Float
                {
                    Spectrum albedo = medium->getSigmaS() / (medium->getSigmaT());
                    Spectrum sigma_t = medium->getSigmaT();
                    Float d = distance(sample.surfacePoint, sample.scatteringPoint);
                    // printf("surface point: %f, scattering point: %f\n",sample.surfacePoint.x, sample.scatteringPoint.x);

                    Float p = medium->getPhaseFunction()->eval(-rayDir, normalize(sample.scatteringPoint - curScatteringP));
                    Ray fakeRay; fakeRay.mint = 0.0f; fakeRay.maxt = distance(sample.scatteringPoint, curScatteringP);
                    Spectrum tr = medium->evalTransmittance(fakeRay);
                    Spectrum retVal = (bssrdf(albedo.getLuminance(), sigma_t.getLuminance(), d) + 1.0f)* (sample.LeDotG) * p * tr;
                    // printf("pHat: %f\n", retVal.getLuminance());
                    // printf("bssrdf: %f, di: %f, LeDotG: %f, p: %f, tr: %f\n",bssrdf(albedo.getLuminance(), sigma_t.getLuminance(), d),sample.di.getLuminance(), sample.LeDotG.getLuminance(), p, tr.getLuminance());
                    return retVal.getLuminance();
                    // return 1.0f;
                };

                auto ris_p = [scatteringJacobi] (const scatteringPointNNEE& sample, const Point& currentSamplePoint, bool requireJ = true) -> Float
                {
                    // printf("ris_p: %f\n", sample.scatteringPdf * scatteringJacobi(sample.sampleCentre, sample.scatteringPoint, currentSamplePoint));
                    if (requireJ)
                        return sample.scatteringPdf * scatteringJacobi(sample.sampleCentre, sample.scatteringPoint, currentSamplePoint);
                    else
                        return sample.scatteringPdf;
                    // return 1.0f;
                };

                auto fulfillMRec = [] (const scatteringPointNNEE& sample, MediumSamplingRecord& mRec) ->void
                {
                    // t p sigmaA sigmaS time medium transmittance
                    // printf("distance %f : %f\n", mRec.t, distance(sample.scatteringPoint, mRec.p));
                    mRec.t = distance(sample.scatteringPoint, mRec.p);
                    // printf("mRec: %f, [%f, %f, %f], [%f]\n", mRec.t, mRec.p.x, mRec.p.y, mRec.p.z, mRec.transmittance.getLuminance());
                    // printf("point %f : %f\n", mRec.p.x, sample.scatteringPoint.x);
                    mRec.p = sample.scatteringPoint;
                    // sigmaA & sigmaS & time & medium remains unchanged
                    Ray pRay;
                    pRay.mint = 0; pRay.maxt = mRec.t;

                    // printf("tr %f : %f\n", mRec.transmittance.getLuminance(), mRec.medium->evalTransmittance(pRay).getLuminance());
                    mRec.transmittance = mRec.medium->evalTransmittance(pRay);

                    // mRec.pdfFailure = 1.0f;
                    // mRec.pdfSuccess = 1.0f;
                    // mRec.pdfSuccessRev = 1.0f;
                };

                Li += nneeContribution;

                Float phasePdf = mSolidAnglePdfArray[0];
                Float phaseVal = mSolidAnglePdfArray[0];

                Float RISWeight = 1.0f;
                Vector rayDir;

                // explicit surface/medium sampling
                const Float estimateSurfaceP = sampleSurfacePdf(mRecordArray[0], rRec, scene);
                Float estimateMultipleScatteringP = 1.0f - estimateSurfaceP;
                bool estimateSurface = rRec.sampler->next1D() < estimateSurfaceP;

                if (!enableRIS || estimateSurface)
                {
                    mRec = mRecordArray[0];
                }
                else
                {
                    // start RIS
                    // printf("RIS begin!\n");
                    Point currentSamplePoint = mRec.p;

                    int sampleIndex = -1;
                    Float weightSum = 0.0f;
                    Float sampleCount = 0.0f;

                    bool select2=true;

                    // sample from proxy distributions
                    // printf("RIS with %d proxy samples!\n",selectedScatteringPoint.size());

                    // TODO: bug fix here. this implementation indroduces bias (or variance?)
                    // for (int i = selectedScatteringPoint.size() - 1; i<selectedScatteringPoint.size();i++)
                    // {
                    //     scatteringPointNNEE& item = *selectedScatteringPoint[i];
                    //     if(item.selected || item.scatteringPoint == mRec.p) continue;
                    //     bool requireJ = (i != selectedScatteringPoint.size() - 1);
                    //     Float weight = pHat(item, currentSamplePoint,ray.d) / ris_p(item, currentSamplePoint, requireJ);
                    //     weightSum += weight;
                    //     sampleCount+=1.0f;

                    //     // if (weight > 200.0f) printf("weight: %f\n",weight);
                    //     if (rRec.sampler->next1D() < weight / weightSum) 
                    //     {
                    //         select2 = false;
                    //         sampleIndex = i;
                    //     }
                    // }

                    // sample from casual distributions
                    // printf("RIS with %d casual samples!\n",previousCasualScatter.size());
                    for (int i = 0; i<previousCasualScatter.size();i++)
                    {
                        scatteringPointNNEE& item = previousCasualScatter[i];
                        if(item.selected || item.scatteringPoint == mRec.p) continue;
                        Float weight = pHat(item, currentSamplePoint,ray.d) / ris_p(item, currentSamplePoint);
                        weightSum += weight;
                        sampleCount+=1.0f;
                        if (rRec.sampler->next1D() <= weight / weightSum) 
                        {
                            select2 = true;
                            sampleIndex = i;
                        }
                    }

                    // printf("select2: %d\n", select2);

                    if (sampleIndex == -1) return Li;

                    // set mRec, phaseVal and RISWeight. phasePdf is not important.
                    if(!select2){
                        RISWeight = weightSum / (pHat(*selectedScatteringPoint[sampleIndex],currentSamplePoint, ray.d) * sampleCount);
                        selectedScatteringPoint[sampleIndex]->selected = true;
                        rayDir = normalize(selectedScatteringPoint[sampleIndex]->scatteringPoint - mRec.p);
                    }
                    else {
                        RISWeight = weightSum / (pHat(previousCasualScatter[sampleIndex],currentSamplePoint, ray.d) * sampleCount);
                        previousCasualScatter[sampleIndex].selected = true;
                        rayDir = normalize(previousCasualScatter[sampleIndex].scatteringPoint - mRec.p);
                    }
                    // printf("set done!\n");
                    // printf("RISWeight: %f\n", RISWeight);
                    phaseVal = mRec.getPhaseFunction()->eval(-ray.d, rayDir);
                    phasePdf = phaseVal;

                    Point originalmRecp = mRec.p;
                    mRec = mRecordArray[0];
                    mRec.p = originalmRecp;
                    if(!select2)
                        fulfillMRec(*selectedScatteringPoint[sampleIndex], mRec);
                    else
                        fulfillMRec(previousCasualScatter[sampleIndex], mRec);

                    // printf("RIS choose next point: %f\n", mRec.p.x);

                    // terminate for corner case (can't find any samples)
                    if (RISWeight == 0.0f) return Li;
                }



                // p_hat = bssrdf * L_e * G

                

                /* ==================================================================== */
                /*                         Multiple scattering                          */
                /* ==================================================================== */

                /* Stop if multiple scattering was not requested */
                if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
                    break;
                rRec.type = RadianceQueryRecord::ERadianceNoEmission;

                if (estimateSurface){
                    // estimate ONLY surface term
                    // throughput *= (phase * transmittance) / (phasePdf * estimateSurfacePdf)
                    // direction is sampled according to phase function
                    // printf("leave medium\n");
                    its = mItsArray[0];
                    ray = Ray(its.p, pRecordArray[0].wo, ray.time);
                    ray.mint = Epsilon;
                    Spectrum transmittance = rRec.medium->evalTransmittance(Ray(ray, 0, its.t), rRec.sampler);
                    throughput *= phaseVal / (phasePdf * estimateSurfaceP);
                    

                    if (!its.isValid()) {
                        if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered)) {
                            Spectrum value = throughput * scene->evalEnvironment(ray);
                            if (rRec.medium)
                                value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                            Li += value;
                        }
                        break;
                    }

                    throughput *= transmittance;
                    rRec.medium = its.getTargetMedium(ray.d);
                    scene->rayIntersect(ray, its);
                    rRec.depth++;

                }
                else {
                    if (!its.isValid()) break; // no problem

                    if (!enableRIS)
                    {
                        ray = Ray(mRec.p, pRecordArray[0].wo, ray.time);
                        ray.mint = Epsilon;
                        throughput *= phaseVal * mRec.sigmaS * mRec.transmittance / (phasePdf * mRec.pdfSuccess * estimateMultipleScatteringP);
                        its = mItsArray[0];
                    }
                    else
                    {
                        ray = Ray(mRec.p, rayDir, ray.time);
                        ray.mint = Epsilon;
                        throughput *= phaseVal * mRec.sigmaS * mRec.transmittance * RISWeight / estimateMultipleScatteringP;
                    }

                }

            } 
            else {

                if (!its.isValid()) {
                    /* If no intersection could be found, possibly return
                       attenuated radiance from a background luminaire */
                    if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                        && (!m_hideEmitters || scattered)) {
                        Spectrum value = throughput * scene->evalEnvironment(ray);
                        if (rRec.medium)
                            value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                        Li += value;
                    }

                    break;
                }

                /* Possibly include emitted radiance if requested */
                if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * its.Le(-ray.d);

                /* Include radiance from a subsurface integrator if requested */
                if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                    Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

                if (rRec.depth >= m_maxDepth && m_maxDepth != -1)
                    break;

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
                      wiDotShN  = Frame::cosTheta(its.wi);
                if (wiDotGeoN * wiDotShN < 0 && m_strictNormals)
                    break;

                /* ==================================================================== */
                /*                          Luminaire sampling                          */
                /* ==================================================================== */

                const BSDF *bsdf = its.getBSDF(ray);
                DirectSamplingRecord dRec(its);

                /* Estimate the direct illumination if this is requested */
                if ((rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
                    (bsdf->getType() & BSDF::ESmooth)) {
                    int interactions = m_maxDepth - rRec.depth - 1;

                    Spectrum value = scene->sampleAttenuatedEmitterDirect(
                            dRec, its, rRec.medium, interactions,
                            rRec.nextSample2D(), rRec.sampler);

                    if (!value.isZero()) {
                        const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                        /* Evaluate BSDF * cos(theta) */
                        BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                        const Spectrum bsdfVal = bsdf->eval(bRec);

                        Float woDotGeoN = dot(its.geoFrame.n, dRec.d);

                        /* Prevent light leaks due to the use of shading normals */
                        if (!bsdfVal.isZero() && (!m_strictNormals ||
                            woDotGeoN * Frame::cosTheta(bRec.wo) > 0)) {
                            /* Calculate prob. of having generated that direction
                               using BSDF sampling */
                            Float bsdfPdf = (emitter->isOnSurface()
                                    && dRec.measure == ESolidAngle)
                                    ? bsdf->pdf(bRec) : (Float) 0.0f;

                            /* Weight using the power heuristic */
                            const Float weight = miWeight(dRec.pdf, bsdfPdf);
                            Li += throughput * value * bsdfVal * weight;
                        }
                    }
                }


                /* ==================================================================== */
                /*                            BSDF sampling                             */
                /* ==================================================================== */

                /* Sample BSDF * cos(theta) */
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                Float bsdfPdf;
                Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
                if (bsdfWeight.isZero())
                    break;

                /* Prevent light leaks due to the use of shading normals */
                const Vector wo = its.toWorld(bRec.wo);
                Float woDotGeoN = dot(its.geoFrame.n, wo);
                if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
                    break;

                /* Trace a ray in this direction */
                ray = Ray(its.p, wo, ray.time);

                /* Keep track of the throughput, medium, and relative
                   refractive index along the path */
                throughput *= bsdfWeight;
                eta *= bRec.eta;
                if (its.isMediumTransition())
                    rRec.medium = its.getTargetMedium(ray.d);

                /* Handle index-matched medium transitions specially */
                if (bRec.sampledType == BSDF::ENull) {
                    if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                        break;
                    rRec.type = scattered ? RadianceQueryRecord::ERadianceNoEmission
                        : RadianceQueryRecord::ERadiance;
                    scene->rayIntersect(ray, its);
                    rRec.depth++;

                    // OUT_DEBUG("null intersection");
                    // corner case for null BSDF
                    if (rRec.medium){
                        auto insideMedium = rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler);
                        if (insideMedium) {
                            throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;
                            enterMedium = true;
                            ray = Ray(mRec.p, ray.d, ray.time);
                            ray.mint = Epsilon;
                            // scene->rayIntersect(ray, its);
                            continue;
                        }
                        else {
                            if (!its.isValid()) {
                                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered)) {
                                    Spectrum value = throughput * scene->evalEnvironment(ray);
                                    if (rRec.medium)
                                        value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                                    Li += value;
                                }
                                break;
                            }
                            throughput *= mRec.transmittance / mRec.pdfFailure;
                            rRec.medium = its.getTargetMedium(ray.d);
                            // printf("rRec.medium == nullptr? %d\n",rRec.medium == nullptr);
                            ray = Ray(its.p, ray.d, ray.time);
                            ray.mint = Epsilon;
                            scene->rayIntersect(ray, its);
                            rRec.depth++;
                        }
                    }
                }

                Spectrum value(0.0f);
                rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
                    m_maxDepth - rRec.depth - 1, ray, its, dRec, value);

                /* If a luminaire was hit, estimate the local illumination and
                   weight using the power heuristic */
                if (!value.isZero() && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                    const Float emitterPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                        scene->pdfEmitterDirect(dRec) : 0;
                    Li += throughput * value * miWeight(bsdfPdf, emitterPdf);
                }

                /* ==================================================================== */
                /*                         Indirect illumination                        */
                /* ==================================================================== */

                /* Stop if indirect illumination was not requested */
                if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                    break;

                rRec.type = RadianceQueryRecord::ERadianceNoEmission;

                // enter medium
                if (rRec.medium){
                    auto insideMedium = rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler);
                    if (insideMedium) {
                        throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;
                        enterMedium = true;
                        ray = Ray(mRec.p, ray.d, ray.time);
                        ray.mint = Epsilon;
                        // scene->rayIntersect(ray, its);
                        continue;
                    }
                    else {
                        if (!its.isValid()) {
                            if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered)) {
                                Spectrum value = throughput * scene->evalEnvironment(ray);
                                if (rRec.medium)
                                    value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                                Li += value;
                            }
                            break;
                        }
                        throughput *= mRec.transmittance / mRec.pdfFailure;
                        rRec.medium = its.getTargetMedium(ray.d);
                        // printf("rRec.medium == nullptr? %d\n",rRec.medium == nullptr);
                        ray = Ray(its.p, ray.d, ray.time);
                        ray.mint = Epsilon;
                        scene->rayIntersect(ray, its);
                        rRec.depth++;
                    }
                }
            }

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }

            scattered = true;
        }
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;
        return Li;
    }

    /**
     * This function is called by the recursive ray tracing above after
     * having sampled a direction from a BSDF/phase function. Due to the
     * way in which this integrator deals with index-matched boundaries,
     * it is necessarily a bit complicated (though the improved performance
     * easily pays for the extra effort).
     *
     * This function
     *
     * 1. Intersects 'ray' against the scene geometry and returns the
     *    *first* intersection via the '_its' argument. (note: *last* intersection in fact)
     *
     * 2. It checks whether the intersected shape was an emitter, or if
     *    the ray intersects nothing and there is an environment emitter.
     *    In this case, it returns the attenuated emittance, as well as
     *    a DirectSamplingRecord that can be used to query the hypothetical
     *    sampling density at the emitter.
     *
     * 3. If current shape is an index-matched medium transition, the
     *    integrator keeps on looking on whether a light source eventually
     *    follows after a potential chain of index-matched medium transitions,
     *    while respecting the specified 'maxDepth' limits. It then returns
     *    the attenuated emittance of this light source, while accounting for
     *    all attenuation that occurs on the wya.
     */
    void rayIntersectAndLookForEmitter(const Scene *scene, Sampler *sampler,
            const Medium *medium, int maxInteractions, Ray ray, Intersection &_its,
            DirectSamplingRecord &dRec, Spectrum &value, Intersection *firstIts = nullptr) const {
        Intersection its2, *its = &_its;
        Spectrum transmittance(1.0f);
        bool surface = false;
        int interactions = 0;

        bool isFirstIts = true;

        while (true) {
            surface = scene->rayIntersect(ray, *its);

            if (firstIts != nullptr && isFirstIts) *firstIts = *its;
            isFirstIts = false;

            if (medium)
                transmittance *= medium->evalTransmittance(Ray(ray, 0, its->t), sampler);

            if (surface && (interactions == maxInteractions ||
                !(its->getBSDF()->getType() & BSDF::ENull) ||
                its->isEmitter())) {
                /* Encountered an occluder / light source */
                break;
            }

            if (!surface)
                break;

            if (transmittance.isZero())
                return;

            if (its->isMediumTransition())
                medium = its->getTargetMedium(ray.d);

            Vector wo = its->shFrame.toLocal(ray.d);
            BSDFSamplingRecord bRec(*its, -wo, wo, ERadiance);
            bRec.typeMask = BSDF::ENull;
            transmittance *= its->getBSDF()->eval(bRec, EDiscrete);

            ray.o = ray(its->t);
            ray.mint = Epsilon;
            its = &its2;

            if (++interactions > 100) { /// Just a precaution..
                Log(EWarn, "rayIntersectAndLookForEmitter(): round-off error issues?");
                return;
            }
        }

        if (surface) {
            /* Intersected something - check if it was a luminaire */
            if (its->isEmitter()) {
                dRec.setQuery(ray, *its);
                value = transmittance * its->Le(-ray.d);
            }
        } else {
            /* Intersected nothing -- perhaps there is an environment map? */
            const Emitter *env = scene->getEnvironmentEmitter();

            if (env && env->fillDirectSamplingRecord(dRec, ray))
                value = transmittance * env->evalEnvironment(RayDifferential(ray));
        }
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "VolumetricPathTracer[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(VolumetricPathTracerNNEE, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricPathTracerNNEE, "Volumetric path tracer with Next-Next Event Estimation (NNEE)");
MTS_NAMESPACE_END
