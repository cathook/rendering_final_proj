#ifndef PBRT_SAMPLERS_SCALLOPED_H_
#define PBRT_SAMPLERS_SCALLOPED_H_


#if defined(_MSC_VER)
#pragma once
#endif


#include <unordered_set>

#include "camera.h"
#include "film.h"
#include "paramset.h"
#include "sampler.h"


class IRegion;


struct Dart {
    Dart(const Point2D &position) : position(position) {}

    const Point2D position;

    std::unordered_set<IRegion*> regions;
};


class IRegion {
public:
    virtual ~IRegion() {}

    Dart *dart() const { return dart_; }

    virtual float area() const = 0;

    virtual bool Eclipse(const Circle &circle, vector<IRegion*> *out) const = 0;

    virtual Point2D SelectPoint(const RNG &rng) const = 0;

protected:
    IRegion(Dart *dart) : dart_(dart) {}

    Dart *dart_;
};


class RegionFactory {
public:
    virtual ~RegionFactory() {}

    virtual vector<IRegion*> CreateRegions(Dart *dart) const = 0;

    static RegionFactory* CreateRegionFactory(float r_ratio);

protected:
    RegionFactory() {}
};


class IRegionSelector {
public:
    virtual ~IRegionSelector() {}

    virtual IRegion* Select() const = 0;

    virtual void Add(IRegion *region) const = 0;

    virtual void Reomve(IRegion *region) const = 0;

protected:
    IRegionSelector() {}
};


class RegionSelectorFactory {
public:
    virtual ~RegionSelectorFactory() {}

    virtual IRegionSelector* CreateRegionSelector() = 0;

    static RegionSelectorFactory* CreateRegionSelectorFactory(bool weighted);

protected:
    RegionSelectorFactory() {}
};


struct DartNode : Dart {
    DartNode(const Point2D &position) : Dart(position) {}

    vector<DartNode*> neighbors;
};


class IDartsNet {
public:
    virtual ~IDartsNet() {}

    virtual vector<DartNode*> GetNeighbors(
            const DartNode *mentor, const Point2D &position) const = 0;

    virtual void AddNeighbor(DartNode *mentor, DartNode* neighbor) const = 0;

protected:
    IDartsNet() {}
};


class ISquareSampler {
public:
    virtual ~ISquareSampler() {}

    virtual void Sample(size_t num_samples, RNG &rng,
                        vector<Point2D> *out) const = 0;

protected:
    ISquareSampler() {}
};


Sampler *CreateScallopedSampler(
        const ParamSet &params, const Film *film, const Camera *camera);

#endif  // PBRT_SAMPLERS_SCALLOPED_H_
