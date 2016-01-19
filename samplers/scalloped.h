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


class IRegion {
public:
    virtual ~IRegion() {}

    size_t dart_id() const { return dart_id_; }

    virtual bool IsOverlap(const Circle &circle) const = 0;

    virtual vector<IRegion*> Eclipse(const Circle &circle) const = 0;

    virtual Point2D ThrowDart() const = 0;

protected:
    IRegion(size_t dart_id) : dart_id_(dart_id) {}

    size_t dart_id_;
};


class RegionFactor {
public:
    virtual ~RegionFactor() {}

    virtual vector<IRegion*> CreateRegions(const Point2D &center) = 0;

    static RegionFactor* CreateRegionFactor(float r_ratio);

protected:
    RegionFactor() {}
};


class IRegionSelector {
public:
    virtual ~IRegionSelector() {}

    virtual const IRegion* Select() const = 0;

    virtual void Add(const IRegion *region) const = 0;

    virtual void Reomve(const IRegion *region) const = 0;

protected:
    IRegionSelector() {}
};


class RegionSelectorFactor {
public:
    virtual ~RegionSelectorFactor() {}

    virtual IRegionSelector* CreateRegionSelector() = 0;

    static RegionSelectorFactor* CreateRegionSelectorFactor(bool weighted);

protected:
    RegionSelectorFactor() {}
};


struct Dart {
    Dart(const Point2D &position) : position(position) {}

    Point2D position;

    std::unordered_set<IRegion*> regions;
};


class IDartsNet {
public:
    virtual ~IDartsNet() {}

    virtual vector<size_t> GetNeighbors(const Point2D &position) const = 0;

    virtual void AddNeighbor(size_t mentor_id, size_t neighbor_id) const = 0;

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
