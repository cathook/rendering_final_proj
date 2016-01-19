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

    virtual ~Dart() {}

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

    virtual IRegion* Select(const RNG &rng) const = 0;

    virtual void Add(IRegion *region) = 0;

    virtual void Remove(IRegion *region) = 0;

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


class IDartsNet {
public:
    virtual ~IDartsNet() {}

    virtual Dart* WrapDart(const Dart &dart) const = 0;

    virtual void DestroyAll(Dart *master) = 0;

    virtual size_t GetNears(Dart *master, const Point2D &position,
                            vector<Dart*> *out) const = 0;

    virtual void AddNeighbor(Dart *master, Dart *neighbor) = 0;

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
