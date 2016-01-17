#ifndef PBRT_SAMPLERS_SCALLOPED_H_
#define PBRT_SAMPLERS_SCALLOPED_H_


#if defined(_MSC_VER)
#pragma once
#endif


#include "camera.h"
#include "film.h"
#include "paramset.h"
#include "sampler.h"


Sampler *CreateScallopedSampler(
        const ParamSet &params, const Film *film, const Camera *camera);

#endif  // PBRT_SAMPLERS_SCALLOPED_H_
