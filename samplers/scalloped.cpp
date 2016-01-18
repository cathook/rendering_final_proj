#include "samplers/scalloped.h"

#include <stdint.h>

#include <memory>

#include "montecarlo.h"
#include "pool.h"


namespace {


class Tmp : public ISquareSampler {
public:
    vector<Point2D> Sample(size_t num_samples, RNG &rng) const {
        vector<Point2D> ret(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            ret[i].x = rng.RandomFloat();
            ret[i].y = rng.RandomFloat();
        }
        return ret;
    }
};


class ScallopedSampler : public Sampler {
public:
    ScallopedSampler(int x_start, int x_end, int y_start, int y_end,
                     int pixel_samples,
                     float s_open, float s_close,
                     std::shared_ptr<const ISquareSampler> sampler) :
            Sampler(x_start, x_end, y_start, y_end,
                    pixel_samples,
                    s_open, s_close),
            x_pos_(x_start),
            y_pos_(y_start),
            time_samples_buf_(new float[pixel_samples]),
            sampler_(sampler) {}

    Sampler *GetSubSampler(int num, int count) {
        int x0, x1, y0, y1;
        ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
        if (x0 == x1 || y0 == y1) {
            return NULL;
        }
        return new ScallopedSampler(x0, x1, y0, y1,
                                    samplesPerPixel,
                                    shutterOpen, shutterClose,
                                    sampler_);
    }

    int RoundSize(int size) const { return size; }

    int MaximumSampleCount() { return samplesPerPixel; }

    int GetMoreSamples(Sample *samples, RNG &rng) {
        if (y_pos_ == yPixelEnd) {
            return 0;
        }

        vector<Point2D> img_xy(sampler_->Sample(samplesPerPixel, rng));
        vector<Point2D> lens_uv(sampler_->Sample(samplesPerPixel, rng));
        StratifiedSample1D(time_samples_buf_.get(), samplesPerPixel, rng, true);

        Shuffle(&lens_uv[0], samplesPerPixel, 1, rng);
        Shuffle(time_samples_buf_.get(), samplesPerPixel, 1, rng);

        for (int i = 0; i < samplesPerPixel; ++i) {
            samples[i].imageX = x_pos_ + img_xy[i].x;
            samples[i].imageY = y_pos_ + img_xy[i].y;
            samples[i].lensU = lens_uv[i].x;
            samples[i].lensV = lens_uv[i].y;
            samples[i].time = time_samples_buf_[i];

            for (uint32_t j = 0; j < samples[i].n2D.size(); ++j) {
                vector<Point2D> tmp(sampler_->Sample(samples[i].n2D[j], rng));
                for (size_t k = 0; k < tmp.size(); ++k) {
                    samples[i].twoD[j][k * 2 + 0] = tmp[k].x;
                    samples[i].twoD[j][k * 2 + 1] = tmp[k].y;
                }
            }

            for (uint32_t j = 0; j < samples[i].n1D.size(); ++j) {
                LatinHypercube(samples[i].oneD[j], samples[i].n1D[j], 1, rng);
            }
        }

        ++x_pos_;
        if (x_pos_ == xPixelEnd) {
            x_pos_ = xPixelStart;
            ++y_pos_;
        }

        return samplesPerPixel;
    }

private:
    int x_pos_;
    int y_pos_;

    std::unique_ptr<float[]> time_samples_buf_;

    std::shared_ptr<const ISquareSampler> sampler_;
};


}  // namespace


Sampler *CreateScallopedSampler(
        const ParamSet &params, const Film *film, const Camera *camera) {
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);

    int pixel_samples = params.FindOneInt("pixelsamples", 1);
    Assert(pixel_samples > 0);

    return new ScallopedSampler(xstart, xend, ystart, yend,
                                pixel_samples,
                                camera->shutterOpen, camera->shutterClose,
                                std::shared_ptr<const ISquareSampler>(new Tmp()));
}
