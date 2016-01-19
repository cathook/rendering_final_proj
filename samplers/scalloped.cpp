#include "samplers/scalloped.h"

#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>

#include "montecarlo.h"
#include "rng.h"
#include "tree_style_array.h"


using std::lower_bound;
using std::pair;
using std::sort;
using std::unique;
using std::unordered_map;


namespace {


class ArcRegion : public IRegion {
public:
    ArcRegion(Dart *dart, float theta0, float theta) :
            IRegion(dart), theta0_(theta0), theta_(theta) {}

    float area() const { return theta_; }

    bool Eclipse(const Circle &circle, vector<IRegion*> *out) const {
        Vector2D c12(circle.center - dart()->position);
        float l2 = c12.LengthSquared();
        float l = sqrtf(l2);
        if (l + 1.f <= circle.radius) {
            out->clear();
            return true;
        }
        float inv_l2 = 1.f / l2;
        float r12 = inv_l2;
        float r22 = circle.radius * circle.radius * inv_l2;

        float theta = acosf((1.f + r12 - r22) * 0.5f * l);
        if (isnan(theta) || theta == 0.f) {
            return false;
        }

        float phi = atan2f(c12.y, c12.x);

        Explode_(phi - theta, phi + theta, out);
        return true;
    }

    Point2D SelectPoint(const RNG &rng) const {
        float angle = theta0_ + rng.RandomFloat() * theta_;
        return dart()->position + Vector2D(1.f, 0.f).Rotate(angle);
    }

private:
    void Explode_(float t1, float t2, vector<IRegion*> *out) const {
        t1 = AdjustToRange_(t1);
        t2 = AdjustToRange_(t2);

        out->clear();
        if (t2 <= t1) {
            swap(t1, t2);
            float tlen = min(t2 - t1, theta0_ + theta_ - t1);
            if (tlen > 0.f) {
                out->emplace_back(
                        new ArcRegion(dart(), AdjustAngle_(t1), tlen));
            }
        } else {
            if (t1 != theta0_) {
                out->emplace_back(new ArcRegion(dart(), theta0_, t1 - theta0_));
            }
            float t3 = theta0_ + theta_;
            if (t2 < t3) {
                out->emplace_back(
                        new ArcRegion(dart(), AdjustAngle_(t2), t3 - t2));
            }
        }
    }

    static float AdjustAngle_(float angle) {
        while (angle < 0.f) {
            angle += 2.f * M_PI;
        }
        while (angle > 2.f * M_PI) {
            angle -= 2.f * M_PI;
        }
        return angle;
    }

    float AdjustToRange_(float angle) const {
        while (angle - theta0_ > theta_) {
            angle -= 2.f * M_PI;
        }
        while (angle < theta0_) {
            angle += 2.f * M_PI;
        }
        return angle;
    }

    float theta0_;
    float theta_;
};


class ArcRegionFactory : public RegionFactory {
public:
    vector<IRegion*> CreateRegions(Dart *dart) const {
        return vector<IRegion*>(1, new ArcRegion(dart, 0, 2.f * M_PI));
    }
};


class TreeStyleArraySquareSampler : public ISquareSampler {
public:
    TreeStyleArraySquareSampler(float size, vector<Point2D> points) :
            size_(size) {
        // Discretize the Y values.
        for (Point2D &p : points) {
            ys_.push_back(p.y);
        }
        sort(ys_.begin(), ys_.end());
        ys_.erase(unique(ys_.begin(), ys_.end()), ys_.end());
        unordered_map<float, size_t> y_index;
        for (size_t i = 0; i < ys_.size(); ++i) {
            y_index[ys_[i]] = i;
        }

        // Discretize the X values.
        vector<pair<size_t, size_t>> d_points;
        d_points.reserve(points.size());
        sort(points.begin(), points.end(), ComparePoint2Ds);
        float x_prev = points[0].x - 1;
        for (Point2D &p : points) {
            if (x_prev != p.x) {
                xs_.emplace_back(p.x);
                x_prev = p.x;
            }
            d_points.emplace_back(xs_.size() - 1, y_index[p.y]);
        }

        // Setups the table and xys index map.
        xys_.resize(xs_.size());
        table_ = new TreeStyleArray2D(xs_.size(), ys_.size());
        for (pair<size_t, size_t> &p : d_points) {
            table_->Update(p.first + 1, p.second + 1, 1);
            xys_[p.first].emplace_back(p.second);
        }

        if (xs_.back() != size_) {
            xs_.push_back(size_);
        }
        if (ys_.back() != size_) {
            ys_.push_back(size_);
        }
    }

    ~TreeStyleArraySquareSampler() {
        delete table_;
    }

    void Sample(size_t num_samples, RNG &rng, vector<Point2D> *out) const {
        Assert(num_samples > 0);

        auto it = xy_max_.insert(pair<size_t, float>(num_samples, size_)).first;
        while (true) {
            // Anchors the bottom-left point and its xy indexes.
            Point2D p0(rng.RandomFloat() * it->second,
                       rng.RandomFloat() * it->second);
            size_t ix0 = BinarySearchIndex(xs_, 0, p0.x);
            size_t iy0 = BinarySearchIndex(ys_, 0, p0.y);

            // Binary searches the suitable square size.
            float lower = 0, upper = min(size_ - p0.x, size_ - p0.y);
            for (int i = 0; i < 50; ++i) {
                // Anchors the top-right xy indexes.
                float mid = (lower + upper) * 0.5;
                size_t ix1 = BinarySearchIndex(xs_, ix0, p0.x + mid);
                size_t iy1 = BinarySearchIndex(ys_, iy0, p0.y + mid);

                size_t n = GetNumPoints(ix0, ix1, iy0, iy1);

                // Updates.
                if (n < num_samples) {
                    lower = mid;
                } else if (n > num_samples) {
                    upper = mid;
                } else {
                    // A solution is found.
                    float size = RandomSquareSize(p0, ix1, iy1, rng);
                    SampleInSquare(ix0, ix1, iy0, iy1, p0, size, out);
                    return;
                }
            }

            it->second = min(p0.x, p0.y);
        }
    }

private:
    static bool ComparePoint2Ds(const Point2D &a, const Point2D &b) {
        return (a.x != b.x ? a.x < b.x : a.y < b.y);
    }

    static size_t BinarySearchIndex(
            const vector<float> &arr, size_t index0, float value) {
        auto it = lower_bound(arr.begin() + index0, arr.end(), value);
        return it - arr.begin();
    }

    size_t GetNumPoints(size_t ix_begin, size_t ix_end,
                        size_t iy_begin, size_t iy_end) const {
        return static_cast<size_t>(table_->Query(ix_end, iy_end)
                                   - table_->Query(ix_end, iy_begin)
                                   - table_->Query(ix_begin, iy_end)
                                   + table_->Query(ix_begin, iy_begin));
    }

    float RandomSquareSize(const Point2D &p0,
                           size_t ix_end, size_t iy_end, RNG &rng) const {
        float size_max = min(xs_[ix_end] - p0.x, ys_[iy_end] - p0.y);
        float size_min = max(xs_[ix_end - 1] - p0.x, ys_[iy_end - 1] - p0.y);
        return size_min + rng.RandomFloat() * (size_max - size_min);
    }

    void SampleInSquare(size_t ix_begin, size_t ix_end,
                        size_t iy_begin, size_t iy_end,
                        const Point2D &p0, float size,
                        vector<Point2D> *out) const {
        float inv_size = 1.f / size;
        size_t counter = 0;
        for (size_t ix = ix_begin; ix < ix_end; ++ix) {
            for (size_t iy : xys_[ix]) {
                if (iy_begin <= iy && iy < iy_end) {
                    Point2D p(xs_[ix], ys_[iy]);
                    out->at(counter) = Point2D(p - p0) * inv_size;
                    ++counter;
                }
            }
        }
    }

    float size_;
    vector<float> xs_;
    vector<float> ys_;
    vector<vector<size_t>> xys_;
    mutable unordered_map<size_t, float> xy_max_;
    TreeStyleArray2D *table_;
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
            sampler_(sampler),
            img_xy_buf_(samplesPerPixel),
            lens_uv_buf_(samplesPerPixel),
            time_samples_buf_(new float[pixel_samples]) {}

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

        sampler_->Sample(samplesPerPixel, rng, &img_xy_buf_);
        sampler_->Sample(samplesPerPixel, rng, &lens_uv_buf_);
        StratifiedSample1D(time_samples_buf_.get(), samplesPerPixel, rng, true);

        Shuffle(&lens_uv_buf_[0], samplesPerPixel, 1, rng);
        Shuffle(time_samples_buf_.get(), samplesPerPixel, 1, rng);

        for (int i = 0; i < samplesPerPixel; ++i) {
            samples[i].imageX = x_pos_ + img_xy_buf_[i].x;
            samples[i].imageY = y_pos_ + img_xy_buf_[i].y;
            samples[i].lensU = lens_uv_buf_[i].x;
            samples[i].lensV = lens_uv_buf_[i].y;
            samples[i].time = time_samples_buf_[i];

            for (uint32_t j = 0; j < samples[i].n2D.size(); ++j) {
                extra_buf_.resize(samples[i].n2D[j]);
                sampler_->Sample(samples[i].n2D[j], rng, &extra_buf_);
                for (size_t k = 0; k < extra_buf_.size(); ++k) {
                    samples[i].twoD[j][k * 2 + 0] = extra_buf_[k].x;
                    samples[i].twoD[j][k * 2 + 1] = extra_buf_[k].y;
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

    std::shared_ptr<const ISquareSampler> sampler_;

    vector<Point2D> img_xy_buf_;
    vector<Point2D> lens_uv_buf_;
    std::unique_ptr<float[]> time_samples_buf_;
    vector<Point2D> extra_buf_;
};


}  // namespace


RegionFactory* RegionFactory::CreateRegionFactory(float r_ratio) {
    if (r_ratio == 1.f) {
        return new ArcRegionFactory();
    }
    Assert(false);
    return NULL;
}


Sampler *CreateScallopedSampler(
        const ParamSet &params, const Film *film, const Camera *camera) {
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);

    float r_ratio = params.FindOneFloat("rratio", 1.f);
    RegionFactory *region_factory = RegionFactory::CreateRegionFactory(r_ratio);
    Assert(region_factory != NULL);

    int pixel_samples = params.FindOneInt("pixelsamples", 1);
    Assert(pixel_samples > 0);

    delete region_factory;
    ///////
    //
    RNG rng;
    vector<Point2D> points;
    for (int i = 0; i < 1000; ++i) {
        points.emplace_back(rng.RandomFloat() * 100, rng.RandomFloat() * 100);
    }
    std::shared_ptr<const ISquareSampler> square_sampler(
            new TreeStyleArraySquareSampler(100, points));

    return new ScallopedSampler(xstart, xend, ystart, yend,
                                pixel_samples,
                                camera->shutterOpen, camera->shutterClose,
                                square_sampler);
}
