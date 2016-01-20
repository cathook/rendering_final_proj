#include "samplers/scalloped.h"

#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "montecarlo.h"
#include "rng.h"
#include "tree_style_array.h"


using std::lower_bound;
using std::pair;
using std::queue;
using std::sort;
using std::unique;
using std::unordered_map;
using std::unordered_set;


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

        return Explode_(phi - theta, phi + theta, out);
    }

    bool Eclipse(const Line2D &line, vector<IRegion*> *out) const {
        Vector2D v(~-Normalize(line.v));

        float theta = acosf((line.p0 - dart()->position).Dot(v));
        if (isnan(theta) || theta == 0.f) {
            return false;
        }

        float phi = atan2f(v.y, v.x);

        return Explode_(phi - theta, phi + theta, out);
    }

    Point2D SelectPoint(const RNG &rng) const {
        float angle = theta0_ + rng.RandomFloat() * theta_;
        return dart()->position + Vector2D(1.f, 0.f).Rotate(angle);
    }

private:
    bool Explode_(float t1, float t2, vector<IRegion*> *out) const {
        t1 = AdjustToRange_(t1);
        t2 = AdjustToRange_(t2);
        if (t2 <= t1) {
            swap(t1, t2);
            if (t1 == theta0_) {
                return false;
            }
            out->clear();
            float tlen = min(t2 - t1, theta0_ + theta_ - t1);
            if (tlen > 0.f) {
                out->emplace_back(
                        new ArcRegion(dart(), AdjustAngle_(t1), tlen));
            }
        } else {
            float t3 = theta0_ + theta_;
            if (t3 <= t1) {
                return false;
            }
            out->clear();
            if (t1 != theta0_) {
                out->emplace_back(new ArcRegion(dart(), theta0_, t1 - theta0_));
            }
            if (t2 < t3) {
                out->emplace_back(
                        new ArcRegion(dart(), AdjustAngle_(t2), t3 - t2));
            }
        }
        return true;
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


class UniformRegionSelector : public IRegionSelector {
public:
    IRegion* Select(const RNG &rng) const {
        size_t index = rng.RandomFloat() * regions_.size();
        return regions_[index];
    }

    void Add(IRegion *region) {
        indexes_[region] = regions_.size();
        regions_.emplace_back(region);
    }

    void Remove(IRegion *region) {
        size_t index = indexes_[region];
        regions_[index] = regions_.back();
        indexes_[regions_.back()] = index;
        regions_.pop_back();
        indexes_.erase(region);
    }

private:
    vector<IRegion*> regions_;
    unordered_map<IRegion*, size_t> indexes_;
};


class UniformRegionSelectorFactory : public RegionSelectorFactory {
public:
    virtual IRegionSelector* CreateRegionSelector() {
        return new UniformRegionSelector();
    }
};


class WeightedRegionSelector : public IRegionSelector {
public:
    WeightedRegionSelector() : weight_sum_(1, 0.f), base_(0), width_(1) {}

    IRegion* Select(const RNG &rng) const {
        float weight = rng.RandomFloat() * weight_sum_[0];
        size_t index = 0;
        while (index < base_) {
            float left_weight = weight_sum_[index * 2 + 1];
            if (weight_sum_[index * 2 + 2] == 0.f || weight <= left_weight) {
                index = index * 2 + 1;
            } else {
                weight -= left_weight;
                index = index * 2 + 2;
            }
        }
        return regions_[index - base_];
    }

    void Add(IRegion *region) {
        if (regions_.size() + 1 > width_) {
            weight_sum_.resize(weight_sum_.size() + width_ * 2, 0);
            for (size_t b = base_ * 2 + 1, b0 = base_; b > 0; b = b0, b0 /= 2) {
                for (size_t i = b0, j = b; i < b; ++i, ++j) {
                    weight_sum_[j] = weight_sum_[i];
                }
                for (size_t i = b0 + (b - b0 + 1) / 2; i < b; ++i) {
                    weight_sum_[i] = 0;
                }
            }
            base_ = base_ * 2 + 1;
            width_ *= 2;
        }

        size_t index = regions_.size();
        indexes_[region] = index;
        regions_.emplace_back(region);
        Update_(base_ + index, region->area());
    }

    void Remove(IRegion *region) {
        Update_(base_ + regions_.size() - 1, -regions_.back()->area());

        size_t index = indexes_[region];
        if (index + 1 != regions_.size()) {
            Update_(base_ + index, -region->area() + regions_.back()->area());
            regions_[index] = regions_.back();
            indexes_[regions_.back()] = index;
        }
        indexes_.erase(region);
        regions_.pop_back();

        if (regions_.size() < width_ / 4) {
            for (size_t b0 = 0, b = 1; b0 < base_; b0 = b, b = b * 2 + 1) {
                for (size_t i = b0, j = b; i < b; ++i, ++j) {
                    weight_sum_[i] = weight_sum_[j];
                }
            }
            weight_sum_.resize(weight_sum_.size() - width_);
            base_ /= 2;
            width_ /= 2;
        }
    }

private:
    void Update_(size_t index, float delta) {
        while (true) {
            weight_sum_[index] += delta;
            if (index == 0) {
                break;
            }
            index = (index - 1) / 2;
        }
    }

    vector<float> weight_sum_;
    size_t base_;
    size_t width_;
    vector<IRegion*> regions_;
    unordered_map<IRegion*, size_t> indexes_;
};


class WeightedRegionSelectorFactory : public RegionSelectorFactory {
public:
    virtual IRegionSelector* CreateRegionSelector() {
        return new WeightedRegionSelector();
    }
};


class BasicDartsNet : public IDartsNet {
public:
    BasicDartsNet(float d1, float d2) :
            max_dist_(d2 * d2 * 4), d1_(d1 * d1), d2_(d2 * d2) {}

    Dart* WrapDart(const Dart &dart) const { return new Node_(dart); }

    void DestroyAll(Dart *master) {
        unordered_set<Node_*> removed_nodes;
        queue<Node_*> qu;
        qu.push(dynamic_cast<Node_*>(master));
        while (!qu.empty()) {
            Node_* node = qu.front();
            qu.pop();
            if (removed_nodes.count(node) == 0) {
                for (Node_* nei : node->neighbors) {
                    qu.push(nei);
                }
                delete node;
                removed_nodes.insert(node);
            }
        }
    }

    size_t GetNears(Dart *master, const Point2D &position,
                    vector<Dart*> *out) const {
        Node_* node = dynamic_cast<Node_*>(master);
        size_t ret = 1;
        out->resize(1);
        out->at(0) = master;
        for (Node_* nei : node->neighbors) {
            float l2 = (nei->position - position).LengthSquared();
            if (l2 > d2_) {
                continue;
            }
            out->emplace_back(nei);
            if (l2 < d1_) {
                swap(out->at(ret), out->back());
                ret += 1;
            }
        }
        return ret;
    }

    void AddNeighbor(Dart *master, Dart *neighbor) {
        Node_* m = dynamic_cast<Node_*>(master);
        Node_* n = dynamic_cast<Node_*>(neighbor);

        unordered_set<Node_*> nodes;
        for (Node_ *nei : m->neighbors) {
            nodes.insert(nei);
            for (Node_ *neii : nei->neighbors) {
                if (neii != m) {
                    nodes.insert(neii);
                }
            }
        }

        Connect_(m, n);
        for (Node_ *node : nodes) {
            TryConnect_(node, n);
        }
    }

private:
    struct Node_ : Dart {
        Node_(const Dart &d) : Dart(d) {}

        vector<Node_*> neighbors;
    };

    void TryConnect_(Node_* a, Node_* b) {
        if ((a->position - b->position).LengthSquared() < max_dist_) {
            Connect_(a, b);
        }
    }

    static void Connect_(Node_* a, Node_* b) {
        a->neighbors.push_back(b);
        b->neighbors.push_back(a);
    }

    float max_dist_;
    float d1_;
    float d2_;
};


class PreSampler {
public:
    PreSampler(float r_ratio,
               float size,
               RegionFactory *region_factory,
               RegionSelectorFactory *region_selector_factory) :
            r_ratio_(r_ratio),
            size_(size),
            region_factory_(region_factory),
            region_selector_(region_selector_factory->CreateRegionSelector()),
            darts_net_(new BasicDartsNet(r_ratio_ + 1.f, r_ratio_ * 2)),
            EclipseAndAddNewRegions_(this),
            EclipseOldRegions_(this) {
        Assert(size > 1.f);
        Assert(1.f <= r_ratio_ && r_ratio_ <= 2.f);

        Point2D p[] = {
            Point2D(0, 0),
            Point2D(size, 0),
            Point2D(size, size),
            Point2D(0, size)
        };
        for (int i = 0; i < 4; ++i) {
            bounds_[i] = Line2D(p[i], p[(i + 1) % 4] - p[i]);
        }
    }

    vector<Point2D> Sample() {
        Dart *old_dart;
        Dart *new_dart;
        vector<Dart*> neighbors;

        vector<Point2D> ret;

        Dart* root = SampleFirstPoint_();
        ret.emplace_back(root->position);

        while (!regions_.empty()) {
            ThrowDart_(&old_dart, &new_dart);
            ret.emplace_back(new_dart->position);

            size_t nct = darts_net_->GetNears(
                    old_dart, new_dart->position, &neighbors);

            EclipseAndAddNewRegions_(
                    region_factory_->CreateRegions(new_dart), neighbors);

            neighbors.resize(nct);
            EclipseOldRegions_(new_dart->position, neighbors);

            darts_net_->AddNeighbor(old_dart, new_dart);
        }

        darts_net_->DestroyAll(root);

        return ret;
    }

private:
    Dart* SampleFirstPoint_() {
        Dart* d = darts_net_->WrapDart(
                Dart(Point2D(rng.RandomFloat(), rng.RandomFloat())));
        for (IRegion* reg : region_factory_->CreateRegions(d)) {
            AddRegion_(reg);
        }
        return d;
    }

    void ThrowDart_(Dart **old_dart, Dart **new_dart) {
        IRegion* reg = region_selector_->Select(rng);
        *old_dart = reg->dart();
        *new_dart = darts_net_->WrapDart(Dart(reg->SelectPoint(rng)));
    }

    class EclipseAndAddNewRegionsHandler_ {
    public:
        EclipseAndAddNewRegionsHandler_(PreSampler *pre_sampler) :
                pre_sampler_(pre_sampler) {}

        void operator()(const vector<IRegion*> &regs,
                        const vector<Dart*> &neighbors) {
            size_t ct = 0;
            for (IRegion *reg : regs) {
                buf1_.resize(1);
                buf1_[0] = reg;
                for (Dart *neighbor : neighbors) {
                    EclipseAndSwap_(
                            Circle(neighbor->position, pre_sampler_->r_ratio_));
                }
                for (size_t i = 0; i < 4; ++i) {
                    EclipseAndSwap_(pre_sampler_->bounds_[i]);
                }
                for (IRegion *new_reg : buf1_) {
                    pre_sampler_->AddRegion_(new_reg);
                }
                ct += buf1_.size();
            }
        }

    private:
        template <typename T>
        void EclipseAndSwap_(const T &obj) {
            buf2_.clear();
            for (IRegion *r : buf1_) {
                if (r->Eclipse(obj, &buf_e_)) {
                    for (IRegion *new_reg : buf_e_) {
                        buf2_.emplace_back(new_reg);
                    }
                    delete r;
                } else {
                    buf2_.emplace_back(r);
                }
            }
            swap(buf1_, buf2_);
        }

        PreSampler* pre_sampler_;

        vector<IRegion*> buf_e_;
        vector<IRegion*> buf1_;
        vector<IRegion*> buf2_;
    };

    class EclipseOldRegionsHandler_ {
    public:
        EclipseOldRegionsHandler_(PreSampler* pre_sampler) :
                pre_sampler_(pre_sampler) {}

        void operator()(const Point2D &p, const vector<Dart*> &neighbors) {
            unordered_set<IRegion*> regions_to_be_removed;
            for (Dart *neighbor : neighbors) {
                for (IRegion *old_reg : neighbor->regions) {
                    if (old_reg->Eclipse(Circle(p, 1.f), &buf_)) {
                        for (IRegion *new_reg : buf_) {
                            pre_sampler_->AddRegion_(new_reg);
                        }
                        regions_to_be_removed.insert(old_reg);
                    }
                }
            }
            for (IRegion* reg : regions_to_be_removed) {
                pre_sampler_->RemoveRegion_(reg);
            }
        }

    private:
        PreSampler* pre_sampler_;

        vector<IRegion*> buf_;
    };

    void AddRegion_(IRegion* reg) {
        reg->dart()->regions.insert(reg);
        regions_.insert(reg);
        region_selector_->Add(reg);
    }

    void RemoveRegion_(IRegion* reg) {
        reg->dart()->regions.erase(reg);
        regions_.erase(reg);
        region_selector_->Remove(reg);
    }

    float r_ratio_;
    float size_;
    RegionFactory *region_factory_;
    IRegionSelector *region_selector_;
    IDartsNet *darts_net_;

    EclipseAndAddNewRegionsHandler_ EclipseAndAddNewRegions_;
    EclipseOldRegionsHandler_ EclipseOldRegions_;

    unordered_set<IRegion*> regions_;
    Line2D bounds_[4];

    RNG rng;
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
            time_samples_buf_(new float[pixel_samples]) {
        Assert(samplesPerPixel > 0);
    }

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


RegionSelectorFactory* RegionSelectorFactory::CreateRegionSelectorFactory(
        bool weighted) {
    if (!weighted) {
        return new UniformRegionSelectorFactory();
    }
    if (weighted) {
        return new WeightedRegionSelectorFactory();
    }
    Assert(false);
    return NULL;
}


Sampler *CreateScallopedSampler(
        const ParamSet &params, const Film *film, const Camera *camera) {
    float r_ratio = params.FindOneFloat("rratio", 1.f);
    float size = params.FindOneFloat("size", 1000.f);
    bool weighted = params.FindOneBool("weighted", false);

    RegionFactory *region_factory = RegionFactory::CreateRegionFactory(r_ratio);
    Assert(region_factory != NULL);

    RegionSelectorFactory *region_selector_factory =
            RegionSelectorFactory::CreateRegionSelectorFactory(weighted);
    Assert(region_selector_factory != NULL);

    vector<Point2D> points(
            PreSampler(r_ratio, size,
                       region_factory, region_selector_factory).Sample());

    delete region_factory;
    delete region_selector_factory;

    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);

    int pixel_samples = params.FindOneInt("pixelsamples", 1);

    std::shared_ptr<const ISquareSampler> square_sampler(
            new TreeStyleArraySquareSampler(size, points));

    return new ScallopedSampler(xstart, xend, ystart, yend,
                                pixel_samples,
                                camera->shutterOpen, camera->shutterClose,
                                square_sampler);
}
