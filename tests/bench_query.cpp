#include <cassert>
#include <cmath>
#include <vector>
#include <set>
#include <string>
#include <tuple>
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <stdint.h>
#include <getopt.h>
#include <sys/time.h>

#include <boost/algorithm/string/replace.hpp>

#include "string_fmt.hpp"
#include "../geotree.hpp"
#include "../geodensity_bounded.hpp"


using namespace std;
using namespace geotools;


struct Duration {
    int64_t sec;
    int64_t usec;

    Duration(int64_t sec, int64_t usec) : sec(sec), usec(usec) {}

    string str(const string &fmt = "%S.%f") const {
        // hacks
        string ans = boost::replace_all_copy(fmt, "%S", strfmt("%ld", this->sec));
        ans = boost::replace_all_copy(ans, "%f", strfmt("%03u", this->usec / 1000));
        return ans;
    }

    double to_seconds() const {
        return double(this->sec) + double(this->usec) / 1000000;
    }

    static Duration US(int64_t usec) {
        return Duration(usec / 1000000, usec % 1000000);
    }
};


struct Time {
    uint64_t seconds;
    uint64_t micro_seconds;

    Time() : seconds(0), micro_seconds(0) {}
    Time(uint64_t seconds, uint64_t micro_seconds = 0)
        : seconds(seconds), micro_seconds(micro_seconds)
    {}

    static Time Now() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return Time(tv.tv_sec, tv.tv_usec);
    }

    string str(const string &fmt = "%Y-%m-%d_%H:%M:%S.%f") const {
        struct tm stm;
        time_t ts = this->seconds;
        gmtime_r(&ts, &stm);

        // hack
        string new_fmt = boost::replace_all_copy(fmt, "%f", strfmt("%03u", this->micro_seconds / 1000));

        vector<char> buf(100, '\0');
        size_t count = 0;
        do {
            count = strftime(buf.data(), buf.size(), new_fmt.c_str(), &stm);
            if (count == 0) {
                buf.resize(buf.size() * 2);
            }
        } while (count == 0);
    
        return string(buf.data());
    }

    Duration operator-(const Time &rhs) const {
        int64_t us = (this->seconds - rhs.seconds) * 1000000 + this->micro_seconds - rhs.micro_seconds;
        return Duration::US(us);
    }
};


void log(const string &prefix, const string &fmt, va_list va) {
    cerr << Time::Now().str() << " " << prefix << vstrfmt(fmt, va) << endl;
}


void info(string fmt, ...) {
    va_list va;
    va_start(va, fmt);
    log("INFO: ", fmt, va);
    va_end(va);
}


void err(string fmt, ...) {
    va_list va;
    va_start(va, fmt);
    log("ERR:  ", fmt, va);
    va_end(va);
}


struct Entry {
    uint32_t uid;
    float lon;
    float lat;
};


vector<Entry> load_file(const string &filename) {
    vector<Entry> ans;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        Entry entry;
        int lon = INT_MAX;
        int lat = INT_MAX;

        int scanned = sscanf(line.c_str(), "%u%d%d", &entry.uid, &lon, &lat);
        assert(scanned == 3);
        entry.lon = lon / 1e6;
        entry.lat = lat / 1e6;
        ans.push_back(entry);
    }

    return ans;
}


struct Arguments {
    string file;
    int32_t split;
    int no_sort;
    set<string> tests;
};


void print_help() {
    // TODO: ...
}


Arguments get_args(int argc, char **argv) {
    Arguments args;
    args.split = 1000;

    while (true) {
        static const struct option long_options[] = {
            {"file",        required_argument, 0, 'f'},
            {"split",       required_argument, 0, 's'},
            {"no-sort",     no_argument,       &args.no_sort, true},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "f:s:", long_options, &option_index);

        switch (c) {
        case 0:
            // flags 
            break;
        case -1:
            while (optind < argc) {
                args.tests.insert(argv[optind++]);
            }
            return args;
        case 'f':
            args.file = optarg;
            break;
        case 's':
            args.split = atol(optarg);
            if (args.split <= 0) {
                err("illegal --split value");
                exit(3);
            }
            break;
        case '?':
            print_help();
            exit(1);
            break;
        default:
            assert(!"Unreachable");
            exit(2);
        }
    }
}


struct DurationLogger {
    Time start;
    string msg;
    uint64_t reqs;

    template <class ...Args>
    DurationLogger(const string &fmt, Args && ...args)
        : start(Time::Now()), msg(strfmt(fmt, std::forward<Args>(args)...)), reqs(0)
    {
         info("[begin] %s", this->msg.c_str());
    };

    ~DurationLogger() {
        Duration d = Time::Now() - this->start;
        if (this->reqs != 0) {
            int rps = this->reqs / d.to_seconds();
            info("[duration:%s][rps:%d] %s", d.str().c_str(), rps, this->msg.c_str());
        } else {
            info("[duration:%s] %s", d.str().c_str(), this->msg.c_str());
        }
    }

    void set_reqs(uint64_t reqs) {
        this->reqs = reqs;
    }
};


typedef GeoTree<uint32_t> Tree;


void bench_geo_density(
    const vector<Entry> &entries, Tree &tree, size_t nearbys, size_t set_count, size_t query_count)
{
    GeoDensityBounded den(0, set_count);
    {
        vector<tuple<int32_t, int32_t, uint32_t>> rs;
        for (size_t i = 0; i < set_count; ++i) {
            const Entry &e = entries[i * 2];
            if (!tree.is_valid(e.lon, e.lat)) {
                continue;
            }

            uint32_t radius = tree.get_nearby_radius_by_count(e.lon, e.lat, nearbys);
            rs.push_back(make_tuple(e.lon, e.lat, radius));
        }

        // set
        {
            DurationLogger dl("GeoDensity.set_radius [nearbys:%zu][set_count:%zu]", nearbys, set_count);
            dl.set_reqs(rs.size());

            for (auto e : rs) {
                int32_t lon, lat;
                uint32_t radius;
                tie(lon, lat, radius) = e;
                den.set_radius(lon, lat, radius);
            }
        }
    }
    info("GeoDensity stats after set: [nearbys:%zu][set_count:%zu] %s", nearbys, set_count, den.stats.repr().c_str());

    vector<tuple<int32_t, int32_t, uint32_t>> lon_lat_radius;
    for (size_t i = 0; i < query_count; ++i) {
        const Entry &e = entries[i * 2 + 1];
        if (!tree.is_valid(e.lon, e.lat)) {
            continue;
        }

        uint32_t radius = tree.get_nearby_radius_by_count(e.lon, e.lat, nearbys);
        if (radius != 0) {
            lon_lat_radius.push_back(make_tuple(e.lon, e.lat, radius));
        }
    }


    // vector<pair<uint32_t, uint32_t>> samples;

    size_t dropped_var = 0;
    vector<double> vars;
    vars.reserve(lon_lat_radius.size());
    // guess
    {
        DurationLogger dl("GeoDensity.guess_radius [nearbys:%zu][query_count:%zu]", nearbys, lon_lat_radius.size());
        dl.set_reqs(lon_lat_radius.size());

        for (auto e : lon_lat_radius) {
            int32_t lon, lat;
            uint32_t radius;
            tie(lon, lat, radius) = e;

            uint32_t est_radius = den.guess_radius(lon, lat);
            double ratio = (double)est_radius / radius;
            double var = abs(log2(ratio));

            if (var > 3.0) {
                dropped_var++;
            } else {
                vars.push_back(var);
            }

            // if (samples.size() < 100) {
            //     samples.push_back(make_pair(est_radius, radius));
            // }
        }
    }

    double avg_var = 0;
    for (size_t i = 0; i < vars.size(); ++i) {
        avg_var = (i * avg_var + vars[i]) / (i + 1);
    }
    info("GeoDensity stats after query: [nearbys:%zu][set_count:%zu][query_count:%zu][avg_var:%f][dropped_var:%zu] %s",
        nearbys, set_count, lon_lat_radius.size(), avg_var, dropped_var, den.stats.repr().c_str());
    // for (auto p : samples) {
    //     info("GeoDensity guess sample: %u/%u", p.first, p.second);
    // }
}


int main(int argc, char **argv) {
    Arguments args = get_args(argc, argv);
    int32_t option = GEO_OPT_NONE;
    if (args.no_sort) {
        option |= GEO_NO_SORT;
    }

    info("start loading file '%s'", args.file.c_str());
    vector<Entry> entries = load_file(args.file);
    info("loaded %zu entries", entries.size());


    // insertion
    Tree tree(args.split);
    for (size_t i = 0; i < entries.size(); ++i) {
        Entry &e = entries[i];
        if (!tree.is_valid(e.lon, e.lat)) {
            err("bad entry: (%f, %f)", e.lon, e.lat);
            continue;
        }

        tree.insert(e.uid, e.lon, e.lat);
    }
    info("inserted %zu unique entries", tree.size());

    const size_t QUERY_RUN = 10000;
    vector<size_t> nearby_counts = {1, 10, 50, 100, 200, 500, 1000};

    // query
    if (args.tests.count("tree")) {
        for (size_t nearbys : nearby_counts) {
            DurationLogger dl(
                "running %zu queries for nearby %zu. [split:%u][opt:%u]",
                QUERY_RUN, nearbys, args.split, option);
            dl.set_reqs(QUERY_RUN);

            for (size_t i = 0; i < QUERY_RUN; ++i) {
                Entry &e = entries[i];
                vector<Tree::Item> nearby = tree.get_nearby(e.lon, e.lat, nearbys, option);
                assert(nearby.size() == nearbys);
            }
        }
    }

    // geodensity
    if (args.tests.count("density")) {
        for (size_t nearbys : nearby_counts) {
            bench_geo_density(entries, tree, nearbys, 10000, 10000);
        }
    }

    return 0;
}
