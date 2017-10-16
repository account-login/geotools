#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <stdint.h>
#include <getopt.h>
#include <sys/time.h>

#include <boost/algorithm/string/replace.hpp>

#include "string_fmt.hpp"
#include "../geotree.hpp"


using namespace std;
using namespace geotree;


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
            // while (optind < argc) {
            //     args.args.push_back(argv[optind++]);
            // }
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


int main(int argc, char **argv) {
    Arguments args = get_args(argc, argv);
    int32_t option = GEO_OPT_NONE;
    if (args.no_sort) {
        option |= GEO_NO_SORT;
    }

    info("start loading file '%s'", args.file.c_str());
    vector<Entry> entries = load_file(args.file);
    info("loaded %zu entries", entries.size());

    typedef GeoTree<uint32_t> Tree;

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

    // query
    const size_t QUERY_RUN = 10000;
    vector<size_t> nearby_counts = {1, 10, 50, 100, 200, 500, 1000};
    for (auto nearbys = nearby_counts.begin(); nearbys != nearby_counts.end(); ++nearbys) {
        DurationLogger dl(
            "running %zu queries for nearby %zu. [split:%u][opt:%u]",
            QUERY_RUN, *nearbys, args.split, option);
        dl.set_reqs(QUERY_RUN);
        for (size_t i = 0; i < QUERY_RUN; ++i) {
            Entry &e = entries[i];
            vector<Tree::Item> nearby = tree.get_nearby(e.lon, e.lat, *nearbys, option);
            assert(nearby.size() == *nearbys);
        }
    }

    return 0;
}
