#include <core/stats.h>
#include <core/stringprint.h>

#include <chrono>
#include <string>
#include <array>

namespace lightfold {

    // Statistics Local Variables
    std::vector<std::function<void(StatsAccumulator&)>>* StatRegisterer::funcs;
    static StatsAccumulator statsAccumulator;

    // For a given profiler state (i.e., a set of "on" bits corresponding to
    // profiling categories that are active), ProfileSample stores a count of
    // the number of times that state has been active when the timer interrupt
    // to record a profiling sample has fired.
    struct ProfileSample {
        std::atomic<uint64_t> profilerState{ 0 };
        std::atomic<uint64_t> count{ 0 };
    };

    // We use a hash table to keep track of the profiler state counts. Because
    // we can't do dynamic memory allocation in a signal handler (and because
    // the counts are updated in a signal handler), we can't easily use
    // std::unordered_map.  We therefore allocate a fixed size hash table and
    // use linear probing if there's a conflict.
    static const int profileHashSize = 256;
    static std::array<ProfileSample, profileHashSize> profileSamples;

    static std::chrono::system_clock::time_point profileStartTime;

#ifdef PBRT_HAVE_ITIMER
    static void ReportProfileSample(int, siginfo_t*, void*);
#endif  // PBRT_HAVE_ITIMER

    // Statistics Definitions
    void ReportThreadStats() {
        static std::mutex mutex;
        std::lock_guard<std::mutex> lock(mutex);
        StatRegisterer::CallCallbacks(statsAccumulator);
    }

    void StatRegisterer::CallCallbacks(StatsAccumulator& accum) {
        for (auto func : *funcs) func(accum);
    }

    void PrintStats(FILE* dest) { statsAccumulator.Print(dest); }

    void ClearStats() { statsAccumulator.Clear(); }

    static void getCategoryAndTitle(const std::string& str, std::string* category,
        std::string* title) {
        const char* s = str.c_str();
        const char* slash = strchr(s, '/');
        if (!slash)
            *title = str;
        else {
            *category = std::string(s, slash - s);
            *title = std::string(slash + 1);
        }
    }

    void StatsAccumulator::Print(FILE* dest) {
        fprintf(dest, "Statistics:\n");
        std::map<std::string, std::vector<std::string>> toPrint;

        for (auto& counter : counters) {
            if (counter.second == 0) continue;
            std::string category, title;
            getCategoryAndTitle(counter.first, &category, &title);
            toPrint[category].push_back(StringPrintf(
                "%-42s               %12" PRIu64, title.c_str(), counter.second));
        }
        for (auto& counter : memoryCounters) {
            if (counter.second == 0) continue;
            std::string category, title;
            getCategoryAndTitle(counter.first, &category, &title);
            double kb = (double)counter.second / 1024.;
            if (kb < 1024.)
                toPrint[category].push_back(StringPrintf(
                    "%-42s                  %9.2f kB", title.c_str(), kb));
            else {
                float mib = kb / 1024.;
                if (mib < 1024.)
                    toPrint[category].push_back(StringPrintf(
                        "%-42s                  %9.2f MiB", title.c_str(), mib));
                else {
                    float gib = mib / 1024.;
                    toPrint[category].push_back(StringPrintf(
                        "%-42s                  %9.2f GiB", title.c_str(), gib));
                }
            }
        }
        for (auto& distributionSum : intDistributionSums) {
            const std::string& name = distributionSum.first;
            if (intDistributionCounts[name] == 0) continue;
            std::string category, title;
            getCategoryAndTitle(name, &category, &title);
            double avg = (double)distributionSum.second /
                (double)intDistributionCounts[name];
            toPrint[category].push_back(
                StringPrintf("%-42s                      %.3f avg [range %" PRIu64
                    " - %" PRIu64 "]",
                    title.c_str(), avg, intDistributionMins[name],
                    intDistributionMaxs[name]));
        }
        for (auto& distributionSum : floatDistributionSums) {
            const std::string& name = distributionSum.first;
            if (floatDistributionCounts[name] == 0) continue;
            std::string category, title;
            getCategoryAndTitle(name, &category, &title);
            double avg = (double)distributionSum.second /
                (double)floatDistributionCounts[name];
            toPrint[category].push_back(
                StringPrintf("%-42s                      %.3f avg [range %f - %f]",
                    title.c_str(), avg, floatDistributionMins[name],
                    floatDistributionMaxs[name]));
        }
        for (auto& percentage : percentages) {
            if (percentage.second.second == 0) continue;
            int64_t num = percentage.second.first;
            int64_t denom = percentage.second.second;
            std::string category, title;
            getCategoryAndTitle(percentage.first, &category, &title);
            toPrint[category].push_back(
                StringPrintf("%-42s%12" PRIu64 " / %12" PRIu64 " (%.2f%%)",
                    title.c_str(), num, denom, (100.f * num) / denom));
        }
        for (auto& ratio : ratios) {
            if (ratio.second.second == 0) continue;
            int64_t num = ratio.second.first;
            int64_t denom = ratio.second.second;
            std::string category, title;
            getCategoryAndTitle(ratio.first, &category, &title);
            toPrint[category].push_back(StringPrintf(
                "%-42s%12" PRIu64 " / %12" PRIu64 " (%.2fx)", title.c_str(), num,
                denom, (double)num / (double)denom));
        }

        for (auto& categories : toPrint) {
            fprintf(dest, "  %s\n", categories.first.c_str());
            for (auto& item : categories.second)
                fprintf(dest, "    %s\n", item.c_str());
        }
    }

    void StatsAccumulator::Clear() {
        counters.clear();
        memoryCounters.clear();
        intDistributionSums.clear();
        intDistributionCounts.clear();
        intDistributionMins.clear();
        intDistributionMaxs.clear();
        floatDistributionSums.clear();
        floatDistributionCounts.clear();
        floatDistributionMins.clear();
        floatDistributionMaxs.clear();
        percentages.clear();
        ratios.clear();
    }

    thread_local uint64_t ProfilerState;
    static std::atomic<bool> profilerRunning{ false };

    void InitProfiler() {
        // Access the per-thread ProfilerState variable now, so that there's no
        // risk of its first access being in the signal handler (which in turn
        // would cause dynamic memory allocation, which is illegal in a signal
        // handler).
        ProfilerState = ProfToBits(Prof::SceneConstruction);
        ClearProfiler();
        profileStartTime = std::chrono::system_clock::now();
        // Set timer to periodically interrupt the system for profiling
        profilerRunning = true;
    }

    static std::atomic<int> profilerSuspendCount{ 0 };

    void SuspendProfiler() { ++profilerSuspendCount; }

    void ResumeProfiler() { --profilerSuspendCount; }

    void ProfilerWorkerThreadInit() { }

    void ClearProfiler() {
        for (ProfileSample& ps : profileSamples) {
            ps.profilerState = 0;
            ps.count = 0;
        }
    }

    void CleanupProfiler() {
        profilerRunning = false;
    }

    static std::string timeString(float pct, std::chrono::system_clock::time_point now) {
        pct /= 100.;  // remap passed value to to [0,1]
        int64_t ns =
            std::chrono::duration_cast<std::chrono::nanoseconds>(now - profileStartTime).count();
        // milliseconds for this category
        int64_t ms = int64_t(ns * pct / 1000000.);
        // Peel off hours, minutes, seconds, and remaining milliseconds.
        int h = ms / (3600 * 1000);
        ms -= h * 3600 * 1000;
        int m = ms / (60 * 1000);
        ms -= m * (60 * 1000);
        int s = ms / 1000;
        ms -= s * 1000;
        ms /= 10;  // only printing 2 digits of fractional seconds
        return StringPrintf("%4d:%02d:%02d.%02d", h, m, s, ms);
    }

    void ReportProfilerResults(FILE* dest) { }

}  // namespace lightfold