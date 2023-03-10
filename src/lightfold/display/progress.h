#pragma once

#include <atomic>
#include <chrono>
#include <string>
#include <thread>

namespace lightfold {

    class ProgressReporter {
    public:
        // ProgressReporter Public Methods
        ProgressReporter(int64_t totalWork, const std::string& title);
        ~ProgressReporter();
        void Update(int64_t num = 1) {
            if (num == 0 ) return;
            workDone += num;
        }
        float ElapsedMS() const {
            std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
            int64_t elapsedMS =
                std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
            return (float)elapsedMS;
        }
        void Done();

    private:
        // ProgressReporter Private Methods
        void PrintBar();

        // ProgressReporter Private Data
        const int64_t totalWork;
        const std::string title;
        const std::chrono::system_clock::time_point startTime;
        std::atomic<int64_t> workDone;
        std::atomic<bool> exitThread;
        std::thread updateThread;
    };

} // namespace lightfold