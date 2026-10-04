#pragma once

#include <condition_variable>
#include <deque>
#include <fstream>
#include <mutex>
#include <string>
#include <vector>

enum class OutputMode {
    Edit,
    Matrix,
    RegionalMatrix
};

class MatchWriterQueue {
public:
    void Push(std::vector<std::string> batch);
    bool Pop(std::vector<std::string>& batch);
    void Close();

private:
    std::mutex mutex_;
    std::condition_variable cv_;
    std::deque<std::vector<std::string>> queue_;
    bool closed_ = false;
};

class MatchWriter {
public:
    MatchWriter(const std::string& path, bool withAlignment, OutputMode mode);
    void WriteHeader();
    void WriteBatch(const std::vector<std::string>& batch);

private:
    std::ofstream out_;
    bool withAlignment_ = false;
    OutputMode mode_ = OutputMode::Edit;
};
