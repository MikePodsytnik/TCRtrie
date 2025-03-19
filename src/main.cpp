#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "Trie.h"

bool CorrectSearchingCheck(const std::string& clonotype, const std::string& sequence, int max_miss_matches) {
    int miss_matches = 0;
    for (int i = 0 ; i < clonotype.size(); ++i) {
        if (clonotype[i] != sequence[i]) ++miss_matches;
        if(miss_matches > max_miss_matches) return false;
    }

    return true;
}

std::set<string> StupidFound(std::vector<std::string>& patterns, const std::string& clonotype, int max_miss_matches) {
    std::set<std::string> result;

    for (auto sequence : patterns) {
        int miss_matches = 0;
        bool should_insert = true;
        for (int i = 0 ; i < clonotype.size(); ++i) {
            if (clonotype[i] != sequence[i]) ++miss_matches;
            if (miss_matches > max_miss_matches) {
                should_insert = false;
                break;
            }
        }
        if (should_insert) result.insert(sequence);
    }

    return result;
}

void load_patterns(const std::string& filename, std::vector<std::string>& patterns) {
    std::ifstream file(filename);
    std::string line;
    while(std::getline(file, line)) {
        patterns.push_back(line);
    }
}

void sample_clonotypes(const std::vector<std::string>& patterns, std::vector<std::string>& clonotypes, size_t sample_size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, patterns.size() - 1);

    clonotypes.clear();
    for (int i = 0; i < sample_size; ++i) {
        clonotypes.push_back(patterns[dis(gen)]);
    }
    std::cout << "Сэмплировано " << sample_size << " клонотипов.\n";
}

int main() {
//    std::vector<size_t> build_counts = {1000000, 900000, 800000, 700000, 600000, 500000, 400000, 300000, 200000, 100000, 75000, 50000, 25000, 10000};
//    std::ofstream results("average_build.csv", std::ios::app);
//    results << "Len" << ", " << "Count" << ", " << "AVG" << '\n';
//    results.close();
//    for (int j = 10; j <= 20; ++j) {
//        std::ofstream results("average_build.csv", std::ios::app);
//        std::vector<std::string> patterns;
//        cout << "Len: " << j << '\n';
//        load_patterns("cdr_len" + to_string(j) + ".txt", patterns);
//        for (int i = 0; i < 1; ++i) {
//            for (auto build_count : build_counts) {
//                std::vector<std::string> clonotypes;
//                sample_clonotypes(patterns, clonotypes, build_count);
//
//                Trie trie(clonotypes);
//                auto avg = trie.CalculateAverageChildrenPerNode();
//                cout << "N: " << build_count << " AVG: " << avg << '\n';
//                results << j << ", " << build_count << ", " << avg << '\n';
//            }
//        }
//        results.close();
//    }
//    std::vector<std::string> patterns;
//    load_patterns("cdr_len15.txt", patterns);
//    std::vector<std::string> clonotypes;
//    sample_clonotypes(patterns, clonotypes, 1000000);
//    Trie trie(clonotypes);
//    std::string clonotype = "CASSLRQETVYGYTF";
//    int max_replacements = 5;
//    set<string> results = trie.Search(clonotype, max_replacements);
//    vector<string> build;
//    for (auto r : results) build.push_back(r);
//    Trie result_trie(build);
//    auto start_time = std::chrono::high_resolution_clock::now();
//    auto result = trie.Search(clonotype, 5);
//    auto end_time = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//    std::cout << elapsed_time.count() << std::endl;
//
//    bool is_correct;
//    for (const std::string& sequence : results) {
//        is_correct = CorrectSearchingCheck(clonotype, sequence, max_replacements);
//    }
//
//    if (!is_correct) {
//        std::cerr << "The search was not performed correctly!" << std::endl;
//    }

//    std::vector<std::string> patterns;
//    load_patterns("cdr3.txt", patterns);
//    std::vector<std::string> clonotypes;
//    sample_clonotypes(patterns, clonotypes, 150);
//    NeighborSearcher searcher = NeighborSearcher(patterns);
//    std::ofstream results("different_length_multithreading.csv", std::ios::app);
//    results << "MissMatches," << "Time";
//    results.close();
//    for (int mm = 1; mm < 5; ++mm) {
//        for (auto clonotype : clonotypes) {
//            if (clonotype.size() >= 21 && clonotype.size() <= 7) continue;
//            cout << clonotype << std::endl;
//            std::ofstream results("different_length_multithreading.csv", std::ios::app);
//            auto start_time = std::chrono::high_resolution_clock::now();
//            searcher.Search(clonotype, mm);
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//            results << mm << ", " << elapsed_time.count() << '\n';
//            cout << ", " << mm << ", " << elapsed_time.count() << '\n';
//            results.close();
//        }
//    }
//    std::vector<size_t> search_counts = {10, 100, 1000, 5000, 10000, 25000, 50000, 100000};
//    for (int mm = 1; mm < 4; ++mm) {
//        for (int i = 0; i < 100000; ++i) {
//            std::ofstream results("no_hyper_threading10^5.csv", std::ios::app);
//            std::vector<std::string> clonotypes;
//            sample_clonotypes(patterns, clonotypes, 50);
//            auto start_time = std::chrono::high_resolution_clock::now();
//            trie.SearchForAll(clonotypes, mm);
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//            results << mm << ", " << elapsed_time.count() << '\n';
//            cout << mm << ", " << elapsed_time.count() << '\n';
//            results.close();
//        }
//    }
//
//    return 0;
    std::vector<std::string> patterns;
    load_patterns("cdr3.txt", patterns);
    Trie trie = Trie(patterns);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto result = trie.SearchForAll(patterns, 3);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
    cout << "All vdj: " << elapsed_time.count() << '\n';

    return 0;
}
//int main() {
//    std::ofstream results("performance_results_2.csv");
//    if (!results.is_open()) {
//        std::cerr << "Не удалось открыть файл для записи результатов.\n";
//        return 1;
//    }
//
//    results << "Length\tMaxReplacements\tTime(ms)\n";
//
//    for (int length = 10; length <= 20; ++length) {
//        std::string filename = "cdr_len" + std::to_string(length) + ".txt";
//        std::vector<std::string> patterns;
//        load_patterns(filename, patterns, 100000);
//
//        Trie trie(patterns);
//        std::vector<std::string> clonotypes;
//        sample_clonotypes(patterns, clonotypes, 500);
//
//        for (auto clonotype : clonotypes) {
//            auto start_time = std::chrono::high_resolution_clock::now();
//            auto result = trie.Search(clonotype, 3);
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//
//            results << length << ',' << elapsed_time.count() << "\n";
//            std::cout << "Length: " << length << ", Size: " << result.size()
//                      << ", Time: " << elapsed_time.count() << " ms\n";
//        }
//    }
//
//    results.close();
//    return 0;
//}

//int main() {
//    std::vector<size_t> build_counts = {1000000};
//    std::vector<int> search_counts = {1000};
//    const size_t build_count = 1000000;
//    std::string filename = "cdr_len15.txt";
//    std::vector<std::string> patterns;
//    load_patterns(filename, patterns, build_count);
//    Trie trie(patterns);
//    std::vector<std::string> clonotypes;
//    int missmatch_count = 1;
//    for (int search_count: search_counts) {
//        sample_clonotypes(patterns, clonotypes, search_count);
//
//        std::string result_filename = "OMP_search" + std::to_string(missmatch_count) + ".csv";
//        std::ofstream results(result_filename, std::ios::app);
//        auto start_time = std::chrono::high_resolution_clock::now();
//        for (auto clonotype: clonotypes) {
//            std::set<std::string> result_set = trie.Search(clonotype, missmatch_count);
//        }
//        auto end_time = std::chrono::high_resolution_clock::now();
//
//        std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//
//        results << search_count << ", " << elapsed_time.count() << "\n";
//        std::cout << "Count of searches: " << search_count << " Time: " << elapsed_time.count() << " ms\n";
//
//        results.close();
//        std::cout << "Результаты для количества замен " << missmatch_count << " записаны в файл: " << result_filename
//                  << "\n";
//    }

//    int missmatch_count = 1;
//    std::string result_filename = "OMP_search_1m.csv";
//    std::ofstream results(result_filename, std::ios::app);
//    if (!results.is_open()) {
//        std::cerr << "Не удалось открыть файл для записи результатов: " << result_filename << "\n";
//    }
//
//    results << "Длина, Количество замен, Количество поисков, Время работы (мс)\n";
//
//    for (int search_count : search_counts) {
//        std::vector<std::string> clonotypes;
//        sample_clonotypes(patterns, clonotypes, search_count);
//
//        for (int i = 0; i < search_count; ++i) {
//            std::string clonotype = clonotypes[i % clonotypes.size()];
//            auto start_time = std::chrono::high_resolution_clock::now();
//            std::set<std::string> results_set = trie.Search(clonotype, missmatch_count);
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//            results << missmatch_count << ", " << search_count << ", " << elapsed_time.count() << "\n";
//            std::cout << to_string(i) + ". Mis count: " << missmatch_count << ", Result count: " << results_set.size()
//                    << ", Time: " << elapsed_time.count() << " ms\n";
//        }
//    }
//
//    return 0;
//}
//    for (int length = 15; length <= 15; ++length) {
//        std::string filename = "cdr_len" + std::to_string(length) + ".txt";
//        std::vector<std::string> patterns;
//        load_patterns(filename, patterns);
//
//        if (patterns.empty()) {
//            std::cerr << "Пропущена длина " << length << " из-за отсутствия данных.\n";
//            continue;
//        }
//
//        Trie ac(patterns);
//        std::cout << "Инициализирован Trie для длины " << length << " с " << patterns.size() << " паттернами.\n";
//
//        std::vector<std::string> clonotypes;
//        int missmatch_count = 1;
//        sample_clonotypes(patterns, clonotypes, 10000);
//
//        std::string result_filename = "stupid_results_len" + std::to_string(length) + ".csv";
//        std::ofstream results(result_filename, std::ios::app);
//        if (!results.is_open()) {
//            std::cerr << "Не удалось открыть файл для записи результатов: " << result_filename << "\n";
//            continue;
//        }
//
//        results << "Длина, Количество замен, Количество поисков, Время работы (мс)\n";
//
//        for (int search_count : search_counts) {
//            auto start_time = std::chrono::high_resolution_clock::now();
//
//            for (int i = 0; i < search_count; ++i) {
//                std::string clonotype = clonotypes[i % clonotypes.size()];
//                std::set<std::string> results_set = ac.Search(clonotype, missmatch_count);
//            }
//
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//
//            results << length << ", " << missmatch_count << ", " << search_count << ", " << elapsed_time.count() << "\n";
//            std::cout << "Length: " << length << ", SearchCount: " << search_count
//                      << ", Time: " << elapsed_time.count() << " ms\n";
//        }
//
//        results.close();
//        std::cout << "Результаты для длины " << length << " записаны в файл: " << result_filename << "\n";
//    }
//
//    return 0;
//}

//int main() {
//    std::vector<int> search_counts = {100, 500, 1000, 2500, 5000, 7500, 10000};
//
//    for (int length = 15; length <= 15; ++length) {
//        std::string filename = "cdr_len" + std::to_string(length) + ".txt";
//        std::vector<std::string> patterns;
//        load_patterns(filename, patterns);
//
//        if (patterns.empty()) {
//            std::cerr << "Пропущена длина " << length << " из-за отсутствия данных.\n";
//            continue;
//        }
//
//        Trie ac(patterns);
//        std::cout << "Инициализирован Trie для длины " << length << " с " << patterns.size() << " паттернами.\n";
//
//        std::vector<std::string> clonotypes;
//        int missmatch_count = 1;
//        sample_clonotypes(patterns, clonotypes, 10000);
//
//        std::string result_filename = "stupid_results_len" + std::to_string(length) + ".csv";
//        std::ofstream results(result_filename, std::ios::app);
//        if (!results.is_open()) {
//            std::cerr << "Не удалось открыть файл для записи результатов: " << result_filename << "\n";
//            continue;
//        }
//
//        results << "Длина, Количество замен, Количество поисков, Время работы (мс)\n";
//
//        for (int search_count : search_counts) {
//            auto start_time = std::chrono::high_resolution_clock::now();
//
//            for (int i = 0; i < search_count; ++i) {
//                std::string clonotype = clonotypes[i % clonotypes.size()];
//                auto result_set = StupidFound(patterns, clonotype, missmatch_count);
////                std::set<std::string> results_set = ac.Search(clonotype, missmatch_count);
//            }
//
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//
//            results << length << ", " << missmatch_count << ", " << search_count << ", " << elapsed_time.count() << "\n";
//            std::cout << "Length: " << length << ", SearchCount: " << search_count
//                      << ", Time: " << elapsed_time.count() << " ms\n";
//        }
//
//        results.close();
//        std::cout << "Результаты для длины " << length << " записаны в файл: " << result_filename << "\n";
//    }
//
//    return 0;
//}

//void load_patterns(const std::string& filename, std::vector<std::string>& patterns) {
//    std::ifstream file(filename);
//    std::string line;
//    while (std::getline(file, line)) {
//        patterns.push_back(line);
//    }
//}
//
//int main() {
//    std::ofstream results("performance_results.txt");
//    if (!results.is_open()) {
//        std::cerr << "Не удалось открыть файл для записи результатов.\n";
//        return 1;
//    }
//
//    results << "Length\tMaxReplacements\tTime(ms)\n";
//
//    for (int length = 12; length <= 20; ++length) {
//        std::string filename = "cdr_len" + std::to_string(length) + ".txt";
//        std::vector<std::string> patterns;
//        load_patterns(filename, patterns);
//
//        Trie ac(patterns);
//        std::string clonotype = "CASSLRQETVYGYTF"; // Пример клото-типа
//
//        for (int max_replacements = 1; max_replacements <= 5; ++max_replacements) {
//            auto start_time = std::chrono::high_resolution_clock::now();
//
//            std::set<std::string> results_set = ac.Search(clonotype, max_replacements);
//
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
//
//            results << length << "\t" << max_replacements << "\t" << elapsed_time.count() << "\n";
//            std::cout << "Length: " << length << ", MaxReplacements: " << max_replacements
//                      << ", Time: " << elapsed_time.count() << " ms\n";
//        }
//    }
//
//    results.close();
//    return 0;
//}

//int main() {
//    std::ifstream inputFile("../../bin/cdr_len15.txt");
//
//    std::string line;
//    std::vector<std::string> patterns;
//    while (std::getline(inputFile, line)) {
//        patterns.push_back(line);
//    }
//
//    Trie ac(patterns);
//    std::string clonotype = "CASSLRQETVYGYTF";
//    int max_replacements = 5;
//    set<string> results = ac.Search(clonotype, max_replacements);
//
//    bool is_correct;
//    for (const std::string& sequence : results) {
//        is_correct = CorrectSearchingCheck(clonotype, sequence, max_replacements);
//    }
//
//    if (!is_correct) {
//        std::cerr << "The search was not performed correctly!" << std::endl;
//    }
//
//    return 0;
//}