#include "CliConfig.h"
#include "RepertoireMatcher.h"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    try {
        CliConfig config = ParseCli(argc, argv);
        RepertoireMatcher matcher(std::move(config));
        return matcher.Run();
    } catch (const std::exception& e) {
        std::cerr << "[tcrtrie] " << e.what() << '\n';
        return 1;
    }
}
