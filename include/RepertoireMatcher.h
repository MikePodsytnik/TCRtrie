#pragma once

#include "CliConfig.h"

class RepertoireMatcher {
public:
    explicit RepertoireMatcher(CliConfig config);
    int Run();

private:
    CliConfig config_;
};
