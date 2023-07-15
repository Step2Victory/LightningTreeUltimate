#include "src/LightningTree.h"
#include <glog/logging.h>
#include <iostream>

int main(int argc, char* argv[]) {
    // Initialize Google’s logging library.
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::GLOG_INFO,"./logs/INFO_");

    auto project_path = std::filesystem::current_path();
    auto lt = LightningTree(project_path / "configs" / "testing.yaml");
    lt.AllParams();
    lt.Info();
    google::ShutdownGoogleLogging();
    return 0;
}