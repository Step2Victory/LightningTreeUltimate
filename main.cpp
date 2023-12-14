#include "src/LightningTree.h"
#include <yaml-cpp/yaml.h>
#include <glog/logging.h>
#include <iostream>
#include <chrono>

int ReadResponse()
{
    int response;
    std::cin >> response;
    return response;
}

int main(int argc, char* argv[]) {

    // Initialize Google’s logging library.
    
    if (argc) {
        google::InitGoogleLogging(argv[0]);
    }
    google::SetLogDestination(google::GLOG_INFO,"./logs/INFO_");

    auto project_path = std::filesystem::current_path().parent_path();
    //auto project_path = ""
    YAML::Node config = YAML::LoadFile(project_path / "configs" / "main.yaml");
    auto path_data = project_path / config["path_to_data_for_python"].as<std::string>();
    auto start = std::chrono::system_clock::now();
    // std::cout<<"Запуск создания обекта LT\n";
    auto lt = LightningTree(project_path / "configs" / config["lt_config"].as<std::string>());
    // std::cout<<"Завершение создание объекта LT\n";
    auto end = std::chrono::system_clock::now();
    lt.AllParams();
    lt.Info();
    
    int n_iter = config["number_of_iterations"].as<int>();
    int input_period = config["input_files_iter_period"].as<int>();
    start = std::chrono::system_clock::now();
    lt.WriteResponse(1);
    auto response = ReadResponse();

    // std::cout<<"Запуск построения молнии\n";
    for (int i = 0; i < n_iter; ++i)
    {
        if(response == 0) return 0;
        try
        {
            lt.NextIter();

            if(i % input_period == 0){
                lt.ReturnFiles(path_data);
                lt.ReturnPhi(path_data);
                lt.Info();
                lt.WriteResponse(1);
                
                response = ReadResponse();
            }
        }
        catch(const std::exception& e)
        {
            LOG(INFO) << e.what() << '\n';
            lt.ReturnFiles(path_data);
            lt.ReturnPhi(path_data);
            lt.Info();
            lt.WriteResponse(0);
            return 0;
        }
        // lt->ReturnFiles(path_data / "vertex_table.txt", path_data / "edge_table.txt", path_data /"q_history_1.txt", path_data /"Q_history.txt");
    }
    end = std::chrono::system_clock::now();
    LOG(INFO) << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << '\n'; 
    // lt->Info();
    // // PrintCurrentState(*lt);
    lt.ReturnFiles(path_data);
    lt.ReturnPhi(path_data);
    lt.Info();
    lt.WriteResponse(0);
    google::ShutdownGoogleLogging();
    return 0;
}