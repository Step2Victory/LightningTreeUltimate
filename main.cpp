#include "src/LightningTree.h"
//#include <glog/logging.h>
#include <iostream>
#include <chrono>

int ReadResponse()
{
    int response;
    std::cin >> response;
    return response;
}

void WriteAnswer(int response, const LightningTree& tree)
{
    tree.Info();
}

int main(int argc, char* argv[]) {
    // Initialize Google’s logging library.
    auto t = std::chrono::system_clock::now();
    //google::InitGoogleLogging(argv[0]);
    //google::SetLogDestination(google::GLOG_INFO,"./logs/INFO_");

    auto project_path = std::filesystem::current_path();
    auto path_data = project_path / "LightningTree_data";
    auto start = std::chrono::system_clock::now();
    //auto lt = LightningTree(project_path / "configs" / "testing.yaml");
    auto lt = LightningTree(
            /*.h = */ 100,
        /* .delta_t =*/ 10e-5,
        /* .r =*/ 0.01,
        /* .R =*/ 50,
        /* .periphery_size =*/ 1,
        /*.q_plus_max =*/ 15e-10,
        /*.q_minus_max =*/ -3e-9,
        /*.Q_plus_s =*/ 4e-4,
        /*.Q_minus_s =*/ 8e-4,
        /*.resistance =*/ 1,
        /*.E_plus =*/ 150000,
        /*.E_minus =*/ 300000,
        /*.alpha =*/ 5e-8,
        /*.beta =*/ 10,
        /*.sigma =*/ 1e-5,
        /*.start_r =*/ {0, 0, 7000},
        /*.end_r =*/ {0, 0, 11000},
        /*degree_probability_growth = */ 2.5,
        /*.seed =*/ 42);
    auto end = std::chrono::system_clock::now();
    lt.AllParams();
    lt.Info();

    int n_iter = 100;
    start = std::chrono::system_clock::now();
    /*WriteAnswer(1, lt);
    auto response = ReadResponse();*/

    for (int i = 0; i < n_iter; ++i)
    {
        //if(response == 0) return 0;
        try
        {
            lt.NextIter();

            if(i % 10 == 0){
                //lt.ReturnFiles(path_data);
                //lt.ReturnPhi(path_data, lt.start_r, lt.end_r);
                lt.Info();
                /*WriteAnswer(1, lt);
                
                response = ReadResponse();*/
            }
        }
        catch(const std::exception& e)
        {
            // std::cerr << e.what() << '\n';
            //lt.ReturnFiles(path_data);
            //lt.ReturnPhi(path_data, lt.start_r, lt.end_r);
            lt.Info();
            //WriteAnswer(0, lt);
            return 0;
        }
        // lt->Info();
        // lt->ReturnFiles(path_data / "vertex_table.txt", path_data / "edge_table.txt", path_data /"q_history_1.txt", path_data /"Q_history.txt");
    }
    end = std::chrono::system_clock::now();
    // std::cout << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << '\n'; 
    // lt->Info();
    // // PrintCurrentState(*lt);
    //lt.ReturnFiles(path_data);
    //lt.ReturnPhi(path_data, lt.start_r, lt.end_r);
    lt.Info();
    //WriteAnswer(0, lt);
    //google::ShutdownGoogleLogging();
    return 0;
}