#include "globalPara.h"

std::string input_file = "namaqualand";//" boulder_smooth"; //"RockSetBr_rotated"; //RockSetBr  test_cube
int Resolution = 128;

std::vector<double> Weights = { 0.8, 1.0, 1.2 };  //mst 分类系数权重：边界-内部，内部-内部，边界-边界
std::vector<double> KT_weights = { 0.5, 0.1, 0.1, 0.3 };  //kernel translucency权重：角度、长度、内外分布、水平垂直
double Isolevel = 0;
double Gauss_level = 0.5;

int PoresNum = 60;
double surface_ratio = 0.6; //表面核占比

int Max_degree = 5;

double Amplitude_min = 1.0;
double Amplitude_max = 1.0;
double Sigma_min = 0.015;
double Sigma_max = 0.033;
double W4sig_max = 3;
double W4sig_min = 25;
double Axis_max_ratio = 1.7;

double BaseLayer = 8; // 15  for Ceramic clay printing    10 for cand2
//250: 0.02, 0.04
//150~200: 0.025, 0.045
// 100: 0.033, 0.05
//120pores 0.02 0.037 0.58 0.38 A
//200pores 0.015 0.033 0.5 0.38 A
// 100 2 4 0.68 0.4 high_rock
//经验公式：sigma =(v/(6.8371*n*w))^(1/3)
// w for max = 5.0     for min = 30

double SmoothT = 15;    //Control the smoothing area size;  smoothT +, the effect -, approaching a regular union.
double Tube_radius_factor = 0.8; // Control the ratio of the pipeline radius relative to the Gaussian kernel radius; the mid_radius_factor + , the channel +.
double Safe_distance_ratio = 0.7;

double Trans_thres = 0.999;  //kernel translucency threshold 
double Adj_dis_thres = 0.35;
bool debug_show = false;
bool standard_show = true;
bool figure_show = false;
bool compare_show = false;
bool scr_figure = false;

bool Iso_kernel = false;
bool Direct_dis = false;

bool Handle_overlap = true;
bool optimize_debug = true;
bool Enable_noise = false;
bool topo_optimize = false;
bool dynamic_change_para = true;

bool fixed_graph_dual_render = true;
bool auto_capture_translucency_graph = true;

double Param2_Amplitude_min = 1.0;
double Param2_Amplitude_max = 1.0;
double Param2_Sigma_min = 0.015;
double Param2_Sigma_max = 0.033;
double Param2_W4sig_max = 18;
double Param2_W4sig_min = 18;
double Param2_Axis_max_ratio = 1.0;
double Param2_Gauss_level = 0.5;
double Param2_SmoothT = 15.0;
double Param2_Tube_radius_factor = 0.7;

float show_degree_x = -97.0;
float show_degree_y = 0.0;
double show_degree_z = -65.0;

std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (std::string::npos == first) return str;
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}


// 辅助函数：解析 vector<double> (例如 "0.8, 1.0, 1.2")
std::vector<double> parseVector(std::string valStr) {
    std::vector<double> vec;
    std::stringstream ss(valStr);
    std::string item;
    while (std::getline(ss, item, ',')) {
        vec.push_back(std::stod(item));
    }
    return vec;
}

bool loadParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open parameter file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        // 跳过注释和空行
        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") continue;

        std::stringstream ss(line);
        std::string key, valStr;

        // 假设格式为 "key = value"
        if (std::getline(ss, key, '=') && std::getline(ss, valStr)) {
            key = trim(key);
            valStr = trim(valStr);

            // 移除行尾可能存在的注释
            size_t commentPos = valStr.find("//");
            if (commentPos != std::string::npos) valStr = valStr.substr(0, commentPos);
            valStr = trim(valStr);

            try {
                if (key == "input_file") input_file = valStr;
                else if (key == "Resolution") Resolution = std::stoi(valStr);
                else if (key == "Weights") Weights = parseVector(valStr);
                else if (key == "KT_weights") KT_weights = parseVector(valStr);
                else if (key == "Isolevel") Isolevel = std::stod(valStr);
                else if (key == "Gauss_level") Gauss_level = std::stod(valStr);
                else if (key == "PoresNum") PoresNum = std::stoi(valStr);
                else if (key == "surface_ratio") surface_ratio = std::stod(valStr);
                else if (key == "Max_degree") Max_degree = std::stoi(valStr);
                else if (key == "Amplitude_min") Amplitude_min = std::stod(valStr);
                else if (key == "Amplitude_max") Amplitude_max = std::stod(valStr);
                else if (key == "Sigma_min") Sigma_min = std::stod(valStr);
                else if (key == "Sigma_max") Sigma_max = std::stod(valStr);
                else if (key == "W4sig_max") W4sig_max = std::stod(valStr);
                else if (key == "W4sig_min") W4sig_min = std::stod(valStr);
                else if (key == "Axis_max_ratio") Axis_max_ratio = std::stod(valStr);
                else if (key == "SmoothT") SmoothT = std::stod(valStr);
                else if (key == "Tube_radius_factor") Tube_radius_factor = std::stod(valStr);
                else if (key == "Safe_distance_ratio") Safe_distance_ratio = std::stod(valStr);
                else if (key == "Trans_thres") Trans_thres = std::stod(valStr);
                else if (key == "Adj_dis_thres") Adj_dis_thres = std::stod(valStr);
                else if (key == "BaseLayer") BaseLayer = std::stod(valStr);
                else if (key == "debug_show") debug_show = (valStr == "true" || valStr == "1");
                else if (key == "standard_show") standard_show = (valStr == "true" || valStr == "1");
                else if (key == "figure_show") figure_show = (valStr == "true" || valStr == "1");
                else if (key == "compare_show") compare_show = (valStr == "true" || valStr == "1");
                else if (key == "Iso_kernel") Iso_kernel = (valStr == "true" || valStr == "1");
                else if (key == "Handle_overlap") Handle_overlap = (valStr == "true" || valStr == "1");
                else if (key == "Direct_dis") Direct_dis = (valStr == "true" || valStr == "1");
                else if (key == "optimize_debug") optimize_debug = (valStr == "true" || valStr == "1");
                else if (key == "Enable_noise") Enable_noise = (valStr == "true" || valStr == "1");
                else if (key == "topo_optimize") topo_optimize = (valStr == "true" || valStr == "1");
                else if (key == "dynamic_change_para") dynamic_change_para = (valStr == "true" || valStr == "1");
                else if (key == "fixed_graph_dual_render") fixed_graph_dual_render = (valStr == "true" || valStr == "1");
                else if (key == "auto_capture_translucency_graph") auto_capture_translucency_graph = (valStr == "true" || valStr == "1");
                else if (key == "Param2_Amplitude_min") Param2_Amplitude_min = std::stod(valStr);
                else if (key == "Param2_Amplitude_max") Param2_Amplitude_max = std::stod(valStr);
                else if (key == "Param2_Sigma_min") Param2_Sigma_min = std::stod(valStr);
                else if (key == "Param2_Sigma_max") Param2_Sigma_max = std::stod(valStr);
                else if (key == "Param2_W4sig_max") Param2_W4sig_max = std::stod(valStr);
                else if (key == "Param2_W4sig_min") Param2_W4sig_min = std::stod(valStr);
                else if (key == "Param2_Axis_max_ratio") Param2_Axis_max_ratio = std::stod(valStr);
                else if (key == "Param2_Gauss_level") Param2_Gauss_level = std::stod(valStr);
                else if (key == "Param2_SmoothT") Param2_SmoothT = std::stod(valStr);
                else if (key == "Param2_Tube_radius_factor") Param2_Tube_radius_factor = std::stod(valStr);
            }
            catch (...) {
                std::cerr << "Warning: Failed to parse value for key: " << key << std::endl;
            }
        }
    }
    file.close();
    return true;
}

void writeVector(std::ofstream& ofs, const std::string& name,
    const std::vector<double>& vec)
{
    ofs << name << " = ";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        ofs << vec[i];
        if (i + 1 < vec.size()) ofs << ", ";
    }
    ofs << "\n";
}

void saveParameters(const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    ofs << "input_file = " << input_file << "\n";
    ofs << "Resolution = " << Resolution << "\n";

    writeVector(ofs, "Weights", Weights);
    writeVector(ofs, "KT_weights", KT_weights);

    ofs << "Isolevel = " << Isolevel << "\n";
    ofs << "Gauss_level = " << Gauss_level << "\n";

    ofs << "PoresNum = " << PoresNum << "\n";
    ofs << "surface_ratio = " << surface_ratio << "\n";

    ofs << "Max_degree = " << Max_degree << "\n";

    ofs << "Amplitude_min = " << Amplitude_min << "\n";
    ofs << "Amplitude_max = " << Amplitude_max << "\n";
    ofs << "Sigma_min = " << Sigma_min << "\n";
    ofs << "Sigma_max = " << Sigma_max << "\n";
    ofs << "W4sig_max = " << W4sig_max << "\n";
    ofs << "W4sig_min = " << W4sig_min << "\n";
    ofs << "Axis_max_ratio = " << Axis_max_ratio << "\n";

    ofs << "SmoothT = " << SmoothT << "\n";
    ofs << "Tube_radius_factor = " << Tube_radius_factor << "\n";
    ofs << "Safe_distance_ratio = " << Safe_distance_ratio << "\n";

    ofs << "Trans_thres = " << Trans_thres << "\n";
    ofs << "Adj_dis_thres = " << Adj_dis_thres << "\n";

    ofs << "BaseLayer = " << BaseLayer << "\n";
    
    ofs << "debug_show = " << std::boolalpha << debug_show << "\n";
    ofs << "standard_show = " << standard_show << "\n";
    ofs << "figure_show = " << figure_show << "\n";
    ofs << "compare_show = " << compare_show << "\n";

    ofs << "Iso_kernel = " << Iso_kernel << "\n";
    ofs << "Direct_dis = " << Direct_dis << "\n";

    ofs << "Handle_overlap = " << Handle_overlap << "\n";
    ofs << "optimize_debug = " << optimize_debug << "\n";
    ofs << "Enable_noise = " << Enable_noise << "\n";
    ofs << "topo_optimize = " << topo_optimize << "\n";
    ofs << "dynamic_change_para = " << dynamic_change_para << "\n";
    ofs << "fixed_graph_dual_render = " << fixed_graph_dual_render << "\n";
    ofs << "auto_capture_translucency_graph = " << auto_capture_translucency_graph << "\n";

    ofs << "Param2_Amplitude_min = " << Param2_Amplitude_min << "\n";
    ofs << "Param2_Amplitude_max = " << Param2_Amplitude_max << "\n";
    ofs << "Param2_Sigma_min = " << Param2_Sigma_min << "\n";
    ofs << "Param2_Sigma_max = " << Param2_Sigma_max << "\n";
    ofs << "Param2_W4sig_max = " << Param2_W4sig_max << "\n";
    ofs << "Param2_W4sig_min = " << Param2_W4sig_min << "\n";
    ofs << "Param2_Axis_max_ratio = " << Param2_Axis_max_ratio << "\n";
    ofs << "Param2_Gauss_level = " << Param2_Gauss_level << "\n";
    ofs << "Param2_SmoothT = " << Param2_SmoothT << "\n";
    ofs << "Param2_Tube_radius_factor = " << Param2_Tube_radius_factor << "\n";

    ofs.close();
}

std::string generateFilename()
{
    std::ostringstream oss;

    oss/* << std::fixed << std::setprecision(3)*/
        << input_file
        << "_" << Resolution
        << "_" << PoresNum
        << "_" << surface_ratio
        << "_" << W4sig_max
        << "_" << W4sig_min
        << "_" << Axis_max_ratio
        << "_" << Safe_distance_ratio
        << "_" << Trans_thres
        << "_" << Adj_dis_thres
        << "_"
        << b2i(Handle_overlap)
        << b2i(optimize_debug)
        << b2i(Enable_noise)
        << b2i(topo_optimize)
        << b2i(dynamic_change_para)
        <<".txt";

    return oss.str();
}
