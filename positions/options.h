#if !defined(__OPTIONS_H__)
#define __OPTIONS_H__

#include <string>

struct Options
{
    std::string input;
    bool normal_inpainting_enabled;
    bool position_inpainting_enabled;
    bool exclusion_enabled;
    double alpha;
    double lambda;
    double noise;
    bool preprocess_normals;
    double epsilon_start;
    double epsilon_finish;
    double epsilon_progression;
    double alpha_pos;
    double beta_pos;
    double alpha_he;
    double beta_he;
};

Options parse_options(int argc, char* argv[]);

#endif
