#include "options.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/errors.hpp>

#include <string>
#include <iostream>

Options parse_options(int argc, char* argv[])
{
    namespace po = boost::program_options;
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::string;

    Options options;
    options.input = "";

    po::options_description po_options("3d-at-positions [options] surface.obj alpha lambda");
    po_options.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<string>(&options.input)->default_value(""), "input surface.obj")
        ("exclusion", po::value<bool>(&options.exclusion_enabled)->default_value(false), "enable feature exclusion in normal regularization (blue channel)")
        ("normal-inpainting", po::value<bool>(&options.normal_inpainting_enabled)->default_value(false), "enable inpainting in normal regularization (green channel)")
        ("normal-preprocess", po::value<bool>(&options.preprocess_normals)->default_value(false), "prefilter normal vector field")
        ("normal-alpha", po::value<double>(&options.alpha)->default_value(1e-1), "normal dimensionless alpha")
        ("normal-lambda", po::value<double>(&options.lambda)->default_value(5e-2), "normal dimensionless lambda")
        ("normal-epsilon-start,s", po::value<double>(&options.epsilon_start)->default_value(1), "normal dimensionless start epsilon")
        ("normal-epsilon-finish,f", po::value<double>(&options.epsilon_finish)->default_value(1e-1), "normal dimensionless finish epsilon")
        ("normal-epsilon-progression,p", po::value<double>(&options.epsilon_progression)->default_value(3), "normal epsilon progression")
        ("position-inpainting", po::value<bool>(&options.position_inpainting_enabled)->default_value(false), "enable inpainting in position regularization (green channel)")
        ("position-alpha", po::value<double>(&options.alpha_pos)->default_value(1), "position dimensionless alpha")
        ("position-beta", po::value<double>(&options.beta_pos)->default_value(5e-2), "position dimensionless beta")
       // ("he-alpha", po::value<double>(&options.alpha_he)->default_value(1e-1), "heshaefer dimensionless alpha")
       // ("he-beta", po::value<double>(&options.beta_he)->default_value(2e-2), "heshaefer dimensionless beta")
        ("noise,k", po::value<double>(&options.noise)->default_value(0), "Noise level")
        ;
    po::positional_options_description positional;
    positional.add("input",1);
    positional.add("alpha",1);
    positional.add("lambda",1);

    try
    {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(po_options).positional(positional).run(), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            cout << po_options;
            std::exit(0);
        }

        if (options.input == "") throw po::required_option("input");

        cout << "#########################" << endl;
        cout << "input=" << options.input << endl;
        cout << "noise=" << options.noise << endl;
        cout << "## normals ##############" << endl;
        cout << "inpainting_enabled=" << options.normal_inpainting_enabled << endl;
        cout << "preprocess=" << options.preprocess_normals << endl;
        cout << "alpha=" << options.alpha << endl;
        cout << "lambda=" << options.lambda << endl;
        cout << "epsilon_start=" << options.epsilon_start << endl;
        cout << "epsilon_finish=" << options.epsilon_finish << endl;
        cout << "epsilon_progression=" << options.epsilon_progression << endl;
        cout << "## positions ##############" << endl;
        cout << "inpainting_enabled=" << options.position_inpainting_enabled << endl;
        cout << "alpha=" << options.alpha_pos << endl;
        cout << "beta=" << options.beta_pos << endl;
        cout << "#########################" << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        cout << po_options;
        std::exit(1);
    }

    return options;
}
