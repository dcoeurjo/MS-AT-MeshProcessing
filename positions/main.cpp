#include <iostream>
using namespace std;

#include <glew.h>
#include "Viewer.h"
#include "Application.h"
#include "DenseMatrix.h"
#include "options.h"
using namespace DDG;

int main( int argc, char** argv )
{
    const Options options = parse_options(argc, argv);

    Viewer viewer;
    viewer.alpha_pos = options.alpha_pos;
    viewer.beta_pos = options.beta_pos;
    viewer.alpha_he = options.alpha_he;
    viewer.beta_he = options.beta_he;
    viewer.alpha = options.alpha;
    viewer.step = 1.0/options.lambda;
    viewer.epsilon_start = options.epsilon_start;
    viewer.epsilon_finish = options.epsilon_finish;
    viewer.epsilon_progression = options.epsilon_progression;
    viewer.normal_inpainting_enabled = options.normal_inpainting_enabled;
    viewer.position_inpainting_enabled = options.position_inpainting_enabled;
    viewer.exclusion_enabled = options.exclusion_enabled;
    viewer.mesh.read( options.input );
    viewer.mesh.noise( options.noise );
    if (options.preprocess_normals)
    {
        cout << "preprocess normals" << endl;

        Application::store_normals(Application::preprocess_normals(Application::compute_normals(viewer.mesh), viewer.mesh), viewer.mesh);
    }
    viewer.init();

    return 0;
}
