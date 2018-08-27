// -----------------------------------------------------------------------------
// libDDG -- Viewer.h
// -----------------------------------------------------------------------------
//
// Viewer provides a graphical user interface (GUI) for inspecting and
// interacting with a Mesh object.  Viewer methods are static in order
// to make them compatible with GLUT callbacks.
//

#ifndef DDG_VIEWER_H
#define DDG_VIEWER_H

#include <glew.h>
#include <glut.h>
#include "Mesh.h"
#include "Camera.h"
#include "Shader.h"
#include "options.h"

namespace DDG
{
    class Viewer
    {
        public:
            static void init( void );
            // displays the viewer until the program ends

            static Mesh mesh;
            // surface mesh visualized by Viewer

        protected:
            // init
            static void initGLUT( void );
            static void initGLSL( void );

            // GLUT callbacks
            static void display( void );
            static void idle( void );
            static void keyboard( unsigned char c, int x, int y );
            static void special( int i, int x, int y );
            static void mouse( int button, int state, int x, int y );
            static void motion( int x, int y );
            static void menu( int value );
            static void view( int value );

            // menu functions
            static void mNoise();
            static void mPreprocessNormals(void);
            static void mProcessVersionA(void);
            static void mProcessVersionB(void);
            static void mProcessVersionC(void);
            static void mProcessHeSchaefer(const bool use_area);
            static void mCopyToReg(void);
            static void mCopyToRaw(void);
            static void mProcessNormals( void );
            static void mProcessPositions( void );
            static void mResetMesh( void );
            static void mWriteMesh( void );
            static void mExit( void );
            static void mWireframe( void );
            static void mDisplayReg( void );
            static void mZoomIn( void );
            static void mZoomOut( void );
            static void mScreenshot( void );
            static void mIncreaseStep( void );
            static void mDecreaseStep( void );
            static void mVectorField( void );
            static void mSaveCamera();
            static void mLoadCamera();

            // unique identifiers for menus
            enum
            {
                menuNoise,
                menuPreprocessNormals,
                menuProcessVersionA,
                menuProcessVersionB,
                menuProcessVersionC,
                menuProcessHeSchaeferCot,
                menuProcessHeSchaeferArea,
                menuCopyToReg,
                menuCopyToRaw,
                menuProcessNormals,
                menuProcessPositions,
                menuResetMesh,
                menuWriteMesh,
                menuExit,
                menuWireframe,
                menuDisplayReg,
                menuZoomIn,
                menuZoomOut,
                menuScreenshot,
                menuIncreaseStep,
                menuDecreaseStep,
                menuVectorField
            };

            // draw routines
            static void setGL( void );
            static void setLighting( void );
            static void setMeshMaterial( void );
            static void callDisplayList( void );
            static void updateDisplayList( void );
            static void drawScene( void );
            static void drawPolygons( void );
            static void drawWireframe( void );
            static void drawVectorField( void );
            static void drawVertices( void );
            static void drawSelectedVertices( void );
            static void drawIsolatedVertices( void );
            static void pickVertex(int x, int y);

            static void storeViewerState( void );
            static void restoreViewerState( void );
            static int windowSize[2];

            static bool renderWireframe;
            // draw wireframe

            static bool renderVectorField;
            static bool renderReg;

            static Camera camera;
            // keeps track of view state

            static GLuint surfaceDL;
            // display list for mesh

            static Shader shader;
            // shader used to determine appearance of surface

        public:
            static double epsilon_start;
            static double epsilon_finish;
            static double epsilon_progression;
            static double alpha_he;
            static double beta_he;
            static double alpha_pos;
            static double beta_pos;
            static double alpha;
            static double step;
            static double noise;
            static double delta;
            static bool normal_inpainting_enabled;
            static bool position_inpainting_enabled;
            static bool exclusion_enabled;
            // fairing time step
    };
}

#endif

