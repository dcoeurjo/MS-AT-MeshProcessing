#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <ctime>

using namespace std;

#include "Viewer.h"
#include "Image.h"
#include "Application.h"

namespace DDG
{
    // declare static member variables
    Mesh Viewer::mesh;
    GLuint Viewer::surfaceDL = 0;
    int Viewer::windowSize[2] = { 512, 512 };
    Camera Viewer::camera;
    Shader Viewer::shader;
    bool Viewer::renderWireframe = false;
    bool Viewer::renderReg = false;
    bool Viewer::renderVectorField = false;
    double Viewer::step = 10;
    double Viewer::delta = 1;
    double Viewer::alpha = 0;
    double Viewer::alpha_pos = 1e-2;
    double Viewer::beta_pos = 1e-1;
    double Viewer::alpha_he = 1e-1;
    double Viewer::beta_he = 2e-2;
    double Viewer::epsilon_start = 1;
    double Viewer::epsilon_finish = .1;
    double Viewer::epsilon_progression = 2;
    bool Viewer::normal_inpainting_enabled = false;
    bool Viewer::position_inpainting_enabled = true;
    bool Viewer::exclusion_enabled = false;

    void Viewer :: init( void )
    {
        restoreViewerState();
        initGLUT();
        GLenum err = glewInit();
        assert( err == GLEW_OK );
        setGL();
        initGLSL();
        updateDisplayList();
        glutMainLoop();
    }

    void Viewer :: initGLUT( void )
    {
        int argc = 0;
        vector< vector<char> > argv(1);

        // initialize window
        glutInitWindowSize( windowSize[0], windowSize[1] );
        glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
        glutInit( &argc, (char**)&argv );
        glutCreateWindow( "DDG" );

        // specify callbacks
        glutDisplayFunc  ( Viewer::display  );
        glutIdleFunc     ( Viewer::idle     );
        glutKeyboardFunc ( Viewer::keyboard );
        glutSpecialFunc  ( Viewer::special  );
        glutMouseFunc    ( Viewer::mouse    );
        glutMotionFunc   ( Viewer::motion   );

        // initialize menus
        int viewMenu = glutCreateMenu( Viewer::view );
        glutSetMenu( viewMenu );
        glutAddMenuEntry( "[f] Wireframe",  menuWireframe     );
        glutAddMenuEntry( "[t] Toggle regularized", menuDisplayReg );
        glutAddMenuEntry( "[k] Zoom In",    menuZoomIn        );
        glutAddMenuEntry( "[l] Zoom Out",   menuZoomOut       );
        glutAddMenuEntry( "[-] Decrease 1/lambda", menuDecreaseStep  );
        glutAddMenuEntry( "[+] Increase 1/lambda", menuIncreaseStep  );
        glutAddMenuEntry( "[v] VectorField",  menuVectorField );

      /*  int processMenu = glutCreateMenu( Viewer::menu );
        glutSetMenu( processMenu );
        glutAddMenuEntry( "Version A", menuProcessVersionA );
        glutAddMenuEntry( "Version B", menuProcessVersionB );
        glutAddMenuEntry( "Version C", menuProcessVersionC );
        glutAddMenuEntry( "[HeSchaefer13] cot", menuProcessHeSchaeferCot );
        glutAddMenuEntry( "[HeSchaefer13] area", menuProcessHeSchaeferArea );
*/
      
      
        int mainMenu = glutCreateMenu( Viewer::menu );
        glutSetMenu( mainMenu );
     //   glutAddSubMenu( "Process", processMenu );
        glutAddMenuEntry( "[c] raw -> reg", menuCopyToReg    );
        glutAddMenuEntry( "[x] reg -> raw", menuCopyToRaw    );
        glutAddMenuEntry( "[a] Preprocess normals", menuPreprocessNormals);
        glutAddMenuEntry( "[n] Process normals", menuProcessNormals    );
        glutAddMenuEntry( "[p] Process positions", menuProcessPositions    );
        glutAddMenuEntry( "[b] Add Noise",  menuNoise     );
        glutAddMenuEntry( "[r] Reset Mesh",       menuResetMesh  );
        glutAddMenuEntry( "[w] Write Mesh",       menuWriteMesh  );
        glutAddMenuEntry( "[s] Screenshot",      menuScreenshot );
        glutAddMenuEntry( "[esc] Exit",           menuExit       );
        glutAddSubMenu( "View", viewMenu );
        glutAttachMenu( GLUT_RIGHT_BUTTON );
    }

    void Viewer :: initGLSL( void )
    {
        for (const std::string& prefix_dir : {"shaders", "../shaders", "."})
        {
            const std::string vertex_shader_filename = prefix_dir+"/vertex.glsl";
            const std::string fragment_shader_filename = prefix_dir+"/fragment.glsl";
            if (!std::ifstream(vertex_shader_filename).good()) continue;
            if (!std::ifstream(fragment_shader_filename).good()) continue;
            shader.loadVertex(vertex_shader_filename.c_str());
            shader.loadFragment(fragment_shader_filename.c_str());
            return;
        }
        cerr << "could not find shaders" << endl;
        std::exit(42);
    }

    void Viewer :: menu( int value )
    {
        switch( value )
        {
            case( menuNoise ):
                mNoise();
                break;
            case( menuPreprocessNormals ):
                mPreprocessNormals();
                break;
            case( menuProcessVersionA ):
                mProcessVersionA();
                break;
            case( menuProcessVersionB ):
                mProcessVersionB();
                break;
            case( menuProcessVersionC ):
                mProcessVersionC();
                break;
            case( menuProcessHeSchaeferArea ):
                mProcessHeSchaefer(true);
                break;
            case( menuProcessHeSchaeferCot ):
                mProcessHeSchaefer(false);
                break;
            case( menuDisplayReg ):
                mDisplayReg();
                break;
            case( menuCopyToReg ):
                mCopyToReg();
                break;
            case( menuCopyToRaw ):
                mCopyToRaw();
                break;
            case( menuProcessNormals ):
                mProcessNormals();
                break;
            case( menuProcessPositions ):
                mProcessPositions();
                break;
            case( menuResetMesh ):
                mResetMesh();
                break;
            case( menuWriteMesh ):
                mWriteMesh();
                break;
            case( menuScreenshot ):
                mScreenshot();
                break;
            case( menuExit ):
                mExit();
                break;
            case( menuIncreaseStep ):
                mIncreaseStep();
                break;
            case( menuDecreaseStep ):
                mDecreaseStep();
                break;
            default:
                break;
        }
    }

    void Viewer :: view( int value )
    {
        switch( value )
        {
            case( menuWireframe ):
                mWireframe();
                break;
            case( menuDisplayReg ):
                mDisplayReg();
                break;
            case( menuZoomIn ):
                mZoomIn();
                break;
            case( menuZoomOut ):
                mZoomOut();
                break;
            case( menuVectorField ):
                mVectorField();
                break;
            default:
                break;
        }
    }

    void Viewer :: keyboard( unsigned char c, int x, int y )
    {
        switch( c )
        {
            case 'f':
                mWireframe();
                break;
            case 'a':
                mPreprocessNormals();
                break;
            case 't':
                mDisplayReg();
                break;
            case 'w':
                mWriteMesh();
                break;
            case 'b':
                mNoise();
                break;
            case 'r':
                mResetMesh();
                break;
            case 's':
                mScreenshot();
                break;
            case ' ':
                Application::run_debug(mesh);
                break;
            case 'c':
                mCopyToReg();
                break;
            case 'x':
                mCopyToRaw();
                break;
            case 'n':
                mProcessNormals();
                break;
            case 'p':
                mProcessPositions();
                break;
            case 'k':
                mZoomIn();
                break;
            case 'l':
                mZoomOut();
                break;
            case 27:
                mExit();
                break;
            case '-':
                mDecreaseStep();
                break;
            case '+':
                mIncreaseStep();
                break;
            case 'v':
                mVectorField();
                break;
            case '1':
                mSaveCamera();
                break;
            case '2':
                mLoadCamera();
                break;
            default:
                cout << "unknown keybinding " << c << endl;
                break;
        }
    }

    void Viewer::mSaveCamera()
    {
        cout << "saving camera" << endl;

        camera.setView();

        const Quaternion drag = camera.pDrag;
        const Quaternion click = camera.pClick;
        const Quaternion last = camera.rLast;

        std::ofstream handle(".at_camera.conf");
        const auto save_quaternion = [&handle](const Quaternion& quaternion)
        {
            handle << quaternion.s << endl;
            handle << quaternion.v.x << endl;
            handle << quaternion.v.y << endl;
            handle << quaternion.v.z << endl;
        };

        //save_quaternion(drag);
        //save_quaternion(click);
        save_quaternion(last);
    }

    void Viewer::mLoadCamera()
    {
        cout << "loading camera" << endl;

        std::ifstream handle(".at_camera.conf");
        const auto load_quaternion = [&handle]()
        {
            double s,vx,vy,vz;
            handle >> s;
            handle >> vx;
            handle >> vy;
            handle >> vz;
            assert( handle.good() );
            return Quaternion(s, Vector(vx,vy,vz));
        };

        camera.pDrag = 1;
        camera.pClick = 1;
        camera.momentum = 1;
        camera.rLast = load_quaternion();
    }

    void Viewer :: special( int i, int x, int y )
    {
        switch( i )
        {
            case GLUT_KEY_UP:
                mZoomIn();
                break;
            case GLUT_KEY_DOWN:
                mZoomOut();
                break;
            case GLUT_KEY_LEFT:
            case GLUT_KEY_RIGHT:
                camera.zoomStop();
                break;
            case 27:
                mExit();
                break;
            default:
                break;
        }
    }

    void Viewer :: mouse( int button, int state, int x, int y )
    {
        if( ( glutGetModifiers() & GLUT_ACTIVE_SHIFT) and state == GLUT_UP )
            pickVertex(x, y);
        else
            camera.mouse( button, state, x, y );
    }

    void Viewer :: motion( int x, int y )
    {
        camera.motion( x, y );
    }

    void Viewer :: idle( void )
    {
        camera.idle();
        glutPostRedisplay();
    }

    void Viewer :: storeViewerState( void )
    {
        ofstream out( ".viewer_state.txt" );

        out << camera.rLast[0] << endl;
        out << camera.rLast[1] << endl;
        out << camera.rLast[2] << endl;
        out << camera.rLast[3] << endl;

        GLint view[4];
        glGetIntegerv( GL_VIEWPORT, view );
        out << view[2] << endl;
        out << view[3] << endl;
    }

    void Viewer :: restoreViewerState( void )
    {
        ifstream in( ".viewer_state.txt" );
        if( !in.is_open() ) return;

        in >> camera.rLast[0];
        in >> camera.rLast[1];
        in >> camera.rLast[2];
        in >> camera.rLast[3];
        in >> windowSize[0];
        in >> windowSize[1];
    }

    void Viewer :: mIncreaseStep( void )
    {
        step += delta;
        std::cout << "1/lambda = " << step << std::endl;
    }

    void Viewer :: mDecreaseStep( void )
    {
        step -= delta;
        if (step <= 0) step = 1;
        std::cout << "1/lambda = " << step << std::endl;
    }

    void Viewer :: mCopyToReg( void )
    {
        for (Face& face : mesh.faces)
            face.normal_regularized = face.normal;

        for (Vertex& vertex : mesh.vertices)
        {
            vertex.position_regularized = vertex.position;
        }

        updateDisplayList();
    }

    void Viewer :: mCopyToRaw( void )
    {
        for (Face& face : mesh.faces)
            face.normal = face.normal_regularized;

        for (Vertex& vertex : mesh.vertices)
        {
            vertex.position = vertex.position_regularized;
        }

        updateDisplayList();
    }

    void Viewer :: mProcessVersionA( void )
    {
        Application::run_version_a(alpha, 1/step, epsilon_start, epsilon_finish, epsilon_progression, alpha_pos, beta_pos, normal_inpainting_enabled, position_inpainting_enabled, exclusion_enabled, mesh);
        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mProcessVersionB( void )
    {
        Application::run_version_b(alpha, 1/step, epsilon_start, epsilon_finish, epsilon_progression, alpha_pos, beta_pos, normal_inpainting_enabled, position_inpainting_enabled, exclusion_enabled, mesh);
        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mProcessVersionC( void )
    {
        Application::run_version_c(alpha, 1/step, epsilon_start, epsilon_finish, epsilon_progression, alpha_pos, beta_pos, normal_inpainting_enabled, position_inpainting_enabled, exclusion_enabled, mesh);
        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mProcessHeSchaefer( const bool use_area )
    {
        Application::run_heschaefer(alpha_he, beta_he, 1e-3, 1e3, sqrt(2), use_area, mesh);
        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mPreprocessNormals()
    {
        Application::store_regularized_normals(Application::preprocess_normals(Application::compute_normals(mesh), mesh), mesh);

        for (Vertex& vertex : mesh.vertices)
        {
            vertex.position_regularized = vertex.position;
        }

        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mProcessNormals( void )
    {
        std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        Application::run_normals(alpha, 1/step, epsilon_start, epsilon_finish, epsilon_progression, normal_inpainting_enabled, exclusion_enabled, mesh);
        end = std::chrono::system_clock::now();
        int elapsed_mseconds = std::chrono::duration_cast<std::chrono::milliseconds>
                              (end-start).count();

        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

        std::cout << "finished computation at " << std::ctime(&end_time)
             << "elapsed time: " << elapsed_mseconds << "ms\n";
        std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

        for (Vertex& vertex : mesh.vertices)
        {
            vertex.position_regularized = vertex.position;
        }

        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mProcessPositions( void )
    {
      std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();
      Application::run_positions(alpha_pos, beta_pos, position_inpainting_enabled, mesh);
      end = std::chrono::system_clock::now();
      int elapsed_mseconds = std::chrono::duration_cast<std::chrono::milliseconds>
                          (end-start).count();

      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
      std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

      std::cout << "finished computation at " << std::ctime(&end_time)
           << "elapsed time: " << elapsed_mseconds << "ms\n";
      std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
        renderReg = true;
        updateDisplayList();
    }

    void Viewer :: mNoise()
    {
        mesh.noise(.025);
        updateDisplayList();
    }

    void Viewer :: mResetMesh( void )
    {
        mesh.reload();
        updateDisplayList();
    }

    void Viewer :: mWriteMesh( void )
    {
        mCopyToRaw();
        mesh.writeMesh( "out_mesh.obj" );
        mesh.writeEdgeTube( "out_tube.obj", 8e-3, 8 );
        mesh.writeEdgeStrips("out_strips.obj");
        mesh.writeEdgeCSV( "out_edge.csv" );
    }

    void Viewer :: mExit( void )
    {
        //storeViewerState();
        exit( 0 );
    }

    void Viewer :: mDisplayReg( void )
    {
        renderReg = !renderReg;
        cout << "regularization " << ( renderReg ? "on" : "off" ) << endl;
        updateDisplayList();
    }

    void Viewer :: mWireframe( void )
    {
        renderWireframe = !renderWireframe;
        updateDisplayList();
    }

    void Viewer :: mVectorField( void )
    {
        renderVectorField = !renderVectorField;
        updateDisplayList();
    }

    void Viewer :: mZoomIn( void )
    {
        camera.zoomIn();
    }

    void Viewer :: mZoomOut( void )
    {
        camera.zoomOut();
    }

    void Viewer :: mScreenshot( void )
    {
        static int index = 0;

        // get window width and height
        GLint view[4];
        glGetIntegerv( GL_VIEWPORT, view );
        int w = view[2];
        int h = view[3];

        // get pixels
        Image image( w, h );
        glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );

        stringstream filename;
        filename << "frames/viewer" << setw(8) << setfill( '0' ) << index << ".tga";
        image.write( filename.str().c_str() );

        index++;
    }

    void Viewer :: display( void )
    {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
      glClearColor( 0.0, 0.0, 0.0, 1.0 );
        GLint viewport[4];
        glGetIntegerv( GL_VIEWPORT, viewport );

        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        glOrtho(0, static_cast<double>(viewport[2]), static_cast<double>(viewport[3]), 0, -1, 1);

        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();

        glBegin(GL_TRIANGLE_STRIP);
        if (renderReg) glColor4f(0,1,0,1);
        else glColor4f(1,0,0,1);
        glVertex3d(10,10,0);
        glVertex3d(20,10,0);
        glVertex3d(10,20,0);
        glVertex3d(20,20,0);
        glEnd();

        glClear( GL_DEPTH_BUFFER_BIT );

        shader.enable();

        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        const double aspect = static_cast<double>(viewport[2])/static_cast<double>(viewport[3]);
        const double fovy = 50.;
        const double clipNear = .01;
        const double clipFar = 1000.;
        gluPerspective( fovy, aspect, clipNear, clipFar );

        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();

        Quaternion    eye = Vector( 0., 0., -2.5*camera.zoom );
        Quaternion center = Vector( 0., 0., 0. );
        Quaternion     up = Vector( 0., 1., 0. );
        gluLookAt(   eye[1],    eye[2],    eye[3],
                center[1], center[2], center[3],
                up[1],     up[2],     up[3] );


        Quaternion r = camera.currentRotation();
        eye = r.conj() * eye * r;
        GLint uniformEye = glGetUniformLocation( shader, "eye" );
        glUniform3f( uniformEye, eye[1], eye[2], eye[3] );

        Quaternion light = Vector( -1., 1., -2. );
        light = r.conj() * light * r;
        GLint uniformLight = glGetUniformLocation( shader, "light" );
        glUniform3f( uniformLight, light[1], light[2], light[3] );

        camera.setView();
        callDisplayList();
        shader.disable();
        glutSwapBuffers();
    }

    void Viewer :: updateDisplayList( void )
    {
        if( surfaceDL )
        {
            glDeleteLists( surfaceDL, 1 );
            surfaceDL = 0;
        }

        surfaceDL = glGenLists( 1 );
        glNewList( surfaceDL, GL_COMPILE );
        setMeshMaterial();
        drawScene();
        glEndList();
    }

    void Viewer :: setGL( void )
    {
        glClearColor( .5, .5, .5, 1. );
        setLighting();
    }

    void Viewer :: setLighting( void )
    {
        GLfloat position[4] = { 20., 30., 40., 0. };
        glLightfv( GL_LIGHT0, GL_POSITION, position );
        glEnable( GL_LIGHT0 );
        glEnable( GL_NORMALIZE );
    }

    void Viewer :: setMeshMaterial( void )
    {
        GLfloat  diffuse[4] = { .8, .5, .3, 1. };
        GLfloat specular[4] = { .3, .3, .3, 1. };
        GLfloat  ambient[4] = { .2, .2, .5, 1. };

        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse  );
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular );
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient  );
        glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 16.      );
    }

    void Viewer :: callDisplayList( void )
    {
        glPushAttrib( GL_ALL_ATTRIB_BITS );
        glEnable( GL_DEPTH_TEST );
        glEnable( GL_LIGHTING );
        glCallList( surfaceDL );
        glPopAttrib();
    }

    void Viewer :: drawScene( void )
    {
        glPushAttrib( GL_ALL_ATTRIB_BITS );

        glEnable( GL_POLYGON_OFFSET_FILL );
        glPolygonOffset( 1., 1. );
        drawPolygons();
        glDisable( GL_POLYGON_OFFSET_FILL );

        if( renderWireframe ) drawWireframe();

        if( renderVectorField ) drawVectorField();

        drawIsolatedVertices();

        drawSelectedVertices();

        glPopAttrib();
    }

    void Viewer :: drawPolygons( void )
    {
        for( FaceCIter f  = mesh.faces.begin();
                f != mesh.faces.end();
                f ++ )
        {
            if( f->isBoundary() ) continue;

            glBegin( GL_POLYGON );

            Vector N = renderReg ? f->normal_regularized : f->normal;
            glNormal3dv( &N[0] );

            HalfEdgeCIter he = f->he;
            do
            {
                /*
                if( not renderWireframe )
                {
                    Vector N = he->vertex->normal();
                    glNormal3dv( &N[0] );
                }
                */

                const Vector color = he->vertex->getColor();
                glColor4f( color.x, color.y, color.z, 1. );
                if (renderReg) glVertex3dv( &he->vertex->position_regularized[0] );
                else glVertex3dv( &he->vertex->position[0] );

                he = he->next;
            }
            while( he != f->he );
            glEnd();
        }
    }

    void Viewer :: drawWireframe( void )
    {
        shader.disable();
        glPushAttrib( GL_ALL_ATTRIB_BITS );

        glDisable( GL_LIGHTING );

        glLineWidth(10);
        glBegin( GL_LINES );
        glColor4f( 1, 1, 0, 1 );
        for( EdgeCIter e  = mesh.edges.begin();
                e != mesh.edges.end();
                e ++ )
        {
            if (e->scalar > .5) continue;
            if (renderReg)
            {
                glVertex3dv( &e->he->vertex->position_regularized[0] );
                glVertex3dv( &e->he->flip->vertex->position_regularized[0] );
            } else {
                glVertex3dv( &e->he->vertex->position[0] );
                glVertex3dv( &e->he->flip->vertex->position[0] );
            }
        }
        glEnd();

        glLineWidth(1);
        glBegin( GL_LINES );
        glColor4f( 0, 0, 0, 1 );
        for( EdgeCIter e  = mesh.edges.begin();
                e != mesh.edges.end();
                e ++ )
        {
            if (e->scalar <= .5) continue;
            if (renderReg)
            {
                glVertex3dv( &e->he->vertex->position_regularized[0] );
                glVertex3dv( &e->he->flip->vertex->position_regularized[0] );
            } else {
                glVertex3dv( &e->he->vertex->position[0] );
                glVertex3dv( &e->he->flip->vertex->position[0] );
            }
        }
        glEnd();

        glPopAttrib();
    }

    void Viewer :: drawIsolatedVertices( void )
    {
        glPushAttrib( GL_ALL_ATTRIB_BITS );

        glPointSize( 5 );
        glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
        glEnable( GL_POINT_SMOOTH );
        glColor3f( 1., 0., 0. );

        glBegin( GL_POINTS );
        for( VertexCIter v  = mesh.vertices.begin();
                v != mesh.vertices.end();
                v ++ )
        {
            if( v->isIsolated() )
            {
                glVertex3dv( &v->position[0] );
            }
        }
        glEnd();

        glPopAttrib();
    }

    void Viewer :: drawVertices( void )
    {
        for( VertexCIter v = mesh.vertices.begin();
                v != mesh.vertices.end();
                v ++ )
        {
            glLoadName(v->index);
            glBegin(GL_POINTS);
            glVertex3dv( &v->position[0] );
            glEnd();
        }
    }

    void Viewer :: drawSelectedVertices( void )
    {
        shader.disable();
        glPushAttrib( GL_ALL_ATTRIB_BITS );

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.5, 0.5, 1 );

        double h = 0.75*mesh.meanEdgeLength();
        for( VertexCIter v = mesh.vertices.begin();
                v != mesh.vertices.end();
                v ++ )
        {
            if( v->tag )
            {
                glPushMatrix();
                glTranslated(v->position.x, v->position.y, v->position.z);
                glutSolidSphere(h, 10, 10);
                glPopMatrix();
            }
        }
        glEnd();

        glPopAttrib();
    }

    void Viewer :: drawVectorField( void )
    {
        shader.disable();
        glPushAttrib( GL_ALL_ATTRIB_BITS );

        double h = 0.25*mesh.meanEdgeLength();

        glDisable( GL_LIGHTING );
        glLineWidth( 2.0 );

        for( FaceCIter f  = mesh.faces.begin();
                f != mesh.faces.end();
                f ++ )
        {
            if( f->isBoundary() ) continue;
            Vector a = f->barycenter();
            Vector b = a + h*f->normal_regularized;
            Vector c = a + h*f->normal;

            /*
            Vector v = b - a;
            Vector v90 = cross(n, v);
            Vector p0 = b;
            Vector p1 = p0 - 0.2 * v - 0.1 * v90;
            Vector p2 = p0 - 0.2 * v + 0.1 * v90;

            Vector w = c - a;
            Vector q0 = c;
            Vector q1 = q0 - .2 * w - .1 * v90;
            Vector q2 = q0 - .2 * w + .1 * v90;
            glColor3f( 0., 1., 0. );
            */

            glBegin( GL_LINES );
            glColor3f( 0., 1., 0. );
            glVertex3dv( &a[0] );
            glVertex3dv( &b[0] );
            glEnd();

            glBegin( GL_LINES );
            glColor3f( 0., 0., 1. );
            glVertex3dv( &a[0] );
            glVertex3dv( &c[0] );
            glEnd();

            /*
            glBegin(GL_TRIANGLES);
            glColor3f( 0., 1., 0. );
            glVertex3dv( &p0[0] );
            glVertex3dv( &p1[0] );
            glVertex3dv( &p2[0] );
            glColor3f( 0., 0., 1. );
            glVertex3dv( &q0[0] );
            glVertex3dv( &q1[0] );
            glVertex3dv( &q2[0] );
            glEnd();
            */
        }

        glPopAttrib();
    }

    void Viewer :: pickVertex(int x, int y)
    {
        int width  = glutGet(GLUT_WINDOW_WIDTH );
        int height = glutGet(GLUT_WINDOW_HEIGHT);
        if( x < 0 || x >= width || y < 0 || y >= height ) return;

        int bufSize = mesh.vertices.size();
        GLuint* buf = new GLuint[bufSize];
        glSelectBuffer(bufSize, buf);

        GLint viewport[4];
        GLdouble projection[16];
        glGetIntegerv( GL_VIEWPORT, viewport );
        glGetDoublev(GL_PROJECTION_MATRIX, projection);

        glRenderMode(GL_SELECT);
        glInitNames();
        glPushName(0);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluPickMatrix(x, viewport[3]-y, 10, 10, viewport);
        glMultMatrixd(projection);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        drawVertices();
        glPopMatrix();

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();

        glMatrixMode(GL_MODELVIEW);
        long hits = glRenderMode(GL_RENDER);

        int index = -1;
        double min_z = 1.0e100;
        for( long i = 0; i < hits; ++i )
        {
            double distance = buf[4*i + 1];
            if( distance < min_z )
            {
                index = buf[4*i + 3];
                min_z = distance;
            }
        }
        delete[] buf;

        if (index >= 0)
        {
            mesh.vertices[index].toggleTag();
            updateDisplayList();
        }
    }
}
