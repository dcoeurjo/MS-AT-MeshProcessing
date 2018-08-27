#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.h"
#include "Real.h"
#include "Utility.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"
#include <array>
#include <set>

namespace DDG
{
    class Application
    {
        public:
            typedef SparseMatrix<Real> Operator;
            typedef DenseMatrix<Real> ScalarField;
            typedef std::array<DenseMatrix<Real>, 3> VectorField;

            static void normalize_feature_scalars(ScalarField& feature_scalars)
            {
                const double feature_min = feature_scalars.minCoeff();
                const double feature_max = feature_scalars.maxCoeff();

                for (int index=0; index<feature_scalars.nRows(); index++)
                {
                    Real& feature = feature_scalars(index);
                    feature = (feature-feature_min)/(feature_max-feature_min);
                }
            }

            static void normalize_normals(VectorField& normals)
            {
                for (int index=0; index<normals[0].nRows(); index++)
                {
                    double norm = 0;
                    for (int dir=0; dir<3; dir++) norm += normals[dir](index)*normals[dir](index);
                    norm = sqrt(norm);
                    if (norm == 0) continue;
                    for (int dir=0; dir<3; dir++) normals[dir](index) /= norm;
                }
            }

            static void solve_normals(const Real& alpha, const Operator& sp0p, const Operator& sp1p, const Operator& ad1,
                                      const Operator& mm, const ScalarField& feature_scalars, const ScalarField& inpainting_scalars,
                                      const VectorField& normals_prescribed, VectorField& normals_estimated)
            {
                Operator uu_diag;
                Diag<Real>::build(mm*feature_scalars, uu_diag);
                Operator inpainting_diag;
                Diag<Real>::build(alpha*inpainting_scalars, inpainting_diag);
                Operator uu_half = uu_diag*ad1;
                Operator uu_operator = inpainting_diag*sp0p + uu_half.transpose()*sp1p*uu_half;
#pragma omp parallel for firstprivate(uu_operator)
                for (int dir=0; dir<3; ++dir)
                {
                    //righthand side alpha.g
                    ScalarField gg_alpha = inpainting_diag*sp0p*normals_prescribed[dir];

                    //solve
                    solveSymmetric(uu_operator, normals_estimated[dir], gg_alpha);
                    cout << "u" << dir << "=" << normals_estimated[dir].minCoeff() << "/" << normals_estimated[dir].maxCoeff() << endl;
                }
            }

            static void solve_feature_scalars(const Real& lambda, const Real& epsilon, const Operator& sp0, const Operator& sp1,
                                              const Operator& sp1p, const Operator& d0, const Operator& ad1, const Operator& mm,
                                              const ScalarField& one_scalars, const VectorField& normals_estimated,
                                              ScalarField& feature_scalars, const ScalarField& exclusion_scalars)
            {
                ScalarField vv_uu_norm_squared(ad1.nRows());
                vv_uu_norm_squared.zero();
                for (int dir=0; dir<3; ++dir)
                {
                    ScalarField foo = ad1*normals_estimated[dir];
                    for (int index=0; index<foo.length(); ++index)
                        foo(index) *= foo(index);
                    vv_uu_norm_squared += foo;
                }
                cout << "vv_uu_norm_squared=" << vv_uu_norm_squared.minCoeff() << "/" << vv_uu_norm_squared.maxCoeff() << endl;
                Operator vv_diag;
                Diag<Real>::build(vv_uu_norm_squared, vv_diag);
                Operator exclusion_diag;
                Diag<Real>::build(exclusion_scalars, exclusion_diag);
                Operator vv_second = mm.transpose()*sp1p*vv_diag*mm;
                Operator vv_half_first = d0;
                Operator vv_operator = Real(lambda/4.0/epsilon)*exclusion_diag*sp0 + Real(lambda*epsilon)*exclusion_diag*vv_half_first.transpose()*sp1*vv_half_first + vv_second;

                ScalarField lambda_over_4_epsilon = Real(lambda/4.0/epsilon)*exclusion_diag*sp0*one_scalars;
                solveSymmetric(vv_operator, feature_scalars, lambda_over_4_epsilon);
                cout << "v=" << feature_scalars.minCoeff() << "/" << feature_scalars.maxCoeff() << endl;
            }

            static ScalarField compute_positions(const Mesh& mesh)
            {
                const int nvertices = mesh.vertices.size();
                ScalarField positions_prescribed(3*nvertices);
                {
                    for (const Vertex& vertex : mesh.vertices)
                        for (int dir=0; dir<3; dir++)
                            positions_prescribed(vertex.index+nvertices*dir) = vertex.position[dir];
                }
                for (int dir=0; dir<3; dir++)
                {
                    cout << "q" << static_cast<char>('x'+dir) << "="
                        << *std::min_element(&positions_prescribed(nvertices*dir), &positions_prescribed(nvertices*(dir+1))) << "/"
                        << *std::max_element(&positions_prescribed(nvertices*dir), &positions_prescribed(nvertices*(dir+1))) << endl;
                }

                return positions_prescribed;
            }

            static void store_feature_scalars(const ScalarField& feature_scalars, Mesh& mesh)
            {
                const double feature_min = feature_scalars.minCoeff();
                const double feature_max = feature_scalars.maxCoeff();

                for (Vertex& vertex : mesh.vertices)
                    vertex.feature_scalar = feature_scalars(vertex.index);

                for (Edge& edge : mesh.edges)
                {
                    edge.scalar = feature_scalars(edge.he->vertex->index) > .5 || feature_scalars(edge.he->flip->vertex->index) > .5 ? 1 : 0;
                    //edge.scalar = (feature_scalars(edge.he->vertex->index)+feature_scalars(edge.he->flip->vertex->index))/2;
                }
            }

            static void store_normals(const VectorField& normals_estimated, Mesh& mesh)
            {
                assert( normals_estimated[0].nRows() == mesh.faces.size() );
                for (Face& face : mesh.faces)
                {
                    const Vector normal_estimated(normals_estimated[0](face.index), normals_estimated[1](face.index), normals_estimated[2](face.index));
                    face.normal = normal_estimated;
                }
            }

            static void store_regularized_normals(const VectorField& normals_estimated, Mesh& mesh)
            {
                assert( normals_estimated[0].nRows() == mesh.faces.size() );
                for (Face& face : mesh.faces)
                {
                    const Vector normal_estimated(normals_estimated[0](face.index), normals_estimated[1](face.index), normals_estimated[2](face.index));
                    face.normal_regularized = normal_estimated;
                }
            }

            static void store_regularized_positions(const ScalarField& positions_regularized, Mesh& mesh)
            {
                const int nvertices = mesh.vertices.size();

                {
                    ScalarField& foo = const_cast<ScalarField&>(positions_regularized);
                    for (int dir=0; dir<3; dir++)
                    {
                        cout << "p" << static_cast<char>('x'+dir) << "="
                            << *std::min_element(&foo(nvertices*dir), &foo(nvertices*(dir+1))) << "/"
                            << *std::max_element(&foo(nvertices*dir), &foo(nvertices*(dir+1))) << endl;
                    }
                }

                for (Vertex& vertex : mesh.vertices)
                    for (int dir=0; dir<3; dir++)
                        vertex.position_regularized[dir] = positions_regularized(vertex.index+nvertices*dir);
            }

            static VectorField preprocess_normals(const VectorField& original_face_normals, const Mesh& mesh)
            {
                VectorField vertex_normals;
                VectorField preprocessed_face_normals;

                for (int dir=0; dir<3; dir++)
                {
                    vertex_normals[dir] = DenseMatrix<Real>(mesh.vertices.size());
                    preprocessed_face_normals[dir] = DenseMatrix<Real>(mesh.faces.size());
                }

                for (const Face& face : mesh.faces)
                {
                    const Vector face_normal = face.normal;
                    for (int dir=0; dir<3; dir++)
                    {
                        HalfEdgeCIter he = face.he;
                        do
                        {
                            vertex_normals[dir](he->vertex->index) += face_normal[dir];
                            he = he->next;
                        }
                        while(he != face.he);
                    }
                }

                normalize_normals(vertex_normals);

                for (const Face& face : mesh.faces)
                {
                    Vector face_normal(0,0,0);

                    for (int dir=0; dir<3; dir++)
                    {
                        HalfEdgeCIter he = face.he;
                        do
                        {
                            face_normal[dir] += vertex_normals[dir](he->vertex->index);
                            he = he->next;
                        }
                        while(he != face.he);
                    }

                    for (int dir=0; dir<3; dir++)
                        preprocessed_face_normals[dir](face.index) = face_normal[dir];
                }

                normalize_normals(preprocessed_face_normals);

                return preprocessed_face_normals;
            }

            static VectorField compute_normals(const Mesh& mesh)
            {
                VectorField normals_prescribed;

                for (int dir=0; dir<3; dir++)
                    normals_prescribed[dir] = DenseMatrix<Real>(mesh.faces.size());

                for (const Face& face : mesh.faces)
                {
                    const Vector normal = face.normal;
                    for (int dir=0; dir<3; dir++)
                        normals_prescribed[dir](face.index) = normal[dir];
                }

                for (int dir=0; dir<3; dir++)
                    cout << "g" << dir << "=" << normals_prescribed[dir].minCoeff() << "/" << normals_prescribed[dir].maxCoeff() << endl;

                return normals_prescribed;
            }

            static void run_normals(const double alpha_normal, const double lambda_normal,
                                    const double epsilon_normal_start, const double epsilon_normal_finish,
                                    const double epsilon_normal_progression,
                                    const bool inpainting_enabled, const bool exclusion_enabled, Mesh& mesh)
            { // at normals
                cout << "////////////////" << endl;
                cout << "inpainting_enabled=" << inpainting_enabled << endl;
                cout << "exclusion_enabled=" << exclusion_enabled << endl;
                cout << "alpha_normal=" << alpha_normal << endl;
                cout << "lambda_normal=" << lambda_normal << endl;
                cout << "epsilon_normal_start=" << epsilon_normal_start << endl;
                cout << "epsilon_normal_finish=" << epsilon_normal_finish << endl;
                cout << "epsilon_normal_progression=" << epsilon_normal_progression << endl;
                cout << "////////////////" << endl;

                const size_t nfaces = mesh.faces.size();
                const size_t nedges = mesh.edges.size();
                const size_t nvertices = mesh.vertices.size();
                const double total_area = mesh.area();
                const double mean_edge_length = mesh.meanEdgeLength();
                //const double mean_diedral_angle = mesh.meanDiedralAngle();
                cout << "nfaces=" << nfaces << endl;
                cout << "nedges=" << nedges << endl;
                cout << "nvertices=" << nvertices << endl;
                cout << "total_area=" << total_area << endl;
                cout << "mean_edge_length=" << mean_edge_length << endl;
                //cout << "mean_diedral_angle=" << mean_diedral_angle << endl;

                cout << endl << "INITIALIZE DEC OPERATORS" << endl;

                Operator d0;
                ExteriorDerivative0Form<Real>::build( mesh, d0 );
                cout << "d0=" << d0.nRows() << "x" << d0.nColumns() << endl;

                Operator ad1;
                {
                    Operator d1;
                    ExteriorDerivative1Form<Real>::build( mesh, d1 );
                    ad1 = d1.transpose();
                }
                cout << "ad1=" << ad1.nRows() << "x" << ad1.nColumns() << endl;

                Operator mm;
                Wedge0FormTo1Form<Real>::build( mesh, mm );
                cout << "mm=" << mm.nRows() << "x" << mm.nColumns() << endl;

                Operator sp0;
                HodgeStar0Form<Real>::build( mesh, false, sp0 );
                for (int index=0; index<sp0.nRows(); index++) sp0(index, index) += total_area*1e-3/mesh.vertices.size();
                cout << "sp0=" << sp0.nRows() << "x" << sp0.nColumns() << endl;

                Operator sp1;
                HodgeStar1Form<Real>::build( mesh, false, sp1 );
                sp1 *= .5;
                cout << "sp1=" << sp1.nRows() << "x" << sp1.nColumns() << endl;

                Operator sp0p;
                HodgeStar2Form<Real>::build( mesh, true, sp0p );
                for (int index=0; index<sp0p.nRows(); index++) sp0p(index, index) += total_area*1e-3/mesh.faces.size();
                cout << "sp0p=" << sp0p.nRows() << "x" << sp0p.nColumns() << endl;

                Operator sp1p;
                HodgeStar1Form<Real>::build( mesh, true, sp1p );
                sp1p *= .5;
                cout << "sp1p=" << sp1p.nRows() << "x" << sp1p.nColumns() << endl;

                cout << endl << "SOLVING NORMAL REGULARIZATION" << endl;

                const VectorField normals_prescribed = compute_normals(mesh);
                VectorField normals_estimated = normals_prescribed;

                /*
                 ScalarField feature_exclusion_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    feature_exclusion_scalars(vertex.index) = exclusion_enabled ? (exclusion_enabled ? 1e1 : 1) : 1;
                 */
                // on an "inpainting" zone we reduce the lambda to increase the "perimeter" of the feature
                 ScalarField feature_exclusion_scalars(nvertices);
                 for (const Vertex& vertex : mesh.vertices)
                 {
                   bool painted = vertex.inpainting ;
                   feature_exclusion_scalars(vertex.index) = inpainting_enabled ? (painted ? 0.1 : 1) : 1;
                 }

                ScalarField inpainting_scalars(nfaces);
                for (const Face& face : mesh.faces)
                {
                    bool painted = face.he->vertex->inpainting || face.he->next->vertex->inpainting || face.he->next->next->vertex->inpainting;
                    inpainting_scalars(face.index) = inpainting_enabled ? ( painted ? 1e-8 : 1 ) : 1;
                }
                ScalarField one_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    one_scalars(vertex.index) = 1;
                ScalarField feature_scalars = one_scalars;
                cout << "feature_scalars=" << feature_scalars.minCoeff() << "/" << feature_scalars.maxCoeff() << endl;

                const Real alpha(alpha_normal/mean_edge_length/mean_edge_length);
                const Real lambda(lambda_normal/mean_edge_length);
                cout << "alpha=" << alpha << endl;
                cout << "lambda=" << lambda << endl;

                ScalarField feature_scalars_last = feature_scalars;
                for (double epsilon_normal=epsilon_normal_start; epsilon_normal>=epsilon_normal_finish; epsilon_normal/=epsilon_normal_progression)
                {
                    cout << "---------------------- epsilon_normal=" << epsilon_normal << endl;
                    const Real epsilon(epsilon_normal*mean_edge_length);
                    for (int iteration=0; iteration<10; iteration++)
                    {
                        cout << "*********** " << iteration << endl;

                        // solving normal field
                        solve_normals(alpha, sp0p, sp1p, ad1, mm, feature_scalars, inpainting_scalars, normals_prescribed, normals_estimated);

                        // solving feature field
                        solve_feature_scalars(lambda, epsilon, sp0, sp1, sp1p, d0, ad1, mm, one_scalars, normals_estimated, feature_scalars, feature_exclusion_scalars);

                        const double feature_scalars_variation_norm = (feature_scalars-feature_scalars_last).norm();
                        cout << "variation_norm=" << feature_scalars_variation_norm << endl;
                        const double feature_scalars_rel_variation_norm = feature_scalars_variation_norm/feature_scalars.norm();
                        cout << "rel_variation_norm=" << feature_scalars_rel_variation_norm << endl;
                        feature_scalars_last = feature_scalars;

                        if (feature_scalars_rel_variation_norm < 1e-2) break;

                    }
                }

                //normalize_normals(normals_estimated);
                store_feature_scalars(feature_scalars, mesh);
                //normalize_feature_scalars(feature_scalars);
                store_regularized_normals(normals_estimated, mesh);

                cout << "FINISHED NORMALS REGULARIZATION!!!!!!!!!" << endl << endl;
            }

            static void run_positions(const double alpha_raw, const double beta_raw, const bool inpainting_enabled, Mesh& mesh)
            {
                cout << "////////////////" << endl;
                cout << "inpainting_enabled=" << inpainting_enabled << endl;
                cout << "alpha_position=" << alpha_raw << endl;
                cout << "beta_position=" << beta_raw << endl;
                cout << "////////////////" << endl;

                const size_t nfaces = mesh.faces.size();
                const size_t nvertices = mesh.vertices.size();
                //const double mean_edge_length = mesh.meanEdgeLength();
                cout << "nfaces=" << nfaces << endl;
                cout << "nvertices=" << nvertices << endl;
                //cout << "mean_edge_length=" << mean_edge_length << endl;

                cout << endl << "SOLVING POSITION REGULARIZATION" << endl;

                //const Real alpha(alpha_raw/mean_edge_length/mean_edge_length);
                //const Real beta(beta_raw/mean_edge_length/mean_edge_length);
                const Real alpha(alpha_raw);
                const Real beta(beta_raw);
                cout << "alpha=" << alpha << endl;
                cout << "beta=" << beta << endl;

                const ScalarField positions_prescribed = compute_positions(mesh);

                Operator aa(3*nvertices, 3*nvertices);
                Operator bb(3*nvertices, 3*nvertices);

                // normal and position coupling
                for (const Face& face : mesh.faces)
                {
                    std::vector<int> vertex_indexes;
                    {
                        HalfEdgeCIter he = face.he;
                        do
                        {
                            vertex_indexes.push_back(he->vertex->index);
                            he = he->next;
                        }
                        while(he != face.he);
                        assert( vertex_indexes.size() == 3 );
                        assert( vertex_indexes[0] != vertex_indexes[1] );
                        assert( vertex_indexes[1] != vertex_indexes[2] );
                        assert( vertex_indexes[2] != vertex_indexes[0] );
                    }

                    //Position of vertices skipped in inpainting areas
                    for (const int& vertex_index_ii : vertex_indexes)
                    {
                        assert( vertex_index_ii < nvertices );
                        //if (inpainting_enabled && mesh.vertices[vertex_index_ii].inpainting) continue;
                        for (const int& vertex_index_jj : vertex_indexes)
                        {
                            assert( vertex_index_jj < nvertices );
                            //if (inpainting_enabled && mesh.vertices[vertex_index_jj].inpainting) continue;
                            const double scale = ( vertex_index_ii == vertex_index_jj ? 2 : -1 );
                            for (int dir_ii=0; dir_ii<3; dir_ii++)
                                for (int dir_jj=0; dir_jj<3; dir_jj++)
                                {
                                    const double normal_factor = face.normal_regularized[dir_ii]*face.normal_regularized[dir_jj];
                                    aa(vertex_index_ii+nvertices*dir_ii, vertex_index_jj+nvertices*dir_jj) += scale*normal_factor;
                                }
                        }
                    }
                }

                // fairness
                for (const Edge& edge : mesh.edges)
                {
                    const bool onBoundary = edge.he->onBoundary || edge.he->flip->onBoundary;
                    if (onBoundary) continue;
                    const int v1_index = edge.he->vertex->index;
                    const int v2_index = edge.he->next->next->vertex->index;
                    const int v3_index = edge.he->flip->vertex->index;
                    const int v4_index = edge.he->flip->next->next->vertex->index;
                    assert( v1_index != v2_index );
                    assert( v1_index != v3_index );
                    assert( v1_index != v4_index );
                    assert( v2_index != v3_index );
                    assert( v2_index != v4_index );
                    assert( v3_index != v4_index );
                    const double feature_scalar = (mesh.vertices[v1_index].feature_scalar+mesh.vertices[v3_index].feature_scalar)/2;
                    const double beta_vv = beta*feature_scalar*feature_scalar;
                    for (int dir=0; dir<3; dir++)
                    {
                        aa(v1_index+nvertices*dir, v1_index+nvertices*dir) += beta_vv;
                        aa(v1_index+nvertices*dir, v2_index+nvertices*dir) -= beta_vv;
                        aa(v1_index+nvertices*dir, v3_index+nvertices*dir) += beta_vv;
                        aa(v1_index+nvertices*dir, v4_index+nvertices*dir) -= beta_vv;

                        aa(v3_index+nvertices*dir, v1_index+nvertices*dir) += beta_vv;
                        aa(v3_index+nvertices*dir, v2_index+nvertices*dir) -= beta_vv;
                        aa(v3_index+nvertices*dir, v3_index+nvertices*dir) += beta_vv;
                        aa(v3_index+nvertices*dir, v4_index+nvertices*dir) -= beta_vv;

                        aa(v2_index+nvertices*dir, v1_index+nvertices*dir) -= beta_vv;
                        aa(v2_index+nvertices*dir, v2_index+nvertices*dir) += beta_vv;
                        aa(v2_index+nvertices*dir, v3_index+nvertices*dir) -= beta_vv;
                        aa(v2_index+nvertices*dir, v4_index+nvertices*dir) += beta_vv;

                        aa(v4_index+nvertices*dir, v1_index+nvertices*dir) -= beta_vv;
                        aa(v4_index+nvertices*dir, v2_index+nvertices*dir) += beta_vv;
                        aa(v4_index+nvertices*dir, v3_index+nvertices*dir) -= beta_vv;
                        aa(v4_index+nvertices*dir, v4_index+nvertices*dir) += beta_vv;
                    }
                }

                // regularization
                // We skip the vertex is in "inpainting" areas
                for (int index=0; index<3*nvertices; index++)
                {
                    const bool painted = mesh.vertices[index%nvertices].inpainting;
                    if (inpainting_enabled && painted) continue;
                    aa(index, index) += alpha;
                    bb(index, index) += alpha;
                }

                // solve
                ScalarField positions_regularized = positions_prescribed;
                {
                    ScalarField right_hand_term = bb*positions_prescribed;
                    solvePositiveDefinite(aa, positions_regularized, right_hand_term);
                }

                store_regularized_positions(positions_regularized, mesh);

                cout << "FINISHED POSITIONS REGULARIZATION!!!!!!!!!" << endl << endl;
            }

            static void run_version_a(const double alpha_normal, const double lambda_normal, const double epsilon_start, const double epsilon_finish, const double epsilon_progression, const double alpha_pos, const double beta_pos, const bool normal_inpainting_enabled, const bool position_inpainting_enabled, const bool exclusion_enabled, Mesh& mesh)
            {
                cout << "###### VERSION A BEGIN ######" << endl;

                for (int iter=0; iter<10; iter++)
                {
                    run_normals(alpha_normal, lambda_normal, epsilon_start, epsilon_finish, epsilon_progression, normal_inpainting_enabled, exclusion_enabled, mesh);
                    run_positions(alpha_pos, beta_pos, position_inpainting_enabled, mesh);

                    double mean_displacement = 0;
                    for (Vertex& vertex : mesh.vertices)
                    {
                        mean_displacement += (vertex.position_regularized-vertex.position).norm();
                        vertex.position = vertex.position_regularized;
                    }
                    mean_displacement /= mesh.vertices.size();
                    cout << "mean_displacement=" << mean_displacement << endl;

                    if (mean_displacement < 2e-3) break;
                }

                cout << "###### VERSION A END ######" << endl;
            }

            static void run_version_b(const double alpha_normal, const double lambda_normal, const double epsilon_start, const double epsilon_finish, const double epsilon_progression, const double alpha_pos, const double beta_pos, const bool normal_inpainting_enabled, const bool position_inpainting_enabled, const bool exclusion_enabled, Mesh& mesh)
            { // at normals
                cout << "###### VERSION B BEGIN ######" << endl;

                cout << "////////////////" << endl;
                cout << "normal_inpainting_enabled=" << normal_inpainting_enabled << endl;
                cout << "position_inpainting_enabled=" << position_inpainting_enabled << endl;
                cout << "exclusion_enabled=" << exclusion_enabled << endl;
                cout << "alpha_normal=" << alpha_normal << endl;
                cout << "lambda_normal=" << lambda_normal << endl;
                cout << "epsilon_start=" << epsilon_start << endl;
                cout << "epsilon_finish=" << epsilon_finish << endl;
                cout << "epsilon_progression=" << epsilon_progression << endl;
                cout << "////////////////" << endl;

                const size_t nfaces = mesh.faces.size();
                const size_t nedges = mesh.edges.size();
                const size_t nvertices = mesh.vertices.size();
                const double total_area = mesh.area();
                const double mean_edge_length = mesh.meanEdgeLength();
                //const double mean_diedral_angle = mesh.meanDiedralAngle();
                cout << "nfaces=" << nfaces << endl;
                cout << "nedges=" << nedges << endl;
                cout << "nvertices=" << nvertices << endl;
                cout << "total_area=" << total_area << endl;
                cout << "mean_edge_length=" << mean_edge_length << endl;
                //cout << "mean_diedral_angle=" << mean_diedral_angle << endl;

                cout << endl << "INITIALIZE DEC OPERATORS" << endl;

                Operator d0;
                ExteriorDerivative0Form<Real>::build( mesh, d0 );
                cout << "d0=" << d0.nRows() << "x" << d0.nColumns() << endl;

                Operator ad1;
                {
                    Operator d1;
                    ExteriorDerivative1Form<Real>::build( mesh, d1 );
                    ad1 = d1.transpose();
                }
                cout << "ad1=" << ad1.nRows() << "x" << ad1.nColumns() << endl;

                Operator mm;
                Wedge0FormTo1Form<Real>::build( mesh, mm );
                cout << "mm=" << mm.nRows() << "x" << mm.nColumns() << endl;

                Operator sp0;
                HodgeStar0Form<Real>::build( mesh, false, sp0 );
                for (int index=0; index<sp0.nRows(); index++) sp0(index, index) += total_area*1e-3/mesh.vertices.size();
                cout << "sp0=" << sp0.nRows() << "x" << sp0.nColumns() << endl;

                Operator sp1;
                HodgeStar1Form<Real>::build( mesh, false, sp1 );
                sp1 *= .5;
                cout << "sp1=" << sp1.nRows() << "x" << sp1.nColumns() << endl;

                Operator sp0p;
                HodgeStar2Form<Real>::build( mesh, true, sp0p );
                for (int index=0; index<sp0p.nRows(); index++) sp0p(index, index) += total_area*1e-3/mesh.faces.size();
                cout << "sp0p=" << sp0p.nRows() << "x" << sp0p.nColumns() << endl;

                Operator sp1p;
                HodgeStar1Form<Real>::build( mesh, true, sp1p );
                sp1p *= .5;
                cout << "sp1p=" << sp1p.nRows() << "x" << sp1p.nColumns() << endl;

                cout << endl << "ENTERING MAIN LOOP" << endl;

                ScalarField feature_exclusion_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    feature_exclusion_scalars(vertex.index) = exclusion_enabled ? (exclusion_enabled ? 1e-1 : 1) : 1;

                ScalarField inpainting_scalars(nfaces);
                for (const Face& face : mesh.faces)
                {
                    bool painted = face.he->vertex->inpainting || face.he->next->vertex->inpainting || face.he->next->next->vertex->inpainting;
                    inpainting_scalars(face.index) = normal_inpainting_enabled ? ( painted ? 1e-8 : 1 ) : 1;
                }
                ScalarField one_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    one_scalars(vertex.index) = 1;
                ScalarField feature_scalars = one_scalars;
                cout << "feature_scalars=" << feature_scalars.minCoeff() << "/" << feature_scalars.maxCoeff() << endl;

                const Real alpha(alpha_normal/mean_edge_length/mean_edge_length);
                const Real lambda(lambda_normal/mean_edge_length);
                cout << "alpha=" << alpha << endl;
                cout << "lambda=" << lambda << endl;

                ScalarField feature_scalars_last = feature_scalars;
                for (double epsilon_normal=epsilon_start; epsilon_normal>=epsilon_finish; epsilon_normal/=epsilon_progression)
                {
                    cout << "---------------------- epsilon_normal=" << epsilon_normal << endl;
                    const Real epsilon(epsilon_normal*mean_edge_length);

                    const VectorField normals_prescribed = compute_normals(mesh);
                    VectorField normals_estimated = normals_prescribed;

                    for (int iteration=0; iteration<50; iteration++)
                    {
                        cout << "*********** " << iteration << endl;

                        // solving normal field
                        solve_normals(alpha, sp0p, sp1p, ad1, mm, feature_scalars, inpainting_scalars, normals_prescribed, normals_estimated);

                        // solving feature field
                        solve_feature_scalars(lambda, epsilon, sp0, sp1, sp1p, d0, ad1, mm, one_scalars, normals_estimated, feature_scalars, feature_exclusion_scalars);

                        const double feature_scalars_variation_norm = (feature_scalars-feature_scalars_last).norm();
                        cout << "variation_norm=" << feature_scalars_variation_norm << endl;
                        const double feature_scalars_rel_variation_norm = feature_scalars_variation_norm/feature_scalars.norm();
                        cout << "rel_variation_norm=" << feature_scalars_rel_variation_norm << endl;
                        feature_scalars_last = feature_scalars;

                        if (feature_scalars_rel_variation_norm < 1e-2) break;

                    }

                    //normalize_normals(normals_estimated);
                    store_feature_scalars(feature_scalars, mesh);
                    //normalize_feature_scalars(feature_scalars);
                    store_regularized_normals(normals_estimated, mesh);

                    run_positions(alpha_pos, beta_pos, position_inpainting_enabled, mesh);
                }

                cout << "FINISHED MAIN LOOP!!!!!!!!!" << endl << endl;

                cout << "###### VERSION B END ######" << endl;
            }

            static void run_version_c(const double alpha_normal, const double lambda_normal, const double epsilon_start, const double epsilon_finish, const double epsilon_progression, const double alpha_pos, const double beta_pos, const bool normal_inpainting_enabled, const bool position_inpainting_enabled, const bool exclusion_enabled, Mesh& mesh)
            { // at normals
                cout << "###### VERSION C BEGIN ######" << endl;

                cout << "////////////////" << endl;
                cout << "normal_inpainting_enabled=" << normal_inpainting_enabled << endl;
                cout << "position_inpainting_enabled=" << position_inpainting_enabled << endl;
                cout << "exclusion_enabled=" << exclusion_enabled << endl;
                cout << "alpha_normal=" << alpha_normal << endl;
                cout << "lambda_normal=" << lambda_normal << endl;
                cout << "epsilon_start=" << epsilon_start << endl;
                cout << "epsilon_finish=" << epsilon_finish << endl;
                cout << "epsilon_progression=" << epsilon_progression << endl;
                cout << "////////////////" << endl;

                const size_t nfaces = mesh.faces.size();
                const size_t nedges = mesh.edges.size();
                const size_t nvertices = mesh.vertices.size();
                const double total_area = mesh.area();
                const double mean_edge_length = mesh.meanEdgeLength();
                //const double mean_diedral_angle = mesh.meanDiedralAngle();
                cout << "nfaces=" << nfaces << endl;
                cout << "nedges=" << nedges << endl;
                cout << "nvertices=" << nvertices << endl;
                cout << "total_area=" << total_area << endl;
                cout << "mean_edge_length=" << mean_edge_length << endl;
                //cout << "mean_diedral_angle=" << mean_diedral_angle << endl;

                cout << endl << "INITIALIZE DEC OPERATORS" << endl;

                Operator d0;
                ExteriorDerivative0Form<Real>::build( mesh, d0 );
                cout << "d0=" << d0.nRows() << "x" << d0.nColumns() << endl;

                Operator ad1;
                {
                    Operator d1;
                    ExteriorDerivative1Form<Real>::build( mesh, d1 );
                    ad1 = d1.transpose();
                }
                cout << "ad1=" << ad1.nRows() << "x" << ad1.nColumns() << endl;

                Operator mm;
                Wedge0FormTo1Form<Real>::build( mesh, mm );
                cout << "mm=" << mm.nRows() << "x" << mm.nColumns() << endl;

                Operator sp0;
                HodgeStar0Form<Real>::build( mesh, false, sp0 );
                for (int index=0; index<sp0.nRows(); index++) sp0(index, index) += total_area*1e-3/mesh.vertices.size();
                cout << "sp0=" << sp0.nRows() << "x" << sp0.nColumns() << endl;

                Operator sp1;
                HodgeStar1Form<Real>::build( mesh, false, sp1 );
                sp1 *= .5;
                cout << "sp1=" << sp1.nRows() << "x" << sp1.nColumns() << endl;

                Operator sp0p;
                HodgeStar2Form<Real>::build( mesh, true, sp0p );
                for (int index=0; index<sp0p.nRows(); index++) sp0p(index, index) += total_area*1e-3/mesh.faces.size();
                cout << "sp0p=" << sp0p.nRows() << "x" << sp0p.nColumns() << endl;

                Operator sp1p;
                HodgeStar1Form<Real>::build( mesh, true, sp1p );
                sp1p *= .5;
                cout << "sp1p=" << sp1p.nRows() << "x" << sp1p.nColumns() << endl;

                cout << endl << "ENTERING MAIN LOOP" << endl;

                ScalarField feature_exclusion_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    feature_exclusion_scalars(vertex.index) = exclusion_enabled ? (exclusion_enabled ? 1e-1 : 1) : 1;

                ScalarField inpainting_scalars(nfaces);
                for (const Face& face : mesh.faces)
                {
                    bool painted = face.he->vertex->inpainting || face.he->next->vertex->inpainting || face.he->next->next->vertex->inpainting;
                    inpainting_scalars(face.index) = normal_inpainting_enabled ? ( painted ? 1e-3 : 1 ) : 1;
                }
                ScalarField one_scalars(nvertices);
                for (const Vertex& vertex : mesh.vertices)
                    one_scalars(vertex.index) = 1;
                ScalarField feature_scalars = one_scalars;
                cout << "feature_scalars=" << feature_scalars.minCoeff() << "/" << feature_scalars.maxCoeff() << endl;

                const Real alpha(alpha_normal/mean_edge_length/mean_edge_length);
                const Real lambda(lambda_normal/mean_edge_length);
                cout << "alpha=" << alpha << endl;
                cout << "lambda=" << lambda << endl;

                ScalarField feature_scalars_last = feature_scalars;
                for (double epsilon_normal=epsilon_start; epsilon_normal>=epsilon_finish; epsilon_normal/=epsilon_progression)
                {
                    cout << "---------------------- epsilon_normal=" << epsilon_normal << endl;
                    const Real epsilon(epsilon_normal*mean_edge_length);

                    for (int iteration=0; iteration<50; iteration++)
                    {
                        cout << "*********** " << iteration << endl;

                        // solving normal field
                        const VectorField normals_prescribed = compute_normals(mesh);
                        VectorField normals_estimated = normals_prescribed;
                        solve_normals(alpha, sp0p, sp1p, ad1, mm, feature_scalars, inpainting_scalars, normals_prescribed, normals_estimated);

                        // solving feature field
                        solve_feature_scalars(lambda, epsilon, sp0, sp1, sp1p, d0, ad1, mm, one_scalars, normals_estimated, feature_scalars, feature_exclusion_scalars);

                        const double feature_scalars_variation_norm = (feature_scalars-feature_scalars_last).norm();
                        cout << "variation_norm=" << feature_scalars_variation_norm << endl;
                        const double feature_scalars_rel_variation_norm = feature_scalars_variation_norm/feature_scalars.norm();
                        cout << "rel_variation_norm=" << feature_scalars_rel_variation_norm << endl;
                        feature_scalars_last = feature_scalars;

                        //normalize_normals(normals_estimated);
                        store_feature_scalars(feature_scalars, mesh);
                        //normalize_feature_scalars(feature_scalars);
                        store_regularized_normals(normals_estimated, mesh);

                        run_positions(alpha_pos, beta_pos, position_inpainting_enabled, mesh);

                        if (feature_scalars_rel_variation_norm < 1e-2) break;
                    }

                }

                cout << "FINISHED MAIN LOOP!!!!!!!!!" << endl << endl;

                cout << "###### VERSION C END ######" << endl;
            }

            static std::array<double, 4> compute_cot_weights(const Edge& edge)
            {
                std::array<double, 4> ww;
                ww[0] = -(edge.he->flip->next->cotan() + edge.he->next->next->cotan());
                ww[1] = edge.he->next->next->cotan() + edge.he->next->cotan();
                ww[2] = -(edge.he->next->cotan() + edge.he->flip->next->next->cotan());
                ww[3] = edge.he->flip->next->next->cotan() + edge.he->flip->next->cotan();
                return ww;
            }

            static std::array<double, 4> compute_area_weights(const Edge& edge, const Vector& v1, const Vector& v2, const Vector& v3, const Vector& v4)
            {
                std::array<double, 4> ww;
                const double area_123 = edge.he->face->area();
                const double area_134 = edge.he->flip->face->area();
                const double area_total = area_123+area_134;
                ww[0] = (area_123*dot(v4-v3, v3-v1) + area_134*dot(v1-v3, v3-v2))/dot(v3-v1, v3-v1)/area_total;
                ww[1] = area_134/area_total;
                ww[2] = (area_123*dot(v3-v1, v1-v4) + area_134*dot(v2-v1, v1-v3))/dot(v3-v1, v3-v1)/area_total;
                ww[3] = area_123/area_total;
                return ww;
            }

            static void run_heschaefer(const double alpha_he, const double lambda_he, const double beta_start, const double beta_finish, const double beta_progression, const bool use_area, Mesh& mesh)
            {
                cout << "###### HE-SCHAEFER L0 BEGIN ######" << endl;
                const double mean_edge_length = mesh.meanEdgeLength();
                const double mean_diedral_angle = mesh.meanDiedralAngle();

                const double alpha_start = alpha_he*mean_edge_length;
                const double lambda = lambda_he*mean_edge_length*mean_edge_length*mean_diedral_angle;

                cout << "////////////////" << endl;
                cout << "mean_edge_length=" << mean_edge_length << endl;
                cout << "mean_diedral_angle=" << mean_diedral_angle << endl;
                cout << "////////////////" << endl;
                cout << "alpha_start=" << alpha_start << endl;
                cout << "use_area=" << use_area << endl;
                cout << "lambda=" << lambda << endl;
                cout << "beta_start=" << beta_start << endl;
                cout << "beta_finish=" << beta_finish << endl;
                cout << "beta_progression=" << beta_progression << endl;
                cout << "////////////////" << endl;
                cout << endl;

                const size_t nfaces = mesh.faces.size();
                const size_t nedges = mesh.edges.size();
                const size_t nvertices = mesh.vertices.size();
                cout << "nfaces=" << nfaces << endl;
                cout << "nedges=" << nedges << endl;
                cout << "nvertices=" << nvertices << endl;

                const ScalarField positions_prescribed = compute_positions(mesh);
                ScalarField positions_regularized = positions_prescribed;
                ScalarField delta(3*nedges); delta.zero();

                Operator rr(3*nedges, 3*nvertices);
                for (const Edge& edge : mesh.edges)
                {
                    const bool onBoundary = edge.he->onBoundary || edge.he->flip->onBoundary;
                    assert( !onBoundary );

                    const int edge_index = edge.index;
                    const int v1_index = edge.he->vertex->index;
                    const int v2_index = edge.he->next->next->vertex->index;
                    const int v3_index = edge.he->flip->vertex->index;
                    const int v4_index = edge.he->flip->next->next->vertex->index;
                    assert( v1_index != v2_index );
                    assert( edge.he->next->vertex->index == v3_index );
                    assert( edge.he->flip->next->vertex->index == v1_index );
                    assert( v1_index != v3_index );
                    assert( v1_index != v4_index );
                    assert( v2_index != v3_index );
                    assert( v2_index != v4_index );
                    assert( v3_index != v4_index );

                    { // fairness
                        for (int dir=0; dir<3; dir++)
                        {
                            rr(edge_index+dir*nedges, v1_index+dir*nvertices) = 1;
                            rr(edge_index+dir*nedges, v2_index+dir*nvertices) = -1;
                            rr(edge_index+dir*nedges, v3_index+dir*nvertices) = 1;
                            rr(edge_index+dir*nedges, v4_index+dir*nvertices) = -1;
                        }
                    }
                }

                Operator idpos(3*nvertices, 3*nvertices);
                for (int index=0; index<3*nvertices; index++)
                    idpos(index, index) = 1;

                Real alpha = alpha_start;
                for (Real beta=beta_start; beta<beta_finish; beta*=beta_progression)
                {
                    cout << "---------------------- beta=" << beta << " epsilon=" << 1/beta << endl;

                    // constructing operators
                    Operator dd(3*nedges, 3*nvertices);
                    Operator rr(3*nedges, 3*nvertices);
                    for (const Edge& edge : mesh.edges)
                    {
                        const bool onBoundary = edge.he->onBoundary || edge.he->flip->onBoundary;
                        assert( !onBoundary );

                        const int edge_index = edge.index;
                        const int v1_index = edge.he->vertex->index;
                        const int v2_index = edge.he->next->next->vertex->index;
                        const int v3_index = edge.he->flip->vertex->index;
                        const int v4_index = edge.he->flip->next->next->vertex->index;
                        assert( v1_index != v2_index );
                        assert( edge.he->next->vertex->index == v3_index );
                        assert( edge.he->flip->next->vertex->index == v1_index );
                        assert( v1_index != v3_index );
                        assert( v1_index != v4_index );
                        assert( v2_index != v3_index );
                        assert( v2_index != v4_index );
                        assert( v3_index != v4_index );

                        { // differential
                            const std::array<double, 4> ww = use_area ? compute_area_weights(edge,
                                mesh.vertices[v1_index].position,
                                mesh.vertices[v2_index].position,
                                mesh.vertices[v3_index].position,
                                mesh.vertices[v4_index].position) : compute_cot_weights(edge);

                            for (int dir=0; dir<3; dir++)
                            {
                                dd(edge_index+dir*nedges, v1_index+dir*nvertices) = ww[0];
                                dd(edge_index+dir*nedges, v2_index+dir*nvertices) = ww[1];
                                dd(edge_index+dir*nedges, v3_index+dir*nvertices) = ww[2];
                                dd(edge_index+dir*nedges, v4_index+dir*nvertices) = ww[3];
                            }
                        }
                    }

                    // solving delta
                    {
                        const double threshold = sqrt(lambda/beta);
                        cout << "threshold=" << threshold << endl;

                        const ScalarField grad_positions = dd*positions_regularized;
                        int delta_norm_0 = 0;
                        delta.zero();
                        for (int edge_index=0; edge_index<nedges; edge_index++)
                        {
                            const Vector grad_position(grad_positions(edge_index), grad_positions(edge_index+nedges), grad_positions(edge_index+2*nedges));
                            const double grad_position_norm = grad_position.norm();

                            if (grad_position_norm < threshold) continue;

                            delta_norm_0 += 1;
                            for (int dir=0; dir<3; dir++) delta(edge_index+dir*nedges) = grad_positions(edge_index+dir*nedges);
                        }
                        cout << "delta_norm_0=" << delta_norm_0 << "/" << nedges << endl;
                    }

                    // solving positions
                    Operator op = idpos + alpha*rr.transpose()*rr + beta*dd.transpose()*dd;
                    ScalarField rhs = positions_prescribed + beta*dd.transpose()*delta;
                    solvePositiveDefinite(op, positions_regularized, rhs);

                    alpha /= 2;
                }

                for (Vertex& vertex : mesh.vertices)
                    vertex.feature_scalar = 1;

                for (Edge& edge : mesh.edges)
                {
                    const Vector delta_edge(delta(edge.index), delta(edge.index+nedges), delta(edge.index+2*nedges));
                    const double delta_edge_norm = delta_edge.norm();

                    edge.scalar = 1;
                    if (delta_edge_norm < 1e-3) continue;

                    edge.scalar = 0;
                    edge.he->vertex->feature_scalar -= .5;
                    edge.he->flip->vertex->feature_scalar -= .5;
                }

                store_regularized_positions(positions_regularized, mesh);

                for (Face& face : mesh.faces)
                {
                    const Vector p0 = face.he->vertex->position_regularized;
                    const Vector p1 = face.he->next->vertex->position_regularized;
                    const Vector p2 = face.he->next->next->vertex->position_regularized;

                    face.normal_regularized = cross( p1-p0, p2-p0 ).unit();
                }

                cout << "###### HE-SCHAEFER L0 END ######" << endl;
            }

            static void run_debug(const Mesh& mesh)
            {
                cout << "area=" << mesh.area() << endl;

                { // primal 0-forms
                    SparseMatrix<Real> sp0;
                    DDG::HodgeStar0Form<Real>::build(mesh, false, sp0);

                    DenseMatrix<Real> ff(mesh.vertices.size());
                    ff.zero(1);

                    const SparseMatrix<Real> scalar_product = ff.transpose().sparse()*sp0*ff.sparse();
                    assert( scalar_product.nRows() == 1 && scalar_product.nColumns() == 1 );
                    cout << "primal_0form_area=" << scalar_product(0,0) << endl;
                }
                { // dual 0-forms
                    SparseMatrix<Real> sp0p;
                    DDG::HodgeStar2Form<Real>::build(mesh, true, sp0p);

                    DenseMatrix<Real> ff(mesh.faces.size());
                    ff.zero(1);

                    const SparseMatrix<Real> scalar_product = ff.transpose().sparse()*sp0p*ff.sparse();
                    assert( scalar_product.nRows() == 1 && scalar_product.nColumns() == 1 );
                    cout << "dual_0form_area=" << scalar_product(0,0) << endl;
                }
                { // primal 1-forms
                    SparseMatrix<Real> sp1;
                    DDG::HodgeStar1Form<Real>::build(mesh, false, sp1);
                    sp1 *= .5;

                    DenseMatrix<Real> ff(mesh.edges.size());
                    for (int index=0; index<ff.nRows(); index++)
                    {
                        const Edge& edge = mesh.edges[index];
                        const double length = (edge.he->vertex->position - edge.he->flip->vertex->position).norm();
                        ff(index) = length;
                    }

                    const SparseMatrix<Real> scalar_product = ff.transpose().sparse()*sp1*ff.sparse();
                    assert( scalar_product.nRows() == 1 && scalar_product.nColumns() == 1 );
                    cout << "primal_1form_area=" << scalar_product(0,0) << endl;
                }
            }

    };

}
#endif
