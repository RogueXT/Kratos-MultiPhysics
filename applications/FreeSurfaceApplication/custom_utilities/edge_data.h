//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//

#if !defined(KRATOS_EDGE_DATA_H_INCLUDED)
#define KRATOS_EDGE_DATA_H_INCLUDED

// we suggest defining the following macro
#define USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

// we suggest defining the following macro
#define USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "free_surface_application.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
    template <unsigned int TDim>
    class EdgesStructureType
    {
    public:
        // component ij of the consistent mass matrix (M = Ni * Nj * dOmega)
        double Mass = 0.0;
        // components kl of the laplacian matrix of edge ij (L = dNi/dxk * dNj/dxl * dOmega)
        // double Laplacian;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> LaplacianIJ;
        // components k of the gradient matrix of edge ij (G = Ni * dNj/dxl * dOmega)
        array_1d<double, TDim> Ni_DNj = ZeroVector(TDim);
        // components k of the transposed gradient matrix of edge ij (GT = dNi/dxl * Nj * dOmega)
        // TRANSPOSED GRADIENT
        array_1d<double, TDim> DNi_Nj = ZeroVector(TDim);

        //*************************************************************************************
        //*************************************************************************************
        // gradient integrated by parts
        // RHSi += DNi_Nj pj + Aboundary * pext ==> RHS += Ni_DNj p_j - DNi_Nj p_i
        // ATTENTION: + Aboundary * pext is NOT included!! it should be included "manually"

        inline void Add_Gp(array_1d<double, TDim> &destination, const double &p_i, const double &p_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination[comp] -= Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i;
        }

        inline void Sub_Gp(array_1d<double, TDim> &destination, const double &p_i, const double &p_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination[comp] += Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i;
        }

        //*************************************************************************************
        //*************************************************************************************
        // gradient
        // RHSi += Ni_DNj[k]*v[k]

        inline void Add_D_v(double &destination,
                            const array_1d<double, TDim> &v_i,
                            const array_1d<double, TDim> &v_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination += Ni_DNj[comp] * (v_j[comp] - v_i[comp]);
        }

        inline void Sub_D_v(double &destination,
                            const array_1d<double, TDim> &v_i,
                            const array_1d<double, TDim> &v_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination -= Ni_DNj[comp] * (v_j[comp] - v_i[comp]);
        }

        //*************************************************************************************
        //*************************************************************************************
        // gradient
        // RHSi += Ni_DNj pj

        inline void Add_grad_p(array_1d<double, TDim> &destination, const double &p_i, const double &p_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination[comp] += Ni_DNj[comp] * (p_j - p_i);
        }

        inline void Sub_grad_p(array_1d<double, TDim> &destination, const double &p_i, const double &p_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination[comp] -= Ni_DNj[comp] * (p_j - p_i);
        }

        //*************************************************************************************
        //*************************************************************************************
        // gradient
        // RHSi += DNi_Nj[k]*v[k]

        inline void Add_div_v(double &destination,
                              const array_1d<double, TDim> &v_i,
                              const array_1d<double, TDim> &v_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination -= Ni_DNj[comp] * v_j[comp] - DNi_Nj[comp] * v_i[comp];
        }

        inline void Sub_div_v(double &destination,
                              const array_1d<double, TDim> &v_i,
                              const array_1d<double, TDim> &v_j)
        {
            for (unsigned int comp = 0; comp < TDim; comp++)
                destination += Ni_DNj[comp] * v_j[comp] - DNi_Nj[comp] * v_i[comp];
        }

        //*************************************************************************************
        //*************************************************************************************
        // gets the trace of the laplacian matrix

        inline void CalculateScalarLaplacian(double &l_ij)
        {
            l_ij = LaplacianIJ(0, 0);
            for (unsigned int comp = 1; comp < TDim; comp++)
                l_ij += LaplacianIJ(comp, comp);
        }

        inline void Add_ConvectiveContribution(array_1d<double, TDim> &destination,
                                               const array_1d<double, TDim> &a_i, const array_1d<double, TDim> &U_i,
                                               const array_1d<double, TDim> &a_j, const array_1d<double, TDim> &U_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION
            double temp = a_i[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] += temp * (U_j[l_comp] - U_i[l_comp]);
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] += aux_j * U_j[l_comp] - aux_i * U_i[l_comp];
#endif
        }

        inline void Sub_ConvectiveContribution(array_1d<double, TDim> &destination,
                                               const array_1d<double, TDim> &a_i, const array_1d<double, TDim> &U_i,
                                               const array_1d<double, TDim> &a_j, const array_1d<double, TDim> &U_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION
            double temp = a_i[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] -= temp * (U_j[l_comp] - U_i[l_comp]);
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] -= aux_j * U_j[l_comp] - aux_i * U_i[l_comp];
#endif
        }

        inline void Sub_ConvectiveContribution(double &destination,
                                               const array_1d<double, TDim> &a_i, const double &phi_i,
                                               const array_1d<double, TDim> &a_j, const double &phi_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION
            double temp = a_i[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];

            destination -= temp * (phi_j - phi_i);
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }
            destination -= aux_j * phi_j - aux_i * phi_i;
#endif
        }

        inline void Add_ConvectiveContribution(double &destination,
                                               const array_1d<double, TDim> &a_i, const double &phi_i,
                                               const array_1d<double, TDim> &a_j, const double &phi_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION
            double temp = a_i[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];

            destination += temp * (phi_j - phi_i);
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }
            destination += aux_j * phi_j - aux_i * phi_i;
#endif
        }

        //*************************************************************************************
        //*************************************************************************************

        inline void CalculateConvectionStabilization_LOW(array_1d<double, TDim> &stab_low,
                                                         const array_1d<double, TDim> &a_i, const array_1d<double, TDim> &U_i,
                                                         const array_1d<double, TDim> &a_j, const array_1d<double, TDim> &U_j)
        {
            double conv_stab = 0.0;
            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
                    conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp, m_comp);
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                stab_low[l_comp] = conv_stab * (U_j[l_comp] - U_i[l_comp]);
        }

        inline void CalculateConvectionStabilization_LOW(double &stab_low,
                                                         const array_1d<double, TDim> &a_i, const double &phi_i,
                                                         const array_1d<double, TDim> &a_j, const double &phi_j)
        {
            double conv_stab = 0.0;
            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
                    conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp, m_comp);
            stab_low = conv_stab * (phi_j - phi_i);
        }
        //*************************************************************************************
        //*************************************************************************************

        inline void CalculateConvectionStabilization_HIGH(array_1d<double, TDim> &stab_high,
                                                          const array_1d<double, TDim> &a_i, const array_1d<double, TDim> &pi_i,
                                                          const array_1d<double, TDim> &a_j, const array_1d<double, TDim> &pi_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION
            double temp = 0.0;
            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                stab_high[l_comp] = -temp * (pi_j[l_comp] - pi_i[l_comp]); // check if the minus sign is correct
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                stab_high[l_comp] = -(aux_j * pi_j[l_comp] - aux_i * pi_i[l_comp]);
#endif
        }

        inline void CalculateConvectionStabilization_HIGH(double &stab_high,
                                                          const array_1d<double, TDim> &a_i, const double &pi_i,
                                                          const array_1d<double, TDim> &a_j, const double &pi_j)
        {

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION
            double temp = 0.0;
            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                temp += a_i[k_comp] * Ni_DNj[k_comp];

            stab_high = -temp * (pi_j - pi_i); // check if the minus sign is correct
#else
            double aux_i = a_i[0] * Ni_DNj[0];
            double aux_j = a_j[0] * Ni_DNj[0];
            for (unsigned int k_comp = 1; k_comp < TDim; k_comp++)
            {
                aux_i += a_i[k_comp] * Ni_DNj[k_comp];
                aux_j += a_j[k_comp] * Ni_DNj[k_comp];
            }

            stab_high = -(aux_j * pi_j - aux_i * pi_i);
#endif
        }
        //*************************************************************************************
        //*************************************************************************************

        inline void Add_StabContribution(array_1d<double, TDim> &destination,
                                         const double tau, const double beta,
                                         const array_1d<double, TDim> &stab_low, const array_1d<double, TDim> &stab_high)
        {
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] += tau * (stab_low[l_comp] - beta * stab_high[l_comp]);
        }

        inline void Add_StabContribution(double &destination,
                                         const double tau, const double beta,
                                         const double &stab_low, const double &stab_high)
        {
            destination += tau * (stab_low - beta * stab_high);
        }

        inline void Sub_StabContribution(array_1d<double, TDim> &destination,
                                         const double tau, const double beta,
                                         const array_1d<double, TDim> &stab_low, const array_1d<double, TDim> &stab_high)
        {
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] -= tau * (stab_low[l_comp] - beta * stab_high[l_comp]);
        }

        inline void Sub_StabContribution(double &destination,
                                         const double tau, const double beta,
                                         const double &stab_low, const double &stab_high)
        {
            destination -= tau * (stab_low - beta * stab_high);
        }

        //*************************************************************************************
        //*************************************************************************************

        inline void Add_ViscousContribution(array_1d<double, TDim> &destination,
                                            const array_1d<double, TDim> &U_i, const double &nu_i,
                                            const array_1d<double, TDim> &U_j, const double &nu_j)
        {
            // calculate scalar laplacian
            double L = 0.0;
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                L += LaplacianIJ(l_comp, l_comp);

            // double nu_avg = 0.5*(nu_i+nu_j);
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] += nu_i * L * (U_j[l_comp] - U_i[l_comp]);
        }

        inline void Sub_ViscousContribution(array_1d<double, TDim> &destination,
                                            const array_1d<double, TDim> &U_i, const double &nu_i,
                                            const array_1d<double, TDim> &U_j, const double &nu_j)
        {
            // calculate scalar laplacian
            double L = 0.0;
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                L += LaplacianIJ(l_comp, l_comp);

            // double nu_avg = 0.5*(nu_i+nu_j);
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                destination[l_comp] -= nu_i * L * (U_j[l_comp] - U_i[l_comp]);
        }
    };

    // class definition of matrices using CSR format

    template <unsigned int TDim, class TSparseSpace>
    class MatrixContainer
    {
    public:
        // name for the self defined structure
        typedef EdgesStructureType<TDim> CSR_Tuple;
        typedef vector<CSR_Tuple> EdgesVectorType;
        // name for row start and column index vectors
        typedef vector<unsigned int> IndicesVectorType;
        // names for separately stored node based values
        typedef vector<double> ValuesVectorType;
        typedef vector<array_1d<double, TDim>> CalcVectorType;

        // constructor and destructor

        MatrixContainer(){};

        ~MatrixContainer(){};

        // functions to return private values

        inline unsigned int GetNumberEdges()
        {
            return mNumberEdges;
        }

        inline EdgesVectorType &GetEdgeValues()
        {
            return mNonzeroEdgeValues;
        }

        inline IndicesVectorType &GetColumnIndex()
        {
            return mColumnIndex;
        }

        inline IndicesVectorType &GetRowStartIndex()
        {
            return mRowStartIndex;
        }

        inline ValuesVectorType &GetLumpedMass()
        {
            return mLumpedMassMatrix;
        }

        inline ValuesVectorType &GetInvertedMass()
        {
            return mInvertedMassMatrix;
        }

        inline CalcVectorType &GetDiagGradient()
        {
            return mDiagGradientMatrix;
        }

        inline ValuesVectorType &GetHmin()
        {
            return mHmin;
        }

        //********************************************************
        // function to size and initialize the vector of CSR tuples

        void ConstructCSRVector(ModelPart &model_part)
        {
            KRATOS_TRY

            // SIZE OF CSR VECTOR

            // defining the number of nodes and edges
            int n_nodes = model_part.Nodes().size();
            // remark: no colouring algorithm is used here (symmetry is neglected)
            //         respectively edge ij is considered different from edge ji
            mNumberEdges = 0;
            // counter to assign and get global nodal index
            int i_node = 0;

            // counting the edges connecting the nodes
            for (typename ModelPart::NodesContainerType::iterator node_it = model_part.NodesBegin(); node_it != model_part.NodesEnd(); node_it++)
            {
                // counting neighbours of each node
                mNumberEdges += (node_it->GetValue(NEIGHBOUR_NODES)).size();
                // DIAGONAL TERMS
                // mNumberEdges++;

                // assigning global index to each node
                node_it->FastGetSolutionStepValue(AUX_INDEX) = static_cast<double>(i_node++);
            }
            // error message in case number of nodes does not coincide with number of indices
            if (i_node != n_nodes)
                KRATOS_WATCH("ERROR - Highest nodal index doesn't coincide with number of nodes!");

            // allocating memory for block of CSR data - setting to zero for first-touch OpenMP allocation
            mNonzeroEdgeValues.resize(mNumberEdges);
            SetToZero(mNonzeroEdgeValues);
            mColumnIndex.resize(mNumberEdges);
            SetToZero(mColumnIndex);
            mRowStartIndex.resize(n_nodes + 1);
            SetToZero(mRowStartIndex);
            mLumpedMassMatrix.resize(n_nodes);
            SetToZero(mLumpedMassMatrix);
            mInvertedMassMatrix.resize(n_nodes);
            SetToZero(mInvertedMassMatrix);
            mDiagGradientMatrix.resize(n_nodes);
            SetToZero(mDiagGradientMatrix);
            mHmin.resize(n_nodes);
            SetToZero(mHmin);

            // INITIALIZING OF THE CSR VECTOR

            // temporary variable as the row start index of a node depends on the number of neighbours of the previous one
            unsigned int row_start_temp = 0;

            int number_of_threads = ParallelUtilities::GetNumThreads();
            std::vector<int> row_partition(number_of_threads);
            OpenMPUtils::DivideInPartitions(model_part.Nodes().size(), number_of_threads, row_partition);

            for (int k = 0; k < number_of_threads; k++)
            {
#pragma omp parallel
                if (OpenMPUtils::ThisThread() == k)
                {
                    // main loop over all nodes
                    for (unsigned int aux_i = static_cast<unsigned int>(row_partition[k]); aux_i < static_cast<unsigned int>(row_partition[k + 1]); aux_i++)
                    {
                        typename ModelPart::NodesContainerType::iterator node_it = model_part.NodesBegin() + aux_i;

                        // getting the global index of the node
                        i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

                        // determining its neighbours
                        GlobalPointersVector<Node<3>> &neighb_nodes = node_it->GetValue(NEIGHBOUR_NODES);

                        // number of neighbours of node i determines row start index for the following node
                        unsigned int n_neighbours = neighb_nodes.size();

                        // reserving memory for work array
                        std::vector<unsigned int> work_array;
                        work_array.reserve(n_neighbours);
                        // DIAGONAL TERMS

                        // nested loop over the neighbouring nodes
                        for (GlobalPointersVector<Node<3>>::iterator neighb_it = neighb_nodes.begin(); neighb_it != neighb_nodes.end(); neighb_it++)
                        {
                            // getting global index of the neighbouring node
                            work_array.push_back(static_cast<unsigned int>(neighb_it->FastGetSolutionStepValue(AUX_INDEX)));
                        }
                        // reordering neighbours following their global indices
                        std::sort(work_array.begin(), work_array.end());

                        // setting current row start index
                        mRowStartIndex[i_node] = row_start_temp;
                        // nested loop over the by now ordered neighbours
                        for (unsigned int counter = 0; counter < n_neighbours; counter++)
                        {
                            // getting global index of the neighbouring node
                            unsigned int j_neighbour = work_array[counter];
                            // calculating CSR index
                            unsigned int csr_index = mRowStartIndex[i_node] + counter;

                            // saving column index j of the original matrix
                            mColumnIndex[csr_index] = j_neighbour;
                        }
                        // preparing row start index for next node
                        row_start_temp += n_neighbours;
                    }
                }
            }
            // adding last entry (necessary for abort criterion of loops)
            mRowStartIndex[n_nodes] = mNumberEdges;

            // INITIALIZING NODE BASED VALUES
            // set the heights to a huge number
            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                mHmin[i_node] = 1e10;
            });

            KRATOS_CATCH("")
        }

        //*********************************
        // function to precalculate CSR data

        void BuildCSRData(ModelPart &model_part)
        {
            KRATOS_TRY

            // PRECALCULATING CSR DATA

            // defining temporary local variables for elementwise addition
            // shape functions
            array_1d<double, TDim + 1> N;
            // shape function derivatives
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> dN_dx;
            // volume
            double volume;
            // weighting factor
            double weighting_factor = 1.0 / static_cast<double>(TDim + 1);
            // elemental matrices
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim + 1> mass_consistent;
            // boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim+1> laplacian;
            array_1d<double, TDim + 1> mass_lumped;
            // global indices of elemental nodes
            array_1d<unsigned int, TDim + 1> nodal_indices;

            array_1d<double, TDim + 1> heights;

            // loop over all elements
            for (typename ModelPart::ElementsContainerType::iterator elem_it = model_part.ElementsBegin(); elem_it != model_part.ElementsEnd(); elem_it++)
            {
                // LOCAL ELEMENTWISE CALCULATIONS

                // getting geometry data of the element
                GeometryUtils::CalculateGeometryData(elem_it->GetGeometry(), dN_dx, N, volume);

                // calculate lenght of the heights of the element
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                {
                    heights[ie_node] = dN_dx(ie_node, 0) * dN_dx(ie_node, 0);
                    for (unsigned int comp = 1; comp < TDim; comp++)
                    {
                        heights[ie_node] += dN_dx(ie_node, comp) * dN_dx(ie_node, comp);
                    }
                    heights[ie_node] = 1.0 / sqrt(heights[ie_node]);
                    // KRATOS_WATCH(heights);
                }

                // setting up elemental mass matrices
                CalculateMassMatrix(mass_consistent, volume);
                noalias(mass_lumped) = ZeroVector(TDim + 1);
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                {
                    for (unsigned int je_node = 0; je_node <= TDim; je_node++)
                    {
                        // mass_consistent(ie_node,je_node) = N(ie_node) * N(je_node) * volume;
                        mass_lumped[ie_node] += mass_consistent(ie_node, je_node);
                    }
                    // mass_lumped[ie_node] = volume * N[ie_node];
                }

                double weighted_volume = volume * weighting_factor;

                // ASSEMBLING GLOBAL DATA STRUCTURE

                // loop over the nodes of the element to determine their global indices
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                    nodal_indices[ie_node] = static_cast<unsigned int>(elem_it->GetGeometry()[ie_node].FastGetSolutionStepValue(AUX_INDEX));

                // assembling global "edge matrices" by adding local contributions
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                {
                    // check the heights and change the value if minimal is found
                    if (mHmin[nodal_indices[ie_node]] > heights[ie_node])
                        mHmin[nodal_indices[ie_node]] = heights[ie_node];

                    for (unsigned int je_node = 0; je_node <= TDim; je_node++)
                    {
                        // remark: there is no edge linking node i with itself!
                        // DIAGONAL TERMS
                        if (ie_node != je_node)
                        {
                            // calculating CSR index from global index
                            unsigned int csr_index = GetCSRIndex(nodal_indices[ie_node], nodal_indices[je_node]);

                            // assigning precalculated element data to the referring edges
                            // contribution to edge mass
                            mNonzeroEdgeValues[csr_index].Mass += mass_consistent(ie_node, je_node);

                            // contribution to edge laplacian
                            boost::numeric::ublas::bounded_matrix<double, TDim, TDim> &laplacian = mNonzeroEdgeValues[csr_index].LaplacianIJ;
                            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                                for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                                    laplacian(l_comp, k_comp) += dN_dx(ie_node, l_comp) * dN_dx(je_node, k_comp) * volume;

                            // contribution to edge gradient
                            array_1d<double, TDim> &gradient = mNonzeroEdgeValues[csr_index].Ni_DNj;
                            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                                // gradient[l_comp] += dN_dx(je_node,l_comp);
                                gradient[l_comp] += dN_dx(je_node, l_comp) * weighted_volume;
                            // TRANSPOSED GRADIENT
                            // contribution to transposed edge gradient
                            array_1d<double, TDim> &transp_gradient = mNonzeroEdgeValues[csr_index].DNi_Nj;
                            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                                transp_gradient[l_comp] += dN_dx(ie_node, l_comp) * weighted_volume;
                        }
                    }
                }

                // assembling node based vectors
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                    // diagonal of the global lumped mass matrix
                    mLumpedMassMatrix[nodal_indices[ie_node]] += mass_lumped[ie_node];
                for (unsigned int ie_node = 0; ie_node <= TDim; ie_node++)
                {
                    // diagonal of the global gradient matrix
                    array_1d<double, TDim> &gradient = mDiagGradientMatrix[nodal_indices[ie_node]];
                    for (unsigned int component = 0; component < TDim; component++)
                        // gradient[component] += dN_dx(ie_node,component);
                        gradient[component] += dN_dx(ie_node, component) * weighted_volume;
                }
            }

            // copy mass matrix to inverted mass matrix
            for (unsigned int inode = 0; inode < mLumpedMassMatrix.size(); inode++)
            {
                mInvertedMassMatrix[inode] = mLumpedMassMatrix[inode];
            }

            // calculating inverted mass matrix (this requires syncronization for MPI paraellelism
            for (unsigned int inode = 0; inode < mInvertedMassMatrix.size(); inode++)
            {
                mInvertedMassMatrix[inode] = 1.0 / mInvertedMassMatrix[inode];
            }

            KRATOS_CATCH("")
        }

        //******************************************
        // function to calculate CSR index of edge ij

        unsigned int GetCSRIndex(unsigned int NodeI, unsigned int NeighbourJ)
        {
            KRATOS_TRY

            // index indicating data position of edge ij
            unsigned int csr_index;
            // searching for coincidence of stored column index and neighbour index j
            for (csr_index = mRowStartIndex[NodeI]; csr_index != mRowStartIndex[NodeI + 1]; csr_index++)
                if (mColumnIndex[csr_index] == NeighbourJ)
                    break;

            // returning CSR index of edge ij
            return csr_index;

            KRATOS_CATCH("")
        }

        //***********************************************
        // function to get pointer to CSR tuple of edge ij

        CSR_Tuple *GetTuplePointer(unsigned int NodeI, unsigned int NeighbourJ)
        {
            KRATOS_TRY

            // index indicating data position of edge ij
            unsigned int csr_index;
            // searching for coincidence of stored column index and neighbour index j
            for (csr_index = mRowStartIndex[NodeI]; csr_index != mRowStartIndex[NodeI + 1]; csr_index++)
                if (mColumnIndex[csr_index] == NeighbourJ)
                    break;

            // returning pointer to CSR tuple of edge ij
            return &mNonzeroEdgeValues[csr_index];

            KRATOS_CATCH("")
        }

        //*******************************
        // function to free dynamic memory

        void Clear()
        {
            KRATOS_TRY

            mNonzeroEdgeValues.clear();
            mColumnIndex.clear();
            mRowStartIndex.clear();
            mInvertedMassMatrix.clear();
            mLumpedMassMatrix.clear();
            mDiagGradientMatrix.clear();
            mHmin.clear();

            KRATOS_CATCH("")
        }
        //****************************
        // functions to access database
        //(note that this is already thought for parallel;
        // for a single processor this could be done in a faster way)

        void FillCoordinatesFromDatabase(CalcVectorType &rDestination, ModelPart::NodesContainerType &rNodes)
        {

            KRATOS_TRY

            // loop over all nodes
            int n_nodes = rNodes.size();
            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

#pragma omp parallel for firstprivate(n_nodes, it_begin)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                //  // 					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));
                unsigned int i_node = i;

                // save value in the destination vector
                for (unsigned int component = 0; component < TDim; component++)
                    (rDestination[i_node])[component] = (*node_it)[component];
            }

            KRATOS_CATCH("");
        }

        //****************************
        // functions to access database
        //(note that this is already thought for parallel;
        // for a single processor this could be done in a faster way)

        void FillVectorFromDatabase(Variable<array_1d<double, 3>> &rVariable, CalcVectorType &rDestination, ModelPart::NodesContainerType &rNodes)
        {

            KRATOS_TRY

            // loop over all nodes

            int n_nodes = rNodes.size();

            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                unsigned int i_node = i;

                // get the requested value in vector form
                array_1d<double, 3> &vector = node_it->FastGetCurrentSolutionStepValue(rVariable, var_pos);
                // save value in the destination vector
                for (unsigned int component = 0; component < TDim; component++)
                    (rDestination[i_node])[component] = vector[component];
            }

            KRATOS_CATCH("");
        }

        void FillOldVectorFromDatabase(Variable<array_1d<double, 3>> &rVariable, CalcVectorType &rDestination, ModelPart::NodesContainerType &rNodes)
        {

            KRATOS_TRY

            // loop over all nodes
            int n_nodes = rNodes.size();

            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                unsigned int i_node = i;

                // get the requested value in vector form
                array_1d<double, 3> &vector = node_it->FastGetSolutionStepValue(rVariable, 1, var_pos);
                // save value in the destination vector
                for (unsigned int component = 0; component < TDim; component++)
                    (rDestination[i_node])[component] = vector[component];
            }

            KRATOS_CATCH("");
        }

        void FillScalarFromDatabase(Variable<double> &rVariable, ValuesVectorType &rDestination, ModelPart::NodesContainerType &rNodes)
        {
            KRATOS_TRY

            // loop over all nodes
            int n_nodes = rNodes.size();

            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                unsigned int i_node = i;

                // get the requested scalar value
                double &scalar = node_it->FastGetCurrentSolutionStepValue(rVariable, var_pos);
                // save value in the destination vector
                rDestination[i_node] = scalar;
            }

            KRATOS_CATCH("");
        }

        void FillOldScalarFromDatabase(Variable<double> &rVariable, ValuesVectorType &rDestination, ModelPart::NodesContainerType &rNodes)
        {
            KRATOS_TRY

            int n_nodes = rNodes.size();
            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                unsigned int i_node = i;

                // get the requested scalar value
                double &scalar = node_it->FastGetSolutionStepValue(rVariable, 1, var_pos);
                // save value in the destination vector
                rDestination[i_node] = scalar;
            }

            KRATOS_CATCH("");
        }

        void WriteVectorToDatabase(Variable<array_1d<double, 3>> &rVariable, CalcVectorType &rOrigin, ModelPart::NodesContainerType &rNodes)
        {
            KRATOS_TRY

            // loop over all nodes
            int n_nodes = rNodes.size();
            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                unsigned int i_node = i;

                // get reference of destination
                array_1d<double, 3> &vector = node_it->FastGetCurrentSolutionStepValue(rVariable, var_pos);
                // save vector in database
                for (unsigned int component = 0; component < TDim; component++)
                    vector[component] = (rOrigin[i_node])[component];
            }

            KRATOS_CATCH("");
        }

        void WriteScalarToDatabase(Variable<double> &rVariable, ValuesVectorType &rOrigin, ModelPart::NodesContainerType &rNodes)
        {
            KRATOS_TRY

            // loop over all nodes
            int n_nodes = rNodes.size();
            ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();

            unsigned int var_pos = it_begin->pGetVariablesList()->Index(rVariable);

#pragma omp parallel for firstprivate(n_nodes, it_begin, var_pos)
            for (int i = 0; i < n_nodes; i++)
            {
                ModelPart::NodesContainerType::iterator node_it = it_begin + i;

                // get the global index of node i
                int i_node = i;

                // get reference of destination
                double &scalar = node_it->FastGetCurrentSolutionStepValue(rVariable, var_pos);
                // save scalar in database
                scalar = rOrigin[i_node];
            }

            KRATOS_CATCH("");
        }

        //*********************************************************************
        // destination = origin1 + value * Minv*origin

        void Add_Minv_value(
            CalcVectorType &destination,
            const CalcVectorType &origin1,
            const double value,
            const ValuesVectorType &Minv_vec,
            const CalcVectorType &origin)
        {
            KRATOS_TRY

            int loop_size = destination.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &dest = destination[i_node];
                const double m_inv = Minv_vec[i_node];
                const array_1d<double, TDim> &origin_vec1 = origin1[i_node];
                const array_1d<double, TDim> &origin_value = origin[i_node];

                double temp = value * m_inv;
                for (unsigned int comp = 0; comp < TDim; comp++)
                    dest[comp] = origin_vec1[comp] + temp * origin_value[comp];
            });

            KRATOS_CATCH("")
        }

        void Add_Minv_value(
            ValuesVectorType &destination,
            const ValuesVectorType &origin1,
            const double value,
            const ValuesVectorType &Minv_vec,
            const ValuesVectorType &origin)
        {
            KRATOS_TRY

            int loop_size = destination.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                double &dest = destination[i_node];
                const double m_inv = Minv_vec[i_node];
                const double &origin_vec1 = origin1[i_node];
                const double &origin_value = origin[i_node];

                double temp = value * m_inv;
                dest = origin_vec1 + temp * origin_value;
            });

            KRATOS_CATCH("")
        }

        //**********************************************************************

        void AllocateAndSetToZero(CalcVectorType &data_vector, int size)
        {
            data_vector.resize(size);
            int loop_size = size;

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &aaa = data_vector[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    aaa[comp] = 0.0;
            });
        }

        void AllocateAndSetToZero(ValuesVectorType &data_vector, int size)
        {
            data_vector.resize(size);
            int loop_size = size;

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                data_vector[i_node] = 0.0;
            });
        }

        //**********************************************************************

        void SetToZero(EdgesVectorType &data_vector)
        {
            int loop_size = data_vector.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                // initializing the CSR vector entries with zero
                data_vector[i_node].Mass = 0.0;
                noalias(data_vector[i_node].LaplacianIJ) = ZeroMatrix(TDim, TDim);
                noalias(data_vector[i_node].Ni_DNj) = ZeroVector(TDim);
                noalias(data_vector[i_node].DNi_Nj) = ZeroVector(TDim);
            });
        }

        void SetToZero(IndicesVectorType &data_vector)
        {
            int loop_size = data_vector.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                data_vector[i_node] = 0.0;
            });
        }

        void SetToZero(CalcVectorType &data_vector)
        {
            int loop_size = data_vector.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &aaa = data_vector[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    aaa[comp] = 0.0;
            });
        }

        void SetToZero(ValuesVectorType &data_vector)
        {
            int loop_size = data_vector.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                data_vector[i_node] = 0.0;
            });
        }

        //**********************************************************************

        void AssignVectorToVector(const CalcVectorType &origin,
                                  CalcVectorType &destination)
        {
            int loop_size = origin.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                const array_1d<double, TDim> &orig = origin[i_node];
                array_1d<double, TDim> &dest = destination[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    dest[comp] = orig[comp];
            });
        }

        void AssignVectorToVector(const ValuesVectorType &origin,
                                  ValuesVectorType &destination)
        {
            int loop_size = origin.size();

            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                destination[i_node] = origin[i_node];
            });
        }

    private:
        // number of edges
        unsigned int mNumberEdges;

        // CSR data vector for storage of the G, L and consistent M components of edge ij
        EdgesVectorType mNonzeroEdgeValues;

        // vector to store column indices of nonzero matrix elements for each row
        IndicesVectorType mColumnIndex;

        // index vector to access the start of matrix row i in the column vector
        IndicesVectorType mRowStartIndex;

        // inverse of the mass matrix ... for parallel calculation each subdomain should contain this correctly calculated (including contributions of the neighbours)
        ValuesVectorType mInvertedMassMatrix;

        // minimum height around one node
        ValuesVectorType mHmin;

        // lumped mass matrix (separately stored due to lack of diagonal elements of the consistent mass matrix)
        ValuesVectorType mLumpedMassMatrix;
        // diagonal of the gradient matrix (separately stored due to special calculations)
        CalcVectorType mDiagGradientMatrix;

        //*******************************************
        // functions to set up elemental mass matrices

        void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3> &mass_consistent, double volume)
        {
            for (unsigned int i_node = 0; i_node <= TDim; i_node++)
            {
                // diagonal terms
                mass_consistent(i_node, i_node) = 0.16666666666666666667 * volume; // 1/6
                // non-diagonal terms
                double temp = 0.08333333333333333333 * volume; // 1/12
                for (unsigned int j_neighbour = i_node + 1; j_neighbour <= TDim; j_neighbour++)
                {
                    // taking advantage of symmetry
                    mass_consistent(i_node, j_neighbour) = temp;
                    mass_consistent(j_neighbour, i_node) = temp;
                }
            }
        }

        void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix<double, 4, 4> &mass_consistent, double volume)
        {
            for (unsigned int i_node = 0; i_node <= TDim; i_node++)
            {
                // diagonal terms
                mass_consistent(i_node, i_node) = 0.1 * volume;
                // non-diagonal terms
                double temp = 0.05 * volume;
                for (unsigned int j_neighbour = i_node + 1; j_neighbour <= TDim; j_neighbour++)
                {
                    // taking advantage of symmetry
                    mass_consistent(i_node, j_neighbour) = temp;
                    mass_consistent(j_neighbour, i_node) = temp;
                }
            }
        }
    };

} // namespace Kratos

#endif // KRATOS_EDGE_DATA_H_INCLUDED defined
