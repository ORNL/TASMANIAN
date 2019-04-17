/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP
#define __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP

// system headers used in many/many places
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "TasmanianConfig.hpp" // contains build options passed down from CMake
#include "tsgUtils.hpp" // contains array wrapper and size_mult for int-to-size_t

//! \file tsgEnumerates.hpp
//! \brief Omnipresent enumerate types.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianEnumerates
//!
//! Contains the enumerate types that are used throughout Tasmanian.
//! This header is included directly or transitively in every other
//! file of the project including \b TasmanianSparseGrids.hpp

/*!
 * \ingroup TasmanianSG
 * \addtogroup SGEnumerates Enumerated types
 *
 * \par Enumerated types
 * Enumerations are used as input to many Tasmanian functions
 * in both internal and external API.
 *
 * \internal
 * \par Hardcoded constants
 * Hardcoding constants is a bad coding practice and should be avoided,
 * but numerical algorithms depend on many small tweaking parameters that require
 * meaningful default values.
 * For example, Tasmanian computes the nodes and abscissas of the Gauss-Legendre rule on the fly,
 * which is done with an iterative eigenvalue decomposition method that needs a stopping criteria.
 * Extra input parameters can be added to \b makeGlobalGrid() which will allow the user to
 * choose the numeric tolerance and maximum number of iterations, but this will also
 * add extra variables that the user has to understand and specify.
 * The goal of Tasmanian is to provide a seamless experience where the user
 * focuses on the desired properties of a quadrature/interpolation rule, and not
 * worry about convergence of a hidden iterative schemes.
 * Therefore, we hardcoded such variables with \a reasonable value,
 * i.e., values that will work well for most use cases.
 * In fringe cases, the values can be adjusted here and Tasmanian can be recompiled
 * with new constants.
 * \endinternal
 */

namespace TasGrid{

//! \internal
//! \brief Describes the relation between two multi-indexes when compared during sorting.
//! \ingroup SGEnumerates

//! The standard C++11 algorithms std::merge and std::binary_search use bools as return
//! types when comparing entries, thus in the case of many repeated entries, two
//! comparisons have to be performed. Repeated entries are encountered in many sparse
//! grids algorithms, especially when dealing with nested one dimensional rules.
//! Thus, for performance reasons, Tasmanian implements merging and searching
//! methods that compare multi-indexes and produce three outcomes in a single comparison.
//!
//! Specifically, see \b SetManipulations::push_merge_map() and \b MultiIndexSet::getSlot()
//! and \b MultiIndexSet::addSortedInsexes()
enum TypeIndexRelation{ // internal for IndexSets, unambiguous set comparison
    //! \brief Indicates that the first multi-index in the call to compare comes before the second.
    type_abeforeb,
    //! \brief Indicates that the second multi-index in the call to compare after the second.
    type_bbeforea,
    //! \brief Indicates that the two multi-indexes are the same.
    type_asameb
};

//! \brief Used by Global Sequence and Fourier grids, indicates the selection criteria.
//! \ingroup SGEnumerates

//! \par Approximation Error
//! The approximation error when interpolating or integrating a target model with
//! a \b Global sparse grid,
//! is bound by a constant multiplied by the best approximation of the model in the
//! set of polynomial exactness. The exactness set consists of all polynomials that
//! can be \a captured exactly (approximated to round-off error) by the sparse grid.
//! The constant is the
//! norm of the sparse grid operator (a.k.a., the Lebesgue constant)
//! and the exactness set is the range.
//!
//! \par
//! The operator norm depends primarily on the underlying one dimensional rule used
//! to construct the grid, and all rules implemented in Tasmanian have reasonable norms.
//! Of primary focus is the target model and constructing a grid with range suitable
//! for either integrating or interpolating the model.
//!
//! \par Categories
//! The types of \b Global grids fall into four categories:
//! - \b Level selection uses the order of the one dimensional rules as criteria
//!   and pays no regard to the polynomial exactness. \b Level types are strongly influenced
//!   by the growth of points in the one dimensional rule and are seldom optimal,
//!   although such grid can be a good initial guess for an iterative refinement process.
//! - \b Total degree polynomial space is good for smooth models, i.e., analytic functions
//!   with convergent Taylors series for every point in the domain.
//! - \b Curved polynomial space is essentially a \b Total degree space with an
//!   added logarithmic correction term. The correction can account for rules with slightly
//!   higher norm growth (e.g., R-Leja rules) or for models with heavily constrained
//!   region of analyticity, i.e., slightly less smooth.
//! - \b Hyperbolic cross section space is suitable for models with more irregular behavior, e.g.,
//!   models with only a finite number of derivatives.
//!
//! \par
//! The four main types can be tuned for either interpolation or quadrature (integration);
//! the corresponding enumerate will have either \b i or \b q in the name.
//! In all cases, Tasmanian will construct the grid with minimum number of points
//! (or fewest basis functions), that will \a capture the corresponding polynomial space.
//!
//! \par Anisotropic Weights
//! All types can be combined with anisotropic weights that adjust the density
//! of the grid in different directions. Below, we denote the anisotropic weights
//! with \f$ \xi_1, \xi_2, \cdots, \xi_d \f$ and \f$ \eta_1, \eta_2, \cdots, \eta_d \f$,
//! where \f$ d \f$ is the number of dimensions in the grid. Note that only curved
//! types use the \f$ \eta \f$ weights.
//!
//! \par
//! The \b TypeDepth also uses the \b depth parameter to define a cut-off point, the
//! \b depth is normalized by multiplying it by the smallest of the \f$ \xi \f$ weights.
//! In the formulas below, \b L stands for the normalized offset.
//!
//! Details regarding sparse grids construction and polynomial spaces can be found in\n
//! M. Stoyanov, C. G. Webster, <a href="https://arxiv.org/pdf/1508.01125.pdf">
//! A Dynamically Adaptive Sparse Grid Method for Quasi-Optimal
//! Interpolation of Multidimensional Analytic Functions</a>,\n also published in
//! <a href="https://www.sciencedirect.com/science/article/pii/S0898122116000031">
//! Computers and Mathematics with Applications 71 (2016) 2449–2465.</a>
//!
//! \par Full Tensor Types
//! Finally, in addition to the aforementioned types, Tasmanian includes the types for
//! dense tensor grids. Only in 2-dimensions Gauss quadrature rules with full tensors
//! are more efficient than a general sparse grid. Performing interpolation or
//! working with more than 2 dimensions, sparse grids will always result in better
//! approximation error per number of nodes, i.e., model evaluations; nevertheless,
//! dense tensor grids naturally fall in the framework and can be useful for testing
//! purposes.
//!
//! \par Fourier Grids
//! The Fourier grids are very similar to Global grids, with the exception that
//! polynomial powers are replaced by frequencies of trigonometric functions.
//! The different selection strategies can be used with frequency exactness in
//! place of polynomial power.
enum TypeDepth{
    //! \brief Null type, should \b never be used for input, indicates an error of some sort.
    type_none,

    //! \brief Ignoring the polynomial space, use rules with index
    //! \f$ \left\{ (i_1, i_2, \cdots, i_d) : \sum_{k=1}^d \xi_k i_k \leq L \right\} \f$
    type_level,

    //! \brief Ignoring the polynomial space, use rules with index
    //! \f$ \left\{ (i_1, i_2, \cdots, i_d) : \sum_{k=1}^d \xi_k i_k + \sum_{k=1}^d \eta_k \log(i_k + 1) \leq L \right\} \f$
    type_curved,

    //! \brief Ignoring the polynomial space, use rules with index
    //! \f$ \left\{ (i_1, i_2, \cdots, i_d) : \prod_{k=1}^d (i_k + 1)^{\xi_k} \leq L \right\} \f$
    type_hyperbolic,

    //! \brief Total degree polynomial space for interpolation, i.e., the span of
    //! \f$ \left\{ \prod_{k=1}^d x^{i_k} : \sum_{k=1}^d \xi_k i_k \leq L \right\} \f$
    type_iptotal,

    //! \brief Total degree polynomial space for quadrature/integration.
    type_qptotal,

    //! \brief Curved polynomial space for interpolation, i.e., the span of
    //! \f$ \left\{ \prod_{k=1}^d x^{i_k} : \sum_{k=1}^d \xi_k i_k +  \sum_{k=1}^d \eta_k \log(i_k + 1) \leq L \right\} \f$
    type_ipcurved,

    //! \brief Curved polynomial space for quadrature/integration.
    type_qpcurved,

    //! \brief Hyperbolic cross section polynomial space for interpolation, i.e., the span of
    //! \f$ \left\{ \prod_{k=1}^d x^{i_k} : \prod_{k=1}^d (i_k + 1)^{\xi_k} \leq L \right\} \f$
    type_iphyperbolic,

    //! \brief Hyperbolic cross section polynomial space for quadrature/integration.
    type_qphyperbolic,

    //! \brief Make a dense tensor grid with rules indexed by \f$ \xi_1, \xi_2, \cdots, \xi_d \f$
    type_tensor,

    //! \brief Make a dense tensor grid with interpolation range that includes
    //! \f$ \left\{ \prod_{k=1}^d x^{i_k} : i_k \leq \xi_k \mbox{ for all } k \right\} \f$
    type_iptensor,

    //! \brief Make a dense tensor grid with quadrature/integration range that includes
    //! \f$ \left\{ \prod_{k=1}^d x^{i_k} : i_k \leq \xi_k \mbox{ for all } k \right\} \f$
    type_qptensor
};

//! \brief Used to specify the one dimensional family of rules that induces the sparse grid.
//! \ingroup SGEnumerates

//! \par One Dimensional Rules
//! A sparse grid is a superposition of tensors of one dimensional rules with varying precision
//! in the different directions. In this context, a \b TypeOneDRule defines a family of one
//! dimensional rules with corresponding abscissas (nodes), quadrature weights, and
//! basis functions.
//!
//! \par Categories
//! The rules fall into 5 categories, corresponding to the 5 main classes of grids.
//! - \b Global using Lagrange polynomials for interpolation and (sometimes) specialized
//!   quadrature weights.
//! - \b Sequence rules are a subset of the \b Global rules, \b Sequence rules are always
//!   nested and each level adds one node to the previous level.
//! - \b Local \b Polynomial rules use uniformly spaced nodes and polynomial basis with
//!   localized support.
//! - \b Wavelet and \b Fourier rules also use uniformly spaced nodes, but with
//!   wavelet and trigonometric basis functions.
//!
//! \par Special Quadrature
//! All rules can be used to approximate integrals of the form \f$ \int_\Gamma f(x) \rho(x) dx \f$
//! where \f$ \Gamma \f$ is the (transformed) domain and \f$ f(x) \f$ is the model.
//! The majority of rules assume that the weight is constant 1, but in some cases
//! more exotic weights are used. Assume that constant weight 1 is used, unless otherwise
//! specified.
//!
//! \par Domain Transformation
//! Most rules are defined on a canonical interval [-1, 1], which can be transformed to an
//! arbitrary [a, b] with \b TasmanianSparseGrid::setDomainTransform(), so long as \b a is less than \b b.
//! Fourier grids use canonical [0, 1] and can be similary transformed to [a, b]. Gauss-Laguerre rule
//! is set to \f$ [0, \infty) \f$ transformed to \f$ [a, \infty) \f$ with scale parameter \b b (see below).
//! Gauss-Hermite is set to \f$ (-\infty, \infty) \f$ with scale and shift parameters.
enum TypeOneDRule{
    //! \brief Null rule, should \b never be used as input (default rule for an empty grid).
    rule_none,
    //! \brief Classic nested rule using Chebyshev nodes with very low Lebesgue constant.
    rule_clenshawcurtis,
    //! \brief Same as \b rule_clenshawcurtis but with modified basis that assumes the model is zero at the boundary.
    rule_clenshawcurtis0,
    //! \brief Similar to \b rule_clenshawcurtis but with nodes strictly in the interior.
    rule_fejer2,
    //! \brief Using Chebyshev nodes with very low Lebesgue constant and slow node growth, but non-nested.
    rule_chebyshev,
    //! \brief Same as \b rule_chebyshev but using only odd levels, partially mitigates the non-nested issues.
    rule_chebyshevodd,
    //! \brief Classic \b sequence rule, moderate Lebesgue constant growth (empirical result only).
    rule_leja,
    //! \brief Same as \b rule_leja but using only odd levels, quadrature is more stable.
    rule_lejaodd,
    //! \brief Classic \b sequence rule based on complex analysis, moderate Lebesgue constant growth (theoretically proven).
    rule_rleja,
    //! \brief Using \b rule_rleja nodes but doubling the nodes every 2 levels, reduces the Lebesgue constant.
    rule_rlejadouble2,
    //! \brief Using \b rule_rleja nodes but doubling the nodes every 4 levels, reduces the Lebesgue constant.
    rule_rlejadouble4,
    //! \brief Same as \b rule_rleja but using only odd levels, quadrature is more stable.
    rule_rlejaodd,
    //! \brief Similar \b sequence to \b rule_rleja but with nodes strictly in the interior.
    rule_rlejashifted,
    //! \brief Same as \b rule_rlejashifted but using only even levels, quadrature is more stable.
    rule_rlejashiftedeven,
    //! \brief Same as \b rule_rlejashifted but doubling the number of nodes per level, which reduced the Lebesgue constant.
    rule_rlejashifteddouble,
    //! \brief A greedy \b sequence rule with nodes placed at the maximum of the Lebesgue function.
    rule_maxlebesgue,
    //! \brief Same as \b rule_maxlebesgue but using only odd levels, quadrature is more stable.
    rule_maxlebesgueodd,
    //! \brief A greedy \b sequence rule with nodes added to minimize the Lebesgue constant.
    rule_minlebesgue,
    //! \brief Same as \b rule_minlebesgue but using only odd levels, quadrature is more stable.
    rule_minlebesgueodd,
    //! \brief A greedy \b sequence rule with nodes added to minimize the norm of the surplus operator.
    rule_mindelta,
    //! \brief Same as \b rule_mindelta but using only odd levels, quadrature is more stable.
    rule_mindeltaodd,
    //! \brief Non-nested rule but optimized for integration.
    rule_gausslegendre,
    //! \brief Same as \b rule_gausslegendre but using only odd levels, partially mitigates the non-nested issues.
    rule_gausslegendreodd,
    //! \brief Nested rule that is optimized for integration, probably the best integration rule in more than 2 dimensions.
    rule_gausspatterson,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{a}^{b} f(x) / \sqrt{(b - x)(x - a)} dx \f$.
    rule_gausschebyshev1,
    //! \brief Same as \b rule_gausschebyshev1 but using only odd levels, partially mitigates the non-nested issues.
    rule_gausschebyshev1odd,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{a}^{b} f(x) \sqrt{(b - x)(x - a)} dx \f$.
    rule_gausschebyshev2,
    //! \brief Same as \b rule_gausschebyshev2 but using only odd levels, partially mitigates the non-nested issues.
    rule_gausschebyshev2odd,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{a}^{b} f(x) (b - x)^\alpha (x - a)^\alpha dx \f$.
    rule_gaussgegenbauer,
    //! \brief Same as \b rule_gaussgegenbauer but using only odd levels, partially mitigates the non-nested issues.
    rule_gaussgegenbauerodd,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{a}^{b} f(x) (b - x)^\alpha (x - a)^\beta dx \f$.
    rule_gaussjacobi,
    //! \brief Same as \b rule_gaussjacobi but using only odd levels, partially mitigates the non-nested issues.
    rule_gaussjacobiodd,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{a}^{\infty} f(x) (x - a)^\alpha \exp \left(-b (x - a) \right) dx \f$.
    rule_gausslaguerre,
    //! \brief Same as \b rule_gausslaguerre but using only odd levels, partially mitigates the non-nested issues.
    rule_gausslaguerreodd,
    //! \brief Non-nested rule optimized for integral of the form \f$ \int_{-\infty}^{\infty} f(x) (x - a)^\alpha \exp \left(-b (x - a)^2 \right) dx \f$.
    rule_gausshermite,
    //! \brief Same as \b rule_gausshermite but using only odd levels, partially mitigates the non-nested issues.
    rule_gausshermiteodd,
    //! \brief User provided rule, nodes and weights must be provided with a separate file.
    rule_customtabulated,
    // Piece-Wise rules
    //! \brief Nested rule with a hierarchy of uniformly distributed nodes and functions with compact support.
    rule_localp,
    //! \brief Variation of \b rule_localp assuming the model is zero at the domain boundary.
    rule_localp0,
    //! \brief Variation of \b rule_localp using increased support in exchange for higher order basis (better for smoother functions).
    rule_semilocalp,
    //! \brief Variation of \b rule_localp focusing nodes on the boundary instead of the interior.
    rule_localpb,
    //! \brief Wavelet basis with uniformly distributed nodes (primarily for internal use).
    rule_wavelet,
    //! \brief Trigonometric basis with uniformly distributed nodes (primarily for internal use).
    rule_fourier
};

//! \brief Refinement strategy for local polynomial and wavelet grids.
//! \ingroup SGEnumerates

//! \par Local Hierarchical Approximation
//! The nodes and basis functions used in local polynomial and wavelet sparse grid form a complex multidimensional hierarchy,
//! where each node is associated with one or more parents and children in each direction. In both cases, the support of the
//! children is (almost always) smaller than the support of the patent and in the polynomial case, the support of the
//! children nodes is strictly contained within the support of the parent.
//!
//! \par
//! The interpolant is expressed as a linear combination of basis functions with coefficients computed so that the
//! interpolant matches the model values at the sparse gird nodes. The coefficients are thus indicators of the local
//! approximation error and (in the asymptotic limit) the coefficients decrease when going form parent to children.
//!
//! \par Convergence Criteria
//! Using the local error estimator, the sparse grid approximation has reached desired tolerance
//! when all nodes with missing children
//! have coefficients with absolute values smaller than the tolerance. This is a necessary and not a sufficient condition
//! an it  heavily relies on the assumption of monotonic decay of the coefficients; however, in practical situations this is a good
//! way to construct an adaptive grid that captures the behavior of the target model while minimizing the number of nodes.
//!
//! \par Adaptive Refinement
//! The idea of refinement is to start with a coarse spars grid and gradually add children of existing nodes, but only if the nodes
//! have a large hierarchical coefficient. The nodes are added until the coefficients of all nodes with with missing children
//! drop below user specified tolerance. The refinement types listed below differ in the way that the parents and children
//! are selected in the adaptive construction, and whether local anisotropy is considered.
//!
//! \par Refinement Types
//! Each nodes is associated with multiple parents corresponding
//! to different directions, even if each one dimensional node has only a single parent. Thus, over several iterations
//! and if the hierarchical coefficients do not decay monotonically, it is possible to have nodes with missing parents
//! for some directions, which can hinter convergence. Furthermore, the hierarchical coefficients do not carry any
//! information regarding potential local anisotropy and a second refinement criteria can be employed, where an interpolant
//! is constructed using only the data associated with the nodes that lay in a single direction. Children and/or parents of
//! a node are added to the grid only when both the isotropic and anisotropic error indicators exceed the user specified tolerance.
//!
//! \par
//! The \b classic refinement strategy adds only the children for the nodes with hierarchical coefficient that exceed the
//! user provided tolerance (local anisotropy and missing parents are ignored). The \b parents \b first and \b fds
//! strategies loop for missing parents and add the parents before adding children, which usually improves stability
//! when compared to adding only children. The \b direction \b selective and \b fds strategies consider local anisotropy
//! which can reduce the number of required model evaluations, but can also lead to decrease in stability.
//!
//! Details regarding adaptive hierarchical sparse grids construction and children-parent relations can be found in: \n
//! M. Stoyanov, <a href="https://link.springer.com/chapter/10.1007/978-3-319-75426-0_8">
//! Adaptive Sparse Grid Construction in a Context of Local Anisotropy and Multiple Hierarchical Parents</a>,\n
//! Sparse Grids and Applications - Miami 2016 pp 175-199.\n
//! Also in: <a href="https://info.ornl.gov/sites/publications/files/Pub45746.pdf">Tech Report ORNL/TM-2013/384.</a>
enum TypeRefinement{
    //! \brief Isotropic refinement using only the children and disregarding missing parents.
    refine_classic,
    //! \brief Isotropic refinement adding children only if the parents are already included.
    refine_parents_first,
    //! \brief Anisotropic refinement using only the children and disregarding missing parents.
    refine_direction_selective,
    //! \brief Anisotropic refinement adding children only if the parents are already included.
    refine_fds,
    //! \brief Null method, should \b never be used as input.
    refine_none
};

//! \brief Types of acceleration for \b TasmanianSparseGrid::evaluateFast() and \b TasmanianSparseGrid::evaluateBatch().
//! \ingroup SGEnumerates

//! \par Evaluation Stages
//! Computing the value of the interpolant at one or more nodes is done in two stages, first a matrix of basis function
//! values is constructed, then the matrix is multiplied by a persistent data matrix (containing model values, or
//! the corresponding hierarchical or Fourier coefficients). The main factors that determine the performance are
//! the number of basis functions (i.e., grid nodes), the number of model outputs, and the size of the batch.
//!
//! \par
//! The BLAS and cuBLAS modes use the CPU (and OpenMP) for the first stage of evaluations, and then the corresponding
//! library for fast linear algebra. Both modes are most efficient when dealing with many outputs.
//! The CUDA mode uses custom kernels for basis evaluation and can be very fast especially when the grid has many nodes,
//! the second stage of the evaluations uses the cuBLAS library. The MAGMA mode is the same as CUDA but using the
//! <a href="http://icl.cs.utk.edu/magma/">UTK MAGMA library.</a> The \b none mode uses only OpenMP (if available), but
//! the implementation minimizes the memory footprint, the matrices are not formed explicitly and thus the mode
//! could be the fastest when working with small grids.
//!
//! \par Defaults
//! If \b Tasmanian_ENABLE_BLAS is enabled in CMake, then \b accel_cpu_blas will be the default acceleration mode,
//! otherwise the default mode is \b accel_none.
//! The \b accel_gpu_default uses the first available in the list: \b accel_gpu_magma, \b accel_gpu_cuda,
//! \b accel_cpu_blas, \b accel_none.
//! When using \b accel_none, \b TasmanianSparseGrid::evaluateFast() is equivalent to \b TasmanianSparseGrid::evaluate()
//! and equivalent to \b TasmanianSparseGrid::evaluateBatch() with \b num_x equal to 1.
//!
//! \par
//! Note that regardless of the acceleration mode, all evaluate methods are required to generate output identical up to
//! rounding error.
//! Thus, the Tasmanian API is designed so that \b TasmanianSparseGrid::enableAcceleration() can be called
//! with any acceleration type regardless whether it has been enabled by CMake.
//! In order to facilitate code portability, Tasmanian implements sane fall-back
//! modes so that "the next best mode" will be automatically selected if the desired mode is not present. Specifically,
//! missing MAGMA will fallback to CUDA, missing CUDA will fallback to BLAS, and missing BLAS will fallback to \b none.
//!
//! \par Memory Management
//! All acceleration modes (other than \b none) require that potentially large matrices are explicitly constructed,
//! which leads to an increase of memory footprint. Furthermore, the GPU modes will copy the matrices to the
//! device memory, which may not be feasible when using weaker GPU devices or working with very large problems.
//!
//! \par
//! Global, Sequence and Fourier grids operate only with dense matrices and require one matrix of size \b getNumPoints()
//! times \b getNumOutputs(), one matrix of size \b getNumOutputs() times \b num_x, and one matrix of size \b getNumPoints()
//! times \b num_x. For \b evaluateFast(), \b num_x is equal to 1. The batch size \b num_x can be adjusted to accommodate
//! memory constraints, but the first stage of evaluations achieved better performance when \b num_x is large (especially
//! when executed on the GPU).
//!
//! \par
//! Local polynomial grids can perform evaluations in either dense or sparse modes. Dense mode has memory requirements
//! identical to the Global case, but dense evaluations are less efficient when working with fewer grid dimensions
//! and few grid nodes.
//! The sparse mode requires two matrices of size \b getNumOutputs() times \b num_x, one matrix of size \b getNumPoints()
//! times \b getNumOutputs(), and one sparse matrix of size \b getNumPoints() times \b num_x. The actual sparsity
//! pattern and the exact memory requirements are hard to predict, but can also be controlled through the batch size,
//! albeit with larger performance penalty for small batches. By default, Tasmanian will automatically select sparse
//! or dense mode based on criteria derived from numerous in-house benchmark tests,
//! but the mode can also be manually selected with \b TasmanianSparseGrid::favorSparseAcceleration().
//!
//! \par
//! In all cases, a matrix of size \b getNumPoints() times \b getNumOutputs() will be stored in device memory and will
//! be kept persistent until the gird is modified (points are added/removed or the object is removed). The persistent
//! memory can be cleared manually by switching to non-GPU acceleration mode.
//!
//! \par Thread Safety
//! The \b accel_none mode and the simple \b TasmanianSparseGrid::evaluate() function are thread safe, i.e., can be called
//! concurrently from multiple threads and together with any other const methods and will always return result that vary
//! by no more than rounding error. Using \b accel_cpu_blas, the \b evaluateBatch() and \b evaluateFast() are thread save
//! to the extent to which the BLAS library is thread save; for example, OpenBLAS must be build with a special flag
//! otherwise concurrent calls to BLAS functions can result in incorrect output.
//!
//! \par
//! The GPU methods use mutable internal data-structures for cuBLAS and MAGAM handles and to load the persistent
//! matrix into GPU memory. Thus, calls to \b evaluateBatch() and \b evaluateFast() are not thread safe when using
//! any of the GPU acceleration modes. However, the mutable data-structures are encapsulated within the corresponding
//! object, thus evaluations for different objects can be performed concurrently and different objects can be associated
//! with different GPU devices without causing any conflicts. Note that this is a feature of Tasmanian and some third
//! party libraries (e.g., OpenBLAS and \b accel_cpu_blas mode) can still cause issues even when different threads
//! manipulate different objects.
enum TypeAcceleration{
    //! \brief Usually the slowest mode, uses only OpenMP multi-threading, but optimized for memory and could be the fastest mode for small problems.
    accel_none,
    //! \brief Default (if available), uses OpenMP to form matrices and BLAS on the CPU for linear algebra.
    accel_cpu_blas,
    //! \brief Equivalent to the first available from \b accel_gpu_magma, \b accel_gpu_cuda, \b accel_cpu_blas, \b accel_none.
    accel_gpu_default,
    //! \brief Using the CPU to form the matrices and the GPU for linear algebra, good when the model number of outputs is large.
    accel_gpu_cublas,
    //! \brief Same as \b accel_gpu_cublas but uses the GPU to form the matrices and is thus faster when using many grid nodes and large batch size.
    accel_gpu_cuda,
    //! \brief Same as \b accel_gpu_cuda but uses the UTK MAGAM library for the linear algebra.
    accel_gpu_magma
};


////////////////////////////////////////////////////////////
//                                                        //
//     ---  Moved from tsgHardcodedConstants.hpp  ---     //
//                                                        //
////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  On the purpose of this section:
//
//
//  Hardcoding constants is a bad coding practice and should be avoided.
//  On the other hand, numerical algorithms depend on many small tweaking parameters.
//  For example, this code computes the nodes and abscissas of Gauss-Legendre rule on the fly,
//  this is done with an iterative eigenvalue decomposition method that needs a stopping criteria,
//  we can add another user specified parameter to the corresponding "makeGlobalGrid()" rule,
//  however, this is an additional tweaking variable that the user has to understand and control.
//  The goal of our code is to provide a most seamless experience to the user,
//  one should think about the desired properties of a quadrature rule or interpolant and not
//  about convergence of an iterative scheme (especially one hidden from top view)
//  Therefore, the tolerance for such convergence criteria needs to be set inside the code to
//  a "reasonable" value, which is a value that would work well for the overwhelming majority
//  of use cases.
//
//  On the other hand, every application is different and the usage of the code will change
//  over time. It is unreasonable to believe that one single value would work for absolutely
//  everyone. Instead of "hiding" hardcoded constant throughout the code, all such constants
//  will be exposed here so the user can make adjustments in compile time.
//
//  Long story short: do not adjust those variables unless you have a good reason.
//
///////////////////////////////////////////////////////////////////////////////////////////////

//! \internal
//! \brief Numerical tolerance for various algorithms.
//! \ingroup SGEnumerates

//! NUM_TOL is used in many places:
//! - as a stopping criteria for various iterative schemes (e.g., finding leja points)
//! - drop criteria for eigenvalue solver related to Gauss rules
//! - comparison between nodes to detect repeated points in non-nested rules (e.g., all odd Chebyshev rules include zero)
//! - determining sparse matrix pattern, entries smaller than NUM_TOL will be ignored (for wavelet grids)
//! - drop criteria in estimating anisotropic coefficients (refinement or just the coefficients) surpluses or Legendre coefficients below 10^3 times NUM_TOL will be ignored
constexpr double TSG_NUM_TOL = 1.E-12;

//! \internal
//! \brief Defines the maximum number of secant method iterations to be used for finding Leja, Lebesgue, and Delta points.
//! \ingroup SGEnumerates

//! This is a simple safeguard criteria to prevent "hanging" in a loop.
constexpr int TSG_MAX_SECANT_ITERATIONS = 1000;

//! \internal
//! \brief Tuning parameter for dense vs sparse evaluations of Local Polynomial Grids when using CPU BLAS.
//! \ingroup SGEnumerates

//! Defines the threshold for switching between sparse and dense version of batch evaluate for Local Polynomial Grids.
//! Generally, because of caching and reuse of data, dense operations have 10x more flops per second than sparse ops;
//! unless the sparse version of the algorithm can cut total flops by 10x (i.e., fill is less than 10%), then it
//! is faster to convert the sparse matrix into a dense one and use dense linear algebra.
//! The difference is noticeable if the number of outputs is sufficiently large;
//! if the outputs are few there is very little caching anyway and hence use the sparse version to save memory.
constexpr size_t TSG_LOCALP_BLAS_NUM_OUTPUTS = 2048;

/*!
 * \internal
 * \ingroup SGEnumerates
 * \brief Half-period of the \b std::sin() and \b std::cos() functions.
 *
 * Borrowed from "math.h", let's hope they don't change the value in a future standard.
 */
constexpr double tsg_pi = 3.14159265358979323846;

}

#endif
