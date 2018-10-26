//! \mainpage Introduction
//! \section Intro Tasmanian - Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN
//!
//! The Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN is a collection of robust
//! libraries for high dimensional integration, interpolation, and random sampling.
//! The code consists of several modules that can be used individually or conjointly.
//!
//!
//! \section SparseGrids Sparse Grids
//! Sparse Grids is a family of algorithms for constructing multidimensional quadrature and interpolation rules
//! using multiple tensor products of one dimensional rules with varying degree of precision.
//! The Tasmanian Sparse Grids Module implements a variety of grids that fall into five major categories:
//! * **Global grids** use polynomials supported over the entire domain of integration or interpolation.
//!     Such grids are suitable for approximating functions with smooth response.
//! * **Sequence grids** work much like Global grids, but use optimized internal data-structures for rules that
//!     are based on sequences of points formed from solving a greedy optimization problem
//! * **Local polynomial grids** use hierarchical piece-wise polynomials with compact support.
//!     Such grids are suitable for approximating functions with localized sharp behavior.
//! * **Wavelet grids** use special functions that form a Riesz basis which allows for more accurate
//!     local error estimations. Compared to Local polynomial grids, the wavelet basis can provide
//!     similar accuracy with significantly fewer points, but the advantage of the Riesz basis could also
//!     be negated from the higher Lebesgue constant near the domain boundary.
//! * **Fourier grids** use trigonometric functions with varying frequency and rely on Discrete (complex)
//!     Fourier Transforms. Such grids are suitable for approximating periodic functions,
//!     since periodicity if preserved in the approximation.
//!
//! \section DREAM DREAM
//! The DiffeRential Evolution Adaptive Metropolis is a method to draw samples from an arbitrary probability
//! distribution defined by an arbitrary non-negative function (not necessarily normalized to integrate to 1).
//! The DREAM approeach is similar to the classical Markov Chain Monte Carlo, but it evolves a large number
//! of chains simultanously which leads to better parallelization and (potentially) faster convergence.
//! In addition, multiple chains allow for better exploration of the probablility domain, which is often
//! advantageous when woring with multi-modal distributions.
//!
//! One of the main applications of DREAM is in the field of Bayesian inference, where samples are drawn
//! from a posterior distribution comprised from a data-informed likelihood and an arbitrary model.
//! The DREAM module of Tasmanian can use Tasmanian Sparse Grids approximation to either the model
//! or the likelihood.
