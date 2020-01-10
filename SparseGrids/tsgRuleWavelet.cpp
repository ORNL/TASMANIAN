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

#include "tsgRuleWavelet.hpp"

namespace TasGrid{

inline int ACCESS_FINE(int i, int level, int depth){ return ((1 << (depth - level - 1)) * (2 * i + 1)); }
inline int ACCESS_COARSE(int i, int level, int depth){ return i * (1 << (depth - level)); }

void RuleWavelet::updateOrder(int ord){
    // Changes the order of the rule to the specified order. If order other than 1 is
    // specified, then the approximation to the wavelets will be recalculated.
    if(order == ord) return;

    // clear is practically free, call it every time
    data = std::vector<std::vector<double>>();
    cachexs = std::vector<double>();

    order = ord;

    if(order == 3){

        data.resize(5); // (xs, level1 (scaling), level2, level3, level4)
        data[0].resize(num_data_points);
        double *xs = data[0].data();

        for(int i = 0; i < num_data_points; i++){
            xs[i] = -1. + 2*(double (i)) / (double (num_data_points - 1));
        }

        cachexs = std::vector<double>((size_t) 4 * (num_data_points - 3));
        for(int i = 0; i < num_data_points - 3; i++){
            double const *s = &(xs[i]);
            double *c = &(cachexs[4*i]);
            double alpha = 1.0 / ((s[0] - s[3]) * (s[0] - s[2]));
            double beta  = 1.0 / ((s[0] - s[3]) * (s[1] - s[3]));
            c[0] =   alpha / (s[0] - s[1]);
            c[1] = - alpha / (s[0] - s[1]) - alpha / (s[1] - s[2]) - beta / (s[1] - s[2]);
            c[2] =   alpha / (s[1] - s[2]) + beta  / (s[1] - s[2]) + beta / (s[2] - s[3]);
            c[3] = - beta  / (s[2] - s[3]);
        }

        // Coefficients derived by solving linear system involving scaling function
        // integrals and moments.
        std::vector<std::vector<double>> coeffs(3);
        coeffs[0] = {0.95958116146167449, 0.27778015867946454, -0.042754937296610125, -0.014317145809288835, // Second level
                    -0.34358723961941895, 0.36254192315551625, 0.20438681551383264, 0.012027193241660438};
        coeffs[1] = {0.90985488901447964, 0.29866296454372104, -0.077377811931657145, 0.017034083534297827, // third level
                    -0.28361250628699103, 0.33723841173046543, 0.24560604387620491, -0.023811935316236606,
                     0.015786802304611155, 0.23056985382971185, 0.31221657974498912, -0.03868102549989539,
                     0.194797093133632, -0.099050236091189195, 0.51520570019738199, -0.073403162960780324};
        coeffs[2] = {0.90985443399962629, 0.29866318097339162, -0.077377917373678995, 0.017034105626588848, // fourth level
                    -0.28361254472311259, 0.33723844527609648, 0.24560602528531161, -0.023811931411878345,
                     0.015786801751471, 0.23056986060092344, 0.31221657434994671, -0.038681024059425327,
                     0.14697942814289319, -0.047340431972365773, 0.51871874565793186, -0.093244980946745618,
                    -0.019628614067688826, 0.25611706552019509, 0.30090457777800012, -0.036629694186987166,
                    -1.0/32.0, 9.0/32.0, 9.0/32.0, -1.0/32.0};

        // Initialize scaling functions
        data[1].resize(3*num_data_points);
        std::fill(data[1].begin(), data[1].end(), 0.0);

        // Point (sparse grid numbering):
        // 1     3     0     4     2
        // X --- X --- X --- X --- X
        // 0     1     2     3     4
        // Point (level 2 coarse indexing)

        // This ordering makes phi1 -> point 0, phi2 -> point 1, phi3 -> point 3
        // Points 2 & 4 can be found by reflection of phi2, phi3, respectively.
        data[1][ACCESS_COARSE(2, 2, iteration_depth)] = 1.0;
        data[1][num_data_points + ACCESS_COARSE(0, 2, iteration_depth)] = 1.0;
        data[1][2*num_data_points + ACCESS_COARSE(1, 2, iteration_depth)] = 1.0;
        cubic_cascade(data[1].data(), 2, iteration_depth);
        cubic_cascade(&(data[1][num_data_points]), 2, iteration_depth);
        cubic_cascade(&(data[1][2*num_data_points]), 2, iteration_depth);

        for(int level = 2; level <= 4; level++){
            // Scaling functions at current level.
            std::vector<double> phi1(num_data_points, 0.0);
            std::vector<double> phi2(num_data_points, 0.0);
            std::vector<double> phi3(num_data_points, 0.0);
            std::vector<double> phi4(num_data_points, 0.0);
            int num_saved = 2 * (level-1); // 2 'unique' functions at level 2, 4 at 3, 6 at 4.
            data[level].resize(num_saved * num_data_points);
            std::fill(data[level].begin(), data[level].end(), 0.0);

            // Initialize first four scaling functions
            phi1[ACCESS_COARSE(0, level, iteration_depth)] = 1.0;
            phi2[ACCESS_COARSE(1, level, iteration_depth)] = 1.0;
            phi3[ACCESS_COARSE(2, level, iteration_depth)] = 1.0;
            phi4[ACCESS_COARSE(3, level, iteration_depth)] = 1.0;
            cubic_cascade(phi1.data(), level, iteration_depth);
            cubic_cascade(phi2.data(), level, iteration_depth);
            cubic_cascade(phi3.data(), level, iteration_depth);
            cubic_cascade(phi4.data(), level, iteration_depth);

            for(int index = 0; index < num_saved; index++){
                // Initialize unlifted wavelet
                double *wavelet = &data[level][index*(num_data_points)];
                wavelet[ACCESS_FINE(index, level, iteration_depth)] = 1;
                cubic_cascade(wavelet, level, iteration_depth);
                double  c1 = coeffs[level-2][4*index],
                        c2 = coeffs[level-2][4*index+1],
                        c3 = coeffs[level-2][4*index+2],
                        c4 = coeffs[level-2][4*index+3];

                if(level > 2 && index >= 2){
                    // Change pointers around to avoid copying data
                    std::vector<double> tmp;
                    std::swap(tmp, phi1);
                    std::swap(phi1, phi2);
                    std::swap(phi2, phi3);
                    std::swap(phi3, phi4);
                    std::swap(phi4, tmp);

                    // Initialize new scaling function
                    std::fill(phi4.begin(), phi4.end(), 0.0);
                    phi4[ACCESS_COARSE(index+2, level, iteration_depth)] = 1;
                    cubic_cascade(phi4.data(), level, iteration_depth);

                }

                for(int i = 0; i < num_data_points; i++){
                    // Lift the wavelet
                    wavelet[i] -= c1 * phi1[i] + c2 * phi2[i] + c3 * phi3[i] + c4 * phi4[i];
                }
            }
        }
    }
}

void RuleWavelet::cubic_cascade(double *y, int starting_level, int in_iteration_depth){
    // Using the cascade algorithm (interpolating subdivision), approximates the values
    // of the desired scaling function or wavelet at 2^(iteration_depth)+1 points.
    // Wavelets are generated by placing a 1 in the appropriate place at a fine point
    // on a given level while scaling functions are generated by placing a 1 in the
    // appropriate place at a coarse point.
    for(int level = starting_level; level < in_iteration_depth; level++){
        int num_pts = (1 << (level));
        int prev_pts = (1 << (level)) + 1;

        auto AC = [&](int i, int l)->int{ return ACCESS_COARSE(i, l, in_iteration_depth); };
        auto AF = [&](int i, int l)->int{ return ACCESS_FINE(i, l, in_iteration_depth); };

        // Boundary predictions
        y[AF(0,level)] += (
                5 *(y[AC(0,level)] + 3*y[AC(1,level)] - y[AC(2,level)])
                + y[AC(3,level)]
                ) / 16.;

        y[AF(num_pts-1,level)] += (
                    5 *(y[AC(prev_pts-1,level)] + 3*y[AC(prev_pts-2,level)] - y[AC(prev_pts-3,level)])
                    + y[AC(prev_pts-4,level)]
                    ) / 16.;

        // Central predictions
        #pragma omp parallel for
        for(int i = 1; i < num_pts-1; i++){
            y[AF(i,level)] += (
                9 * (y[AC(i,level)] + y[AC(i+1,level)])
                  - (y[AC(i-1,level)] + y[AC(i+2,level)])
                ) / 16.;
        }
    }
}

int RuleWavelet::getOrder() const{
    return order;
}

int RuleWavelet::getNumPoints(int level) const{
    // Returns the number of points on a given level.
    if(order == 1){
        return (1 << (level + 1)) + 1;
    }else{
        return (1 << (level + 2)) + 1;
    }
}

const char * RuleWavelet::getDescription() const{
    if (order == 1){
        return "First-Order Wavelet Basis";
    }else{
        return "Third-Order Wavelet Basis";
    }
}

int RuleWavelet::getLevel(int point) const{
    // Returns the level to which the given node belongs.
    if(order == 1){
        return (point <= 2) ? 0 : Maths::intlog2(point - 1);
    }else{
        return (point < 5) ? 0 : Maths::intlog2(point - 1) - 1;
    }
}
void RuleWavelet::getChildren(int point, int &first, int &second) const{
    // Returns the children of the given node in first and second. If the node has only a
    // single child, then second is set to -1.
    if (order == 1){
        if (point >= 3){ // most likely case
            first = 2*point-1;
            second = 2*point;
        }else if (point < 2 ){ // second most likely case
            first = 3;
            second = (point == 0) ? 4 : -1;
        }else{ // point == 2
            first = 4;
            second = -1;
        }
    }else if (order == 3){
        if (point >= 3){
            first = 2*point-1;
            second = 2*point;
        }else if (point == 0){
            first = 6;
            second = 7;
        }else if (point == 1){
            first = 5;
            second = -1;
        }else{ // point == 2
            first = 8;
            second = -1;
        }
    }
}
int RuleWavelet::getParent(int point) const{
    // Returns the parent of a given node.
    // Miro: this is a hack, -1 indicates no parent, >= 0 indicates the parent, -2 indicates all nodes on level 0
    if (order == 1){
        if(point <= 2) return -1;
        if(point <= 4) return -2;
    }else{
        if(point <= 4) return -1;
        if(point <= 8) return -2;
    }
    return (point+1)/2;
}

double RuleWavelet::getNode(int point) const {
    // Returns the x-coordinate in the canonical domain associated with the given wavelet.
    if (point == 0) return  0.0;
    if (point == 1) return -1.0;
    if (point == 2) return  1.0;
    return ((double)(2*point - 1)) / ((double) Maths::int2log2(point - 1)) - 3.0;
}

double RuleWavelet::getWeight(int point) const{
    // Returns the integral of the given wavelet.
    if (order == 1){
        return (point == 0) ? 1.0 : (point <= 2) ? 0.5 : 0.0;
    }else if(order == 3){
        if (point > 4){
            return 0.0;
        }
        switch (point){
        case 1:
        case 2:
            return 1.680555555555555691e-01;
        case 3:
        case 4:
            return 6.611111111111110938e-01;
        case 0:
            return 3.416666666666666741e-01;
        }
    }
    return 0.0;
}

double RuleWavelet::eval(int point, double x) const{
    // Evaluates a wavelet designated by point at coordinate x.
    if(order == 1){
        // Level 0
        if (point < 3){
            double w = 1.0 - std::abs(x - getNode(point));
            return (w < 0.0) ? 0.0 : w;
        }
        // Level 1+
        return eval_linear(point, x);
    }
    else if(order == 3){
        return eval_cubic(point, x);
    }
    return 0.;
}

inline double RuleWavelet::eval_cubic(int point, double x) const{
    // Evaluates a third order wavelet at a given point x.
    if (point < 5){ // Scaling functions
        if (point == 2){ // Reflect across y-axis
            point = 1;
            x = -x;
        }else if(point == 4){
            point = 3;
            x = -x;
        }
        const double *phi = &(data[1][((point+1)/2) * num_data_points]);
        return interpolate(phi, x);
    }
    int l = Maths::intlog2(point - 1);

    if(l == 2){
        if (point > 6){
            // i.e. 7 or 8
            // These wavelets are reflections across the y-axis of 6 & 5, respectively
            x = -x;
            point = 13 - point;
        }
        point -= 5;
        return interpolate(&data[2][point*num_data_points],x);
    }else if(l == 3){
        if (point > 12){
            // i.e. 13, 14, 15, 16
            // These wavelets are reflections of 12, 11, 10, 9, respectively
            x = -x;
            point = 25 - point;
        }
        point -= 9;
        return interpolate(&data[3][point*num_data_points],x);
    }
    // Standard lifted wavelets.
    int subindex = (point - 1) % (1 << l);
    double scale = pow(2,l-4);
    // Left Boundary
    if (subindex < 5){
        return interpolate(&data[4][subindex*num_data_points],scale * (x + 1.) - 1.);
    }
    // Right Boundary
    if ((1 << l) - 1 - subindex < 5){
        return interpolate(&data[4][((1 << l) - subindex - 1)*num_data_points],scale * (1. - x) - 1.);
    }
    // Center
    double shift = 0.125 * (double (subindex - 5));
    return interpolate(&data[4][5*num_data_points], scale * (x + 1.) -1. - shift);
}

void RuleWavelet::getShiftScale(int point, double &scale, double &shift) const{
    if (point < 3){
        scale = getNode(point);
        shift = -1.0;
    }else{
        int l = Maths::intlog2(point - 1);
        int subindex = (point - 1) % (1 << l);
        scale = std::pow(2,l-2);
        if (subindex == 0){
            shift = -2.0;
        }else if (subindex == (1 << l) - 1){
            shift = -3.0;
        }else{
            shift = 0.5 * (double (subindex - 1));
        }
    }
}
double RuleWavelet::getSupport(int point) const{
    if (order == 1){
        if (point <= 2) return 1.0; // linear-polynomial points
        return 3.0  / ( 4.0 * std::pow(2.0, Maths::intlog2(point - 1) - 2) );
    }else{
        if (point <= 8) return 2.0; // cubic-polynomial points and level 1 have global support
        // the strange 4.2 compensates for numerical errors in computing the cubic-cascade
        return 4.2 / (3.0 * std::pow(2.0, Maths::intlog2(point - 1) - 3));
    }
}

inline double RuleWavelet::eval_linear(int point, double x) const{
    // Given a wavelet designated by point and a value x, evaluates the wavelet at x.

    // Standard Lifted Wavelets
    int l = Maths::intlog2(point - 1);
    int subindex = (point - 1) % (1 << l);
    double scale = std::pow(2,l-2);

    // Left Boundary
    if (subindex == 0){
        return linear_boundary_wavelet(scale * (x + 1.) - 1.);
    }
    // Right Boundary
    else if (subindex == (1 << l) - 1){
        return linear_boundary_wavelet(scale * (1. - x) - 1.);
    }
    double shift = 0.5 * (double (subindex - 1));
    return linear_central_wavelet(scale * (x + 1) - 1. - shift);
}

inline double RuleWavelet::linear_boundary_wavelet(double x) const{
    // Evaluates the first order boundary wavelet with support on [-1, 0].
    if (std::abs(x + 0.5) > 0.5) return 0.0;

    if ((x <= -0.75)){
        return 0.75 * (7.0 * x + 6.0);
    }else if (x <= -0.5){
        return -0.25 * (11.0 * x + 6.0);
    }else{
        return 0.25 * x;
    }
}

inline double RuleWavelet::linear_central_wavelet(double x) const {
    // Evaluates the first order central wavelet with support on [-1, .5].
    if (std::abs(x + 0.25) > 0.75) return 0.0;

    if ((x <= -0.5)){
        return -0.5 * x - 0.5;
    }else if (x >= 0.0){
        return 0.5 * x - 0.25;
    }else if (x <= -0.25){
        return 4.0 * x + 1.75;
    }else{
        return -4.0 * x - 0.25;
    }
}

inline int RuleWavelet::find_index(double x) const{
    // Finds an interval such that x_i <= x < x_i+1 using bisection search and returns i.
    if (x > 1. || x < -1.){
        return -1;
    }
    // Bisection search
    const double *xs = data[0].data();
    int low = 0;
    int high = num_data_points-1;
    while(high - low > 1){
        int test = (high + low)/2;
        if (x < xs[test]){
            high = test;
        }
        else{
            low = test;
        }
    }
    return low;
}

inline double RuleWavelet::interpolate(const double *y, double x) const{
    // For a given x value and dataset y, calculates the value of the interpolating
    // polynomial of given order going through the nearby points.
    constexpr int interpolation_order = 3;
    int idx = find_index(x);

    if (idx == -1){
        // Outside of table
        return 0.0;
    }

    // Neville's Algorithm
    if (idx < interpolation_order/2){
        idx = interpolation_order/2;
    }else if(num_data_points - idx - 1 < (interpolation_order+1)/2){
        idx = num_data_points - 1 - (interpolation_order+1)/2;
    }

    size_t start = (size_t)( idx - interpolation_order / 2 );
    double const *xs = &(data[0][start]);
    double const *ps = &(y[start]);

    double const dx0 = x - xs[0];
    double const dx1 = x - xs[1];
    double const dx2 = x - xs[2];
    double const dx3 = x - xs[3];

    double const *c = &(cachexs[4*start]);

    return   c[0] * ps[0] * dx1 * dx2 * dx3
           + c[1] * ps[1] * dx0 * dx2 * dx3
           + c[2] * ps[2] * dx0 * dx1 * dx3
           + c[3] * ps[3] * dx0 * dx1 * dx2;
}

} // namespace TasGrid
