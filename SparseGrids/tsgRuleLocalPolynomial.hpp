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

#ifndef __TSG_RULE_LOCAL_POLYNOMIAL_HPP
#define __TSG_RULE_LOCAL_POLYNOMIAL_HPP

#include "tsgCoreOneDimensional.hpp"

namespace TasGrid{

#ifndef __TASMANIAN_DOXYGEN_SKIP

namespace RuleLocal {
    // effective rule type
    enum class erule {
        pwc, localp, semilocalp, localp0, localpb
    };

    inline erule getEffectiveRule(int order, TypeOneDRule rule) {
        if (order == 0) return erule::pwc;
        switch(rule) {
            case rule_semilocalp: return erule::semilocalp;
            case rule_localp0: return erule::localp0;
            case rule_localpb: return erule::localpb;
            default:
                return erule::localp;
        }
    }
    inline TypeOneDRule getRule(erule effective_rule) {
        switch(effective_rule) {
            case erule::pwc:
            case erule::localp:     return rule_localp;
            case erule::semilocalp: return rule_semilocalp;
            case erule::localp0:    return rule_localp0;
            default: //case erule::localpb:
                return rule_localpb;
        };
    }

    template<erule effective_rule>
    int getNumPoints(int level) {
        switch(effective_rule) {
            case erule::pwc: {
                int n = 1;
                while (level-- > 0) n *= 3;
                return n;
            }
            case erule::localp:
            case erule::semilocalp:
                return (level == 0) ? 1 : ((1 << level) + 1);
            case erule::localp0:
                return (1 << (level+1)) -1;
            default: // case erule::localpb:
                return ((1 << level) + 1);
        };
    }

    template<erule effective_rule>
    int getMaxNumKids() { return (effective_rule == erule::pwc) ? 4 : 2; }
    template<erule effrule>
    int getMaxNumParents() {
        return ((effrule == erule::pwc or effrule == erule::semilocalp or effrule == erule::localpb) ? 2 : 1);
    }

    template<erule effective_rule>
    int getParent(int point) {
        switch(effective_rule) {
            case erule::pwc:
                return (point == 0) ? -1 : point / 3;
            case erule::localp:
            case erule::semilocalp: {
                int dad = (point + 1) / 2;
                if (point < 4) dad--;
                return dad;
            }
            case erule::localp0:
                return (point == 0) ? -1 : (point - 1) / 2;
            default: // case erule::localpb:
                return (point < 2) ? -1 : ((point + 1) / 2);
        };
    }

    template<erule effective_rule>
    int getStepParent(int point) {
        if (effective_rule == erule::pwc){
            int i3l3 = Maths::int3log3(point);
            if (point == i3l3/3) return -1;
            if (point == i3l3-1) return -1;
            int mod3 = point % 3;
            int mod2 = point % 2;
            if (mod3 == 2 and mod2 == 0) return point / 3 + 1;
            if (mod3 == 0 and mod2 == 1) return point / 3 - 1;
            return -1;
        }else{
            if (effective_rule == erule::semilocalp){
                switch(point) {
                    case 3: return 2;
                    case 4: return 1;
                    default:
                        return -1;
                };
            }else if (effective_rule == erule::localpb){
                return (point == 2) ? 0 : -1;
            }
            return -1;
        }
    }

    template<erule effective_rule>
    int getKid(int point, int kid_number) {
        switch(effective_rule) {
            case erule::pwc:
                if (point == 0) return (kid_number == 0) ? 1 : (kid_number==1) ? 2 : -1;
                if (kid_number == 3){
                    int i3l3 = Maths::int3log3(point);
                    if (point == i3l3/3) return -1;
                    if (point == i3l3-1) return -1;
                    return (point % 2 == 0) ? 3*point + 3 : 3*point - 1;
                }
                return 3*point + kid_number;
            case erule::localp:
            case erule::semilocalp:
                if (kid_number == 0){
                    switch(point) {
                        case 0: return 1;
                        case 1: return 3;
                        case 2: return 4;
                        default:
                            return 2 * point - 1;
                    };
                } else {
                    switch(point) {
                        case 0: return 2;
                        case 1:
                        case 2: return -1;
                        default:
                            return 2 * point;
                    };
                }
            case erule::localp0:
                return 2 * point + ((kid_number == 0) ? 1 : 2);
            default: // case erule::localpb:
                switch(point) {
                    case 0:
                    case 1:
                        return (kid_number == 0) ? 2 : -1;
                    default:
                        return 2*point - ((kid_number == 0) ? 1 : 0);
                };
        };
    }

    template<erule effective_rule>
    double getNode(int point) {
        switch(effective_rule) {
            case erule::pwc:
                return -2.0 + (1.0 / ((double) Maths::int3log3(point))) * (3*point + 2 - point % 2);
            case erule::localp:
            case erule::semilocalp:
                switch(point) {
                    case 0: return  0.0;
                    case 1: return -1.0;
                    case 2: return  1.0;
                    default:
                        return ((double)(2*point - 1)) / ((double) Maths::int2log2(point - 1)) - 3.0;
                };
            case erule::localp0:
                return ((double)(2*point + 3) ) / ((double) Maths::int2log2(point + 1) ) - 3.0;
            default: // case erule::localpb:
                switch(point) {
                    case 0: return -1.0;
                    case 1: return  1.0;
                    case 2: return  0.0;
                    default:
                        return ((double)(2*point - 1)) / ((double) Maths::int2log2(point - 1)) - 3.0;
                };
        };
    }

    template<erule effective_rule>
    int getLevel(int point) {
        switch(effective_rule) {
            case erule::pwc: {
                int level = 0;
                while(point >= 1){ point /= 3; level += 1; }
                return level;
            }
            case erule::localp:
            case erule::semilocalp:
                return (point == 0) ? 0 : (point == 1) ? 1 : (Maths::intlog2(point - 1) + 1);
            case erule::localp0:
                return Maths::intlog2(point + 1);
            default: // case erule::localpb:
                return (point <= 1) ? 0 : (Maths::intlog2(point - 1) + 1);
        };
    }

    template<erule effective_rule>
    double getSupport(int point) {
        switch(effective_rule) {
            case erule::pwc:
                return 1.0 / (double) Maths::int3log3(point);
            case erule::localp:
            case erule::semilocalp:
                return (point == 0) ? 1.0 : 1.0 / ((double) Maths::int2log2(point - 1));
            case erule::localp0:
                return 1.0 / ((double) Maths::int2log2(point + 1));
            default: // case erule::localpb:
                return (point <= 1) ? 2.0 : 1.0 / ((double) Maths::int2log2(point - 1));
        };
    }

    template<erule effective_rule>
    double scaleDiffX(int point) {
        switch(effective_rule) {
            case erule::pwc:
                return 0.0;
            case erule::localp:
                return (point <= 2) ? 1.0 : static_cast<double>(Maths::int2log2(point - 1));
            case erule::semilocalp:
                return static_cast<double>(Maths::int2log2(point - 1));
            case erule::localp0:
                return (point == 0) ? 1.0 : static_cast<double>(Maths::int2log2(point + 1));
            default: // case erule::localpb:
                switch(point) {
                    case 0:
                    case 1: return 0.5;
                    case 2: return 1.0;
                    default:
                        return static_cast<double>(Maths::int2log2(point - 1));
                };
        };
    }

    template<erule effective_rule>
    double scaleX(int point, double x) {
        switch(effective_rule) {
            case erule::pwc:
                return 0.0; // should never be called
            case erule::localp:
                switch(point) {
                    case 0: return x;
                    case 1: return (x + 1.0);
                    case 2: return (x - 1.0);
                    default:
                        return ((double) Maths::int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
                };
            case erule::semilocalp:
                return ((double) Maths::int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
            case erule::localp0:
                return ((double) Maths::int2log2(point + 1) * (x + 3.0) - 3.0 - (double) (2*point));
            default: // case erule::localpb:
                switch(point) {
                    case 0: return (x + 1.0) / 2.0;
                    case 1: return (x - 1.0) / 2.0;
                    case 2: return x;
                    default:
                        return ((double) Maths::int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
                };
        };
    }

    template<erule effective_rule>
    double evalPWQuadratic(int point, double x) {
        if (effective_rule == erule::localp){
            switch(point) {
                case 1: return 1.0 - x;
                case 2: return 1.0 + x;
                default:
                    return (1.0 - x) * (1.0 + x);
            };
        }else if (effective_rule == erule::localpb){
            switch(point) {
                case 0: return 1.0 - x;
                case 1: return 1.0 + x;
                default:
                    return (1.0 - x) * (1.0 + x);
            };
        }
        return (1.0 - x) * (1.0 + x);
    }
    template<erule effective_rule>
    double evalPWCubic(int point, double x) {
        if (effective_rule == erule::localp){
            switch(point) {
                case 0: return 1.0;
                case 1: return 1.0 - x;
                case 2: return 1.0 + x;
                case 3:
                case 4: return (1.0 - x) * (1.0 + x);
                default:
                    return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
            };
        }else if (effective_rule == erule::localp0){
            if (point == 0) return (1.0 - x) * (1.0 + x);
        }else if (effective_rule == erule::localpb){
            switch(point) {
                case 0: return 1.0 - x;
                case 1: return 1.0 + x;
                case 2: return (1.0 - x) * (1.0 + x);
                default:
                    return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
            };
        }
        return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
    }

    template<erule effective_rule>
    double evalPWPower(int max_order, int point, double x) {
        // use the cubic implementation until we have enough points
        if (effective_rule == erule::localp)     if (point <= 8) return evalPWCubic<effective_rule>(point, x); // if order is cubic or less, use the hard-coded functions
        if (effective_rule == erule::semilocalp) if (point <= 4) return evalPWCubic<effective_rule>(point, x);
        if (effective_rule == erule::localpb)    if (point <= 4) return evalPWCubic<effective_rule>(point, x);
        if (effective_rule == erule::localp0)    if (point <= 2) return evalPWCubic<effective_rule>(point, x);
        int level = getLevel<effective_rule>(point);
        int most_turns = 1;
        double value = (1.0 - x)*(1.0 + x);
        double phantom_distance = 1.0;
        int max_ancestors = [&]()->int{ // maximum number of ancestors to consider, first set the possible max then constrain by max_order
            switch(effective_rule) {
                case erule::pwc:        return 0; // should never happen
                case erule::localp:     return max_ancestors = level-2;
                case erule::semilocalp: return max_ancestors = level-1;
                case erule::localpb:    return max_ancestors = level-1;
                default:                return max_ancestors = level;
            };
        }();
        if (max_order > 0) max_ancestors = std::min(max_ancestors, max_order - 2); // use the minimum of the available ancestors or the order restriction

        for(int j=0; j < max_ancestors; j++){
            // Lagrange polynomial needs to be constructed using normalized x (basis support (-1, 1)).
            // The support of the basis is equal to 2, thus we use units of "half-support"
            // The first two nodes are used in the initialization of value, those are the nearest ancestors (in spacial distance, not hierarchy level).
            // The other nodes are "phantoms" that lay strictly outside of the support [-1, 1] (i.e., not on the edge)
            // The walking distance more than doubles and is an odd number (due to the half-support)
            // The walk through the ancestors can take a left or right turn at each step,
            //   most_turns is the total number of turns possible for the current ancestor, turns (or most_turns - 1 - turns) is the actual number of turns
            // Every time we turn, we backtrack and we lose 2 units, the phantom node is at maximum distance minus 2 time the number of turns (to the left or right)
            most_turns *= 2;
            phantom_distance = 2.0 * phantom_distance + 1.0;
            int turns = (effective_rule == erule::localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            double node = (turns < most_turns / 2) ? (phantom_distance - 2.0 * ((double) turns)) : (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
            value *= - ( x - node ) / node;
        }
        return value;
    }

    template<erule effective_rule>
    double evalSupport(int max_order, int point, double x, bool &isSupported) {
        switch(effective_rule) {
            case erule::pwc: {
                    double distance = std::abs(x - getNode<effective_rule>(point));
                    double support  = getSupport<effective_rule>(point);
                    isSupported = (distance <= 2.0 * support);
                    return (distance <= support) ? 1.0 : 0.0;
                }
            case erule::localp:
                isSupported = true;
                if (point == 0) {
                    return 1.0;
                } else {
                    double xn = scaleX<effective_rule>(point, x);
                    if (std::abs(xn) <= 1.0) {
                        switch(max_order) {
                            case 1: return 1.0 - std::abs(xn);
                            case 2: return evalPWQuadratic<effective_rule>(point, xn);
                            case 3: return evalPWCubic<effective_rule>(point, xn);
                            default:
                                return evalPWPower<effective_rule>(max_order, point, xn);
                        };
                    } else {
                        isSupported = false;
                        return 0.0;
                    }
                }
            case erule::semilocalp:
                isSupported = true;
                switch(point) {
                    case 0: return 1.0;
                    case 1: return 0.5 * x * (x - 1.0);
                    case 2: return 0.5 * x * (x + 1.0);
                    default: {
                        double xn = scaleX<effective_rule>(point, x);
                        if (std::abs(xn) <= 1.0) {
                            switch(max_order) {
                                case 1: return 1.0 - std::abs(xn);
                                case 2: return evalPWQuadratic<effective_rule>(point, xn);
                                case 3: return evalPWCubic<effective_rule>(point, xn);
                                default:
                                    return evalPWPower<effective_rule>(max_order, point, xn);
                            };
                        } else {
                            isSupported = false;
                            return 0.0;
                        }
                    }
                };
            case erule::localp0:
            default: { // case erule::localpb:
                    double xn = scaleX<effective_rule>(point, x);
                    if (std::abs(xn) <= 1.0) {
                        isSupported = true;
                        switch(max_order) {
                            case 1: return 1.0 - std::abs(xn);
                            case 2: return evalPWQuadratic<effective_rule>(point, xn);
                            case 3: return evalPWCubic<effective_rule>(point, xn);
                            default:
                                return evalPWPower<effective_rule>(max_order, point, xn);
                        };
                    } else {
                        isSupported = false;
                        return 0.0;
                    }
                }
        };
    }

    template<erule effective_rule>
    double evalRaw(int max_order, int point, double x) {
        switch(effective_rule) {
            case erule::pwc:
                return (std::abs(x - getNode<effective_rule>(point)) <= getSupport<effective_rule>(point)) ? 1.0 : 0.0;
            case erule::localp:
                if (point == 0) {
                    return 1.0;
                } else {
                    double xn = scaleX<effective_rule>(point, x);
                    if (std::abs(xn) <= 1.0) {
                        switch(max_order) {
                            case 1: return 1.0 - std::abs(xn);
                            case 2: return evalPWQuadratic<effective_rule>(point, xn);
                            case 3: return evalPWCubic<effective_rule>(point, xn);
                            default:
                                return evalPWPower<effective_rule>(max_order, point, xn);
                        };
                    } else {
                        return 0.0;
                    }
                }
            case erule::semilocalp:
                switch(point) {
                    case 0: return 1.0;
                    case 1: return 0.5 * x * (x - 1.0);
                    case 2: return 0.5 * x * (x + 1.0);
                    default: {
                        double xn = scaleX<effective_rule>(point, x);
                        if (std::abs(xn) <= 1.0) {
                            switch(max_order) {
                                case 1: return 1.0 - std::abs(xn);
                                case 2: return evalPWQuadratic<effective_rule>(point, xn);
                                case 3: return evalPWCubic<effective_rule>(point, xn);
                                default:
                                    return evalPWPower<effective_rule>(max_order, point, xn);
                            };
                        } else {
                            return 0.0;
                        }
                    }
                };
            case erule::localp0:
            default: { // case erule::localpb:
                    double xn = scaleX<effective_rule>(point, x);
                    if (std::abs(xn) <= 1.0) {
                        switch(max_order) {
                            case 1: return 1.0 - std::abs(xn);
                            case 2: return evalPWQuadratic<effective_rule>(point, xn);
                            case 3: return evalPWCubic<effective_rule>(point, xn);
                            default:
                                return evalPWPower<effective_rule>(max_order, point, xn);
                        };
                    } else {
                        return 0.0;
                    }
                }
        };
    }

    template<erule effective_rule>
    double diffPWQuadratic(int point, double x) {
        if (effective_rule == erule::localp) {
            switch(point) {
                case 1: return -1.0;
                case 2: return  1.0;
                default:
                    return -2.0 * x;
            };
        } else if (effective_rule == erule::localpb) {
            switch(point) {
                case 0: return -1.0;
                case 1: return  1.0;
                default:
                    return -2.0 * x;
            };
        }
        return -2.0 * x;
    }
    template<erule effective_rule>
    double diffPWCubic(int point, double x) {
        if (effective_rule == erule::localp) {
            switch(point) {
                case 0: return  0.0;
                case 1: return -1.0;
                case 2: return  1.0;
                case 3:
                case 4: return  -2.0 * x;
                default:
                    return (point % 2 == 0) ? 1.0 / 3.0 - x * (x + 2.0) : -1.0 / 3.0 + x * (x - 2.0);
            };
        } else if (effective_rule == erule::localpb) {
            switch(point) {
                case 0: return -1.0;
                case 1: return  1.0;
                case 2: return -2.0 * x;
                default:
                    return (point % 2 == 0) ? 1.0 / 3.0 - x * (x + 2.0) : -1.0 / 3.0 + x * (x - 2.0);
            };
        } else if (effective_rule == erule::localp0) {
            return (point == 0) ? -2.0 * x : ((point % 2 == 0) ? 1.0 / 3.0 - x * (x + 2.0) : -1.0 / 3.0 + x * (x - 2.0));
        }
        return (point % 2 == 0) ? 1.0 / 3.0 - x * (x + 2.0) : -1.0 / 3.0 + x * (x - 2.0);
    }

    template<erule effrule>
    double diffPWPower(int max_order, int point, double x) {
        // use the cubic implementation until we have enough points
        if (effrule == erule::localp     and point <= 8) return diffPWCubic<effrule>(point, x);
        if (effrule == erule::semilocalp and point <= 4) return diffPWCubic<effrule>(point, x);
        if (effrule == erule::localpb    and point <= 4) return diffPWCubic<effrule>(point, x);
        if (effrule == erule::localp0    and point <= 2) return diffPWCubic<effrule>(point, x);
        int level = getLevel<effrule>(point);
        int max_ancestors = [&]()->int {
            switch(effrule) {
                case erule::pwc:        return 0;
                case erule::localp:     return level - 2;
                case erule::semilocalp:
                case erule::localpb:    return level - 1;
                default:                return level;
            };
        }();

        if (max_order > 0) max_ancestors = std::min(max_ancestors, max_order - 2);

        // This lambda captures most_turns and phantom_distance by reference, uses those as internal state variables, and on each
        // call it returns the next normalized ancestor node.
        int most_turns = 1;
        double phantom_distance = 1.0;
        auto update_and_get_next_node = [&]() {
            most_turns *= 2;
            phantom_distance = 2.0 * phantom_distance + 1.0;
            int turns = (effrule == erule::localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            return (turns < most_turns / 2) ?
                    (phantom_distance - 2.0 * ((double) turns)) :
                    (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
        };

        // This lambda is the inverse transform of update_and_get_next_node() above, and returns the previous descendant node.
        auto rollback_and_get_prev_node = [&]() {
            most_turns /= 2;
            phantom_distance = 0.5 * (phantom_distance - 1.0);
            int turns = (effrule == erule::localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            return (turns < most_turns / 2) ?
                    (phantom_distance - 2.0 * ((double) turns)) :
                    (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
        };

        // Does not include the additional factor (1-x) * (1+x) and the Lagrange coefficient.
        std::vector<double> left_prods(max_ancestors);
        left_prods[0] = 1.0;
        double node = update_and_get_next_node();
        double coeff = 1.0 / (-node);
        for(int j=1; j<max_ancestors; j++) {
            left_prods[j] = left_prods[j-1] * (x - node);
            node = update_and_get_next_node();
            coeff *= 1.0 / (-node);
        }
        double right_prod = 1.0;
        double derivative = left_prods[max_ancestors-1];
        for (int j=max_ancestors-2; j>=0; j--) {
            right_prod *= x - node;
            derivative += right_prod * left_prods[j];
            node = rollback_and_get_prev_node();
        }

        // Adjust for the additional factor (1-x) * (1+x) and the Lagrange coefficient.
        derivative = derivative * (1.0 - x) * (1.0 + x) + right_prod * (x - node) * (-2.0) * x;
        derivative *= coeff;

        return derivative;
    }

    template<erule effrule>
    double diffRaw(int max_order, int point, double x) {
        switch(effrule) {
            case erule::pwc: return 0.0;
            case erule::localp: {
                if (point == 0) return 0.0;
                double xn = scaleX<effrule>(point, x);
                double an = scaleDiffX<effrule>(point);
                switch(max_order) {
                    case 1: return ((xn >= 0 ? -1.0 : 1.0) * an);
                    case 2: return an * diffPWQuadratic<effrule>(point, xn);
                    case 3: return an * diffPWCubic<effrule>(point, xn);
                    default:
                        return an * diffPWPower<effrule>(point, xn);
                };
            }
            case erule::semilocalp:
                switch(point) {
                    case 0: return 0.0;
                    case 1: return x - 0.5;
                    case 2: return x + 0.5;
                    default: {
                        double xn = scaleX<effrule>(point, x);
                        double an = scaleDiffX<effrule>(point);
                        switch(max_order) {
                            case 2: return an * diffPWQuadratic<effrule>(point, xn);
                            case 3: return an * diffPWCubic<effrule>(point, xn);
                            default:
                                return an * diffPWPower<effrule>(point, xn);
                        };
                    }
                };
            case erule::localp0: {
                double xn = scaleX<effrule>(point, x);
                double an = scaleDiffX<effrule>(point);
                switch(max_order) {
                    case 1: return (x == 1.0 and point == 0) ? -1.0 : ((xn >= 0 ? -1.0 : 1.0) * an);
                    case 2: return an * diffPWQuadratic<effrule>(point, xn);
                    case 3: return an * diffPWCubic<effrule>(point, xn);
                    default:
                        return an * diffPWPower<effrule>(point, xn);
                };
            }
            default: { // case erule::localpb:
                double xn = scaleX<effrule>(point, x);
                double an = scaleDiffX<effrule>(point);
                switch(max_order) {
                    case 1: return ((xn >= 0 ? -1.0 : 1.0) * an);
                    case 2: return an * diffPWQuadratic<effrule>(point, xn);
                    case 3: return an * diffPWCubic<effrule>(point, xn);
                    default:
                        return an * diffPWPower<effrule>(point, xn);
                };
            }
        };
    }

    template<erule effrule>
    double diffSupport(int max_order, int point, double x, bool &isSupported) {
        switch(effrule) {
            case erule::pwc:
                isSupported = false;
                return 0.0;
            case erule::localp: {
                if (point == 0) { isSupported = true; return 0.0; }
                double xn = scaleX<effrule>(point, x);
                isSupported = (-1.0 <= xn and xn < 1.0) or (x == 1.0 and xn == 1.0);
                if (isSupported) {
                    double an = scaleDiffX<effrule>(point);
                    switch(max_order) {
                        case 1: return (x == 1.0 and point == 2) ? an : (((xn >= 0 ? -1.0 : 1.0) * an));
                        case 2: return an * diffPWQuadratic<effrule>(point, xn);
                        case 3: return an * diffPWCubic<effrule>(point, xn);
                        default:
                            return an * diffPWPower<effrule>(max_order, point, xn);
                    };
                } else {
                    return 0.0;
                }
            }
            case erule::semilocalp:
                switch(point) {
                    case 0: { isSupported = true; return 0.0; }
                    case 1: { isSupported = true; return x - 0.5; }
                    case 2: { isSupported = true; return x + 0.5; }
                    default: {
                        double xn = scaleX<effrule>(point, x);
                        isSupported = (-1.0 <= xn and xn < 1.0) or (x == 1.0 and xn == 1.0);
                        if (isSupported) {
                            double an = scaleDiffX<effrule>(point);
                            switch(max_order) {
                                case 2: return an * diffPWQuadratic<effrule>(point, xn);
                                case 3: return an * diffPWCubic<effrule>(point, xn);
                                default:
                                    return an * diffPWPower<effrule>(max_order, point, xn);
                            };
                        } else {
                            return 0.0;
                        }
                    }
                };
            case erule::localp0: {
                double xn = scaleX<effrule>(point, x);
                isSupported = (-1.0 <= xn and xn < 1.0) or (x == 1.0 and xn == 1.0);
                if (isSupported) {
                    double an = scaleDiffX<effrule>(point);
                    switch(max_order) {
                        case 1: return (x == 1.0 and point == 0) ? -1.0 : ((xn >= 0 ? -1.0 : 1.0) * an);
                        case 2: return an * diffPWQuadratic<effrule>(point, xn);
                        case 3: return an * diffPWCubic<effrule>(point, xn);
                        default:
                            return an * diffPWPower<effrule>(max_order, point, xn);
                    };
                } else {
                    return 0.0;
                }
            }
            default: { // case erule::localpb:
                double xn = scaleX<effrule>(point, x);
                isSupported = (-1.0 <= xn and xn < 1.0) or (x == 1.0 and xn == 1.0);
                if (isSupported) {
                    double an = scaleDiffX<effrule>(point);
                    switch(max_order) {
                        case 1: return ((xn >= 0 ? -1.0 : 1.0) * an);
                        case 2: return an * diffPWQuadratic<effrule>(point, xn);
                        case 3: return an * diffPWCubic<effrule>(point, xn);
                        default:
                            return an * diffPWPower<effrule>(max_order, point, xn);
                    };
                } else {
                    return 0.0;
                }
            }
        };
    }

    template<erule effrule>
    double getArea(int max_order, int point, std::vector<double> const &w, std::vector<double> const &x) {
        switch(effrule) {
            case erule::pwc:
                return 2.0 * getSupport<effrule>(point);
            case erule::localp:
                switch(point) {
                    case 0: return 2.0;
                    case 1:
                    case 2: return 0.5;
                    default:
                        switch(max_order) {
                            case 1: return getSupport<effrule>(point);
                            case 2:
                            case 3: return (4.0/3.0) * getSupport<effrule>(point);
                            default:
                                if (point <= 8) return (4.0/3.0) * getSupport<effrule>(point);
                                break;
                        };
                        break;
                };
                break;
            case erule::semilocalp:
                switch(point) {
                    case 0: return 2.0;
                    case 1:
                    case 2: return 1.0/3.0;
                    default:
                        switch(max_order) {
                            case 2:
                            case 3: return (4.0/3.0) * getSupport<effrule>(point);
                            default:
                                if (point <= 4) return (4.0/3.0) * getSupport<effrule>(point);
                                break;
                        };
                        break;
                };
                break;
            case erule::localp0:
                switch(max_order) {
                    case 1: return getSupport<effrule>(point);
                    case 2:
                    case 3: return (4.0/3.0) * getSupport<effrule>(point);
                    default:
                        if (point <= 2) return (4.0/3.0) * getSupport<effrule>(point);
                        break;
                };
                break;
            default: // case erule::localpb:
                if (point <= 1) {
                    return 1.0;
                } else {
                    switch(max_order) {
                        case 1: return getSupport<effrule>(point);
                        case 2:
                        case 3: return (4.0/3.0) * getSupport<effrule>(point);
                        default:
                            if (point <= 4) return (4.0/3.0) * getSupport<effrule>(point);
                            break;
                    };
                }
                break;
        };
        double sum = 0.0;
        for(size_t i=0; i<w.size(); i++) sum += w[i] * evalPWPower<effrule>(max_order, point, x[i]);
        return sum * getSupport<effrule>(point);
    }

    template<erule effrule>
    void van_matrix(int max_order, int num_rows, std::vector<int> &pntr, std::vector<int> &indx, std::vector<double> &vals) {
        int max_level = getLevel<effrule>(num_rows);

        if (effrule == erule::pwc) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                case 1:
                    indx = {0, 0, 1};
                    vals = {1.0, 1.0, 1.0};
                    pntr = {0, 1, 3};
                    break;
                case 2:
                    indx = {0, 0, 1, 0, 2};
                    vals = {1.0, 1.0, 1.0, 1.0, 1.0};
                    pntr = {0, 1, 3, 5};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    for(auto i : std::array<int, 5>{0, 0, 1, 0, 2})
                        indx.push_back(i);
                    for(size_t i=0; i<5; i++)
                        vals.push_back(1.0);
                    pntr[1] = 1;
                    pntr[2] = 3;
                    {
                        std::vector<int> ancestors;
                        ancestors.reserve(max_level);
                        for(int r=3; r<num_rows; r++) {
                            pntr[r] = static_cast<int>(indx.size());
                            ancestors.clear();
                            int kid = r;
                            int dad = kid / 3;
                            while(dad != 0) {
                                if (dad % 2 == 0) {
                                    if (kid != 3 * dad)
                                        ancestors.push_back(dad);
                                } else {
                                    if (kid != 3 * dad + 2)
                                        ancestors.push_back(dad);
                                }
                                kid = dad;
                                dad = kid / 3;
                            }
                            indx.push_back(0);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            for(size_t i=0; i<ancestors.size() + 2; i++)
                                vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        } else if (effrule == erule::localp) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    indx.push_back(0);
                    vals.push_back(1.0);
                    {
                        std::vector<int> ancestors;
                        std::vector<double> ancestors_vals;
                        ancestors.reserve(max_level);
                        ancestors_vals.reserve(max_level);
                        for(int r=1; r<num_rows; r++) {
                            double x = getNode<effrule>(r);
                            pntr[r] = static_cast<int>(indx.size());
                            ancestors.clear();
                            ancestors_vals.clear();
                            int kid = r;
                            int dad = (kid + 1) / 2;
                            if (kid < 4) dad--;
                            while(dad != 0) {
                                ancestors.push_back(dad);
                                ancestors_vals.push_back( evalRaw<effrule>(max_order, dad, x) );
                                kid = dad;
                                dad = (kid + 1) / 2;
                                if (kid < 4) dad--;
                            }
                            indx.push_back(0);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            vals.push_back(1.0);
                            vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                            vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        } else if (effrule == erule::semilocalp) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                case 1:
                    indx = {0, 0, 1};
                    vals = {1.0, 1.0, 1.0};
                    pntr = {0, 1, 3};
                    break;
                case 2:
                    indx = {0, 0, 1, 0, 2};
                    vals = {1.0, 1.0, 1.0, 1.0, 1.0};
                    pntr = {0, 1, 3, 5};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    for(auto i : std::array<int, 5>{0, 0, 1, 0, 2})
                        indx.push_back(i);
                    for(int i=0; i<5; i++)
                        vals.push_back(1.0);
                    pntr[1] = 1;
                    pntr[2] = 3;
                    {
                        std::vector<int> ancestors;
                        std::vector<double> ancestors_vals;
                        ancestors.reserve(max_level);
                        ancestors_vals.reserve(max_level);
                        for(int r=3; r<num_rows; r++) {
                            pntr[r] = static_cast<int>(indx.size());
                            double x = getNode<effrule>(r);
                            ancestors.clear();
                            ancestors_vals.clear();
                            int kid = r;
                            int dad = (kid + 1) / 2;
                            while(dad > 2) {
                                ancestors.push_back(dad);
                                ancestors_vals.push_back( evalRaw<effrule>(max_order, dad, x) );
                                kid = dad;
                                dad = (kid + 1) / 2;
                            }
                            indx.push_back(0);
                            indx.push_back(1);
                            indx.push_back(2);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            vals.push_back(1.0);
                            vals.push_back( evalRaw<effrule>(max_order, 1, x) );
                            vals.push_back( evalRaw<effrule>(max_order, 2, x) );
                            vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                            vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        } else if (effrule == erule::localp0) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    indx.push_back(0);
                    vals.push_back(1.0);
                    {
                        std::vector<int> ancestors;
                        std::vector<double> ancestors_vals;
                        ancestors.reserve(max_level);
                        ancestors_vals.reserve(max_level);
                        for(int r=1; r<num_rows; r++) {
                            double x = getNode<effrule>(r);
                            pntr[r] = static_cast<int>(indx.size());
                            ancestors.clear();
                            ancestors_vals.clear();
                            int kid = r;
                            int dad = (kid - 1) / 2;
                            while(dad != 0) {
                                ancestors.push_back(dad);
                                ancestors_vals.push_back( evalRaw<effrule>(max_order, dad, x) );
                                kid = dad;
                                dad = (kid - 1) / 2;
                            }
                            indx.push_back(0);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            vals.push_back( evalRaw<effrule>(max_order, 0, x) );
                            vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                            vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        } else { // if (effrule == erule::localpb) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                case 1:
                    indx = {0, 1};
                    vals = {1.0, 1.0};
                    pntr = {0, 1, 2};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    indx.push_back(0);
                    indx.push_back(1);
                    vals.push_back(1.0);
                    vals.push_back(1.0);
                    pntr[1] = 1;
                    {
                        std::vector<int> ancestors;
                        std::vector<double> ancestors_vals;
                        ancestors.reserve(max_level);
                        ancestors_vals.reserve(max_level);
                        for(int r=2; r<num_rows; r++) {
                            pntr[r] = static_cast<int>(indx.size());
                            double x = getNode<effrule>(r);
                            ancestors.clear();
                            ancestors_vals.clear();
                            int kid = r;
                            int dad = (kid + 1) / 2;
                            while(dad > 1) {
                                ancestors.push_back(dad);
                                ancestors_vals.push_back( evalRaw<effrule>(max_order, dad, x) );
                                kid = dad;
                                dad = (kid + 1) / 2;
                            }
                            indx.push_back(0);
                            indx.push_back(1);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            vals.push_back( evalRaw<effrule>(max_order, 0, x) );
                            vals.push_back( evalRaw<effrule>(max_order, 1, x) );
                            vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                            vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        }
    }

} // namespace RuleLocal

#endif // __TASMANIAN_DOXYGEN_SKIP

}

#endif
