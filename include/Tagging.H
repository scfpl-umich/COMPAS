// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef TAGGING_H
#define TAGGING_H

#include <ProblemICBC.H>
#include <AMReX_Array4.H>
#include <unordered_map>
#include <string>
#include <limits>
#include <AMReX_Print.H>

#include <unordered_map>
#include <string>

struct VarInfo {
    int index;        // component index in its MultiFab
    bool is_derived;  // false = state, true = derived
};



class VariableRegistry
{
public:
    static const VariableRegistry& instance()
    {
        static VariableRegistry registry;
        return registry;
    }

    bool has(const std::string& name) const
    {
        return m_map.find(name) != m_map.end();
    }

    const VarInfo& get(const std::string& name) const
    {
        auto it = m_map.find(name);
        if (it == m_map.end()) {
            amrex::Abort("VariableRegistry: unknown variable '" + name + "'");
        }
        return it->second;
    }

private:
    std::unordered_map<std::string, VarInfo> m_map;

    VariableRegistry()
    {
        // 1. Register STATE variables first (authoritative)
        for (int i = 0; i < vnamesC.size(); ++i) {
            m_map[vnamesC[i]] = VarInfo{i, false};
        }

        // 2. Register DERIVED variables only if not already present
        for (int i = 0; i < vnamesD.size(); ++i) {
            const std::string& name = vnamesD[i];

            if (m_map.find(name) != m_map.end()) {
                amrex::Print()
                    << "VariableRegistry: '" << name
                    << "' exists in both vnamesC and vnamesD; "
                    << "using state variable.\n";
                continue;
            }

            m_map[name] = VarInfo{i, true};
        }
    }
};


/** \brief Class that manages the refinement indicators for different variables
 * @todo Add the ability to manually refine in space
 */
struct TagVar{
public:
    /** Default constructor */
    TagVar(){};
    /** Default destructor */
    ~TagVar(){};

    /** String name of the variable */
    std::string name;
    /** Index of the variable in the MultiFab (not considering HO modes) */
    int index; 
    /** Cells are tagged on the levels below this one (default: every level) */
    int max_level = std::numeric_limits<int>::max();
    /** Indicates whether the variable needs to be calculated */
    bool derived = false;
    /** Tells ErrorEst whether or not to check for large gradients */
    bool do_grad = false;
    /** Tells ErrorEst whether or not to check for small values */
    bool do_value_less = false;
    /** Tells ErrorEst whether or not to check for large values */
    bool do_value_greater = false;
    /** Tells ErrorEst whether or not to check for ranges of values */
    bool do_value_range = false;

    /** Refinement thresholds in the cell-average gradient */
    amrex::Vector<amrex::Real> grad;
    /** Refinement thresholds for when the state is less than these values */
    amrex::Vector<amrex::Real> value_less;
    /** Refinement thresholds for when the state is greater than these values */
    amrex::Vector<amrex::Real> value_greater;
    /** Refinement thresholds for when the state is in between than these pairs of values */
    amrex::Vector<amrex::Real> value_range;

};



AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_value_greater (int i, int j, int k, int v,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real> const& state,
             amrex::Real value_greater, char tagval)
{
/** \brief Calculate the maximum magnitude in cell-averages for tagging
*/
    if (state(i,j,k,v) > value_greater){
        tag(i,j,k) = tagval;
    }
    
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_value_less (int i, int j, int k, int v,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real> const& state,
             amrex::Real value_less, char tagval)
{
/** \brief Calculate the minimum magnitude in cell-averages for tagging
*/
    if (state(i,j,k,v) < value_less){
        tag(i,j,k) = tagval;
    }

}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_value_range (int i, int j, int k, int v,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real> const& state,
             amrex::Real value_low,
             amrex::Real value_high,
             char tagval)
{
/** \brief Calculate the minimum magnitude in cell-averages for tagging
*/
    if (state(i,j,k,v) > value_low && state(i,j,k,v) < value_high){
        tag(i,j,k) = tagval;
    }

}



AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_grad (int i, int j, int k, int v,
            amrex::Array4<char> const& tag,
            amrex::Array4<amrex::Real> const& state,
            amrex::Real const grad_threshold,
            char const tagval) noexcept
{
/** \brief Calculate the maximum gradient in cell-averages for tagging
*/
    amrex::Real ax = amrex::Math::abs(state(i-1,j,k,v) - state(i,j,k,v));
    amrex::Real ay = amrex::Math::abs(state(i,j-1,k,v) - state(i,j,k,v));
    ax = amrex::max(ax, amrex::Math::abs(state(i,j,k,v) - state(i+1,j,k,v)));
    ay = amrex::max(ay, amrex::Math::abs(state(i,j,k,v) - state(i,j+1,k,v)));

    amrex::Real az;
    if (AMREX_SPACEDIM == 2) {
        az = 0.0;
    } else {
        az = amrex::Math::abs(state(i,j,k-1,v) - state(i,j,k,v));
        az = amrex::max(az, amrex::Math::abs(state(i,j,k,v) - state(i,j,k+1,v)));
    }

    if (amrex::max(ax,ay,az) >= grad_threshold) {
        tag(i,j,k) = tagval;
    }
}


#ifndef SPATIO_TEMPORAL_TAG_USER_FUNC_
#define SPATIO_TEMPORAL_TAG_USER_FUNC_
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void spatio_temporal_tag (int & TagVal, int lev, amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real t, Parm const& parm){};
#else
// The case defines its own spatio_temporal_tag, so ErrorEst calls it
#ifndef SPACE_TIME_REF
#define SPACE_TIME_REF
#endif
#endif

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
space_time (int i,
            int j,
            int k,
            int lev,
            amrex::Real t,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
            amrex::Array4<char> const& tag,
            Parm const& parm) noexcept
{
    int TagVal = -1;
    amrex::Real z = 0.0; 
#if (AMREX_SPACEDIM > 2)
    z = prob_lo[2] + (0.5+k) * dx[2];
#endif
    amrex::Real y = prob_lo[1] + (0.5+j) * dx[1];
    amrex::Real x = prob_lo[0] + (0.5+i) * dx[0];

    spatio_temporal_tag(TagVal,lev,x,y,z,t,parm);
    if (TagVal > -1){
        tag(i,j,k) = TagVal;
    }

}



#endif
